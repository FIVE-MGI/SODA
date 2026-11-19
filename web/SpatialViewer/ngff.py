r"""OME-Zarr (OME-NGFF) images and NGFF coordinate transformations, read with zarr 2 (zarr v2) or zarrv3.py (zarr v3).

Used by `prepare.py omezarr` (an OME-Zarr image on its own) and `prepare.py spatialdata`
(the images/ of a SpatialData store, which are OME-Zarr images too) in stdformats.py.

    open_image(path)            -> NgffImage: the multiscale image at `path` (a .zarr folder, a
                                   bioformats2raw layout's series, or one image group of a SpatialData store)
    NgffImage.source(...)       -> OmeZarrSource, a pyramid.Source over level 0 (one channel, one band
                                   of rows at a time; never the whole image)
    element_transforms(attrs)   -> {coordinate system: 3x3 affine} of a SpatialData element
    to_affine(t, axes)          -> one NGFF / SpatialData transformation as a 3x3 affine on (x, y)

What is read (NGFF 0.1 to 0.4 in the zarr v2 layout, read with zarr 2; NGFF 0.5 in the zarr v3 layout,
read with zarrv3.py, its metadata under attributes["ome"]):
    multiscales[0]: axes (names, types, units; 0.3 strings; none in 0.1 / 0.2 = t, c, z, y, x),
                    datasets (level 0 = the first, its scale gives the pixel size when the x / y axes carry a
                    unit), the multiscale-level coordinateTransformations
    omero.channels: label (channel name) and color (hex); window is not used (SODA sets its own)
    bioformats2raw.layout: the series folders 0, 1, ...; --series picks one (default 0)
Refused, with the reason: a zarr v3 array outside the subset zarrv3.py reads (the codec or data type is
named), HCS plates / wells, an image with more than one time point, an axis of an unknown type with more
than one entry, and labels (segmentation masks) given as the image.

Transformations (to_affine): identity, scale, translation, affine and sequence (of those) on the
x / y axes. A transformation that mixes x / y with another axis (c, z), mapAxis, rotation given
by name, and anything else is refused by name: SODA does not guess a frame.
"""
from __future__ import annotations

import json
import threading
from pathlib import Path

import numpy as np

from pyramid import Source

UNIT_UM = {"micrometer": 1.0, "micron": 1.0, "um": 1.0, "µm": 1.0, "nanometer": 1e-3, "nm": 1e-3,
           "millimeter": 1e3, "mm": 1e3, "centimeter": 1e4, "cm": 1e4, "meter": 1e6, "m": 1e6,
           "angstrom": 1e-4, "picometer": 1e-6}
_RGB_NAMES = ({"r", "g", "b"}, {"red", "green", "blue"}, {"0", "1", "2"})
_RGB_COLORS = {"ff0000", "00ff00", "0000ff"}


class NgffError(SystemExit):
    """A store SODA cannot read, with the reason (a SystemExit: prepare.py stops with the message)."""


# ---------------------------------------------------------------------------------------
# zarr store helpers (v2 metadata files are plain JSON)
# ---------------------------------------------------------------------------------------
def read_attrs(folder) -> dict:
    """The attributes of a zarr group / array folder ({} when there are none): the .zattrs of zarr v2, the
    zarr.json "attributes" of zarr v3, where NGFF 0.5's "ome" block is lifted to the top (so multiscales,
    omero, bioformats2raw.layout, plate, image-label read as in 0.4; its version becomes "version")."""
    p = Path(folder) / ".zattrs"
    if not p.exists():
        if (Path(folder) / "zarr.json").exists():
            import zarrv3
            a = zarrv3.attrs_of(folder)
            ome = a.get("ome")
            if isinstance(ome, dict):
                a = {**{k: v for k, v in a.items() if k != "ome"}, **ome}
            return a
        return {}
    try:
        return json.loads(p.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}


def zarr_version(folder) -> int | None:
    """2 (.zgroup / .zarray), 3 (zarr.json) or None (not a zarr node)."""
    f = Path(folder)
    if (f / ".zgroup").exists() or (f / ".zarray").exists():
        return 2
    if (f / "zarr.json").exists():
        return 3
    return None


def open_array(path):
    """A zarr array folder, either format: zarr 2's Array for v2, zarrv3.Array for v3 (the same reading API:
    shape, dtype, ndim, chunks, a[...])."""
    if zarr_version(path) == 3:
        import zarrv3
        return zarrv3.open_array(path)
    import zarr
    return zarr.open_array(str(path), mode="r")


def open_group(path):
    """A zarr group folder, either format (zarr 2's Group for v2, zarrv3.Group for v3)."""
    if zarr_version(path) == 3:
        import zarrv3
        return zarrv3.open_group(path)
    import zarr
    return zarr.open_group(str(path), mode="r")


def is_group(n) -> bool:
    import zarr
    import zarrv3
    return isinstance(n, (zarr.hierarchy.Group, zarrv3.Group))


def _axis(a) -> dict:
    if isinstance(a, str):
        return {"name": a.lower(), "type": {"t": "time", "c": "channel", "z": "space", "y": "space", "x": "space"}.get(a.lower(), ""),
                "unit": None}
    return {"name": str(a.get("name", "")).lower(), "type": str(a.get("type") or "").lower(), "unit": a.get("unit")}


def axis_names(axes) -> list[str]:
    return [_axis(a)["name"] for a in axes or []]


# ---------------------------------------------------------------------------------------
# transformations
# ---------------------------------------------------------------------------------------
def to_affine(t: dict, axes_in=None, axes_out=None, what: str = "transformation") -> np.ndarray:
    """One NGFF / SpatialData coordinate transformation as a 3x3 affine on (x, y, 1).

    `axes_in` / `axes_out`: the axis names the transformation's vectors run over (taken from its own
    "input" / "output" when it has them, else these, else x, y)."""
    def names(side, fallback):
        v = t.get(side)
        if isinstance(v, dict) and v.get("axes"):
            return axis_names(v["axes"])
        return list(fallback) if fallback else ["x", "y"]
    ain, aout = names("input", axes_in), names("output", axes_out or axes_in)
    typ = str(t.get("type", "")).lower()
    M = np.eye(3)
    if typ == "identity":
        return M
    if typ == "scale":
        s = [float(v) for v in t.get("scale") or []]
        if len(s) != len(ain) or "x" not in ain or "y" not in ain:
            raise NgffError(f"[ngff] {what}: a scale of {len(s)} values for the axes {ain}: not understood")
        M[0, 0], M[1, 1] = s[ain.index("x")], s[ain.index("y")]
        return M
    if typ == "translation":
        s = [float(v) for v in t.get("translation") or []]
        if len(s) != len(ain) or "x" not in ain or "y" not in ain:
            raise NgffError(f"[ngff] {what}: a translation of {len(s)} values for the axes {ain}: not understood")
        M[0, 2], M[1, 2] = s[ain.index("x")], s[ain.index("y")]
        return M
    if typ == "affine":
        A = np.asarray(t.get("affine"), np.float64)
        if A.ndim != 2 or A.shape[1] != len(ain) + 1 or A.shape[0] not in (len(aout), len(aout) + 1) \
                or not {"x", "y"} <= set(ain) or not {"x", "y"} <= set(aout):
            raise NgffError(f"[ngff] {what}: an affine of shape {list(A.shape)} for the axes {ain} -> {aout}: not understood")
        ix, iy = ain.index("x"), ain.index("y")
        ox, oy = aout.index("x"), aout.index("y")
        others = [i for i in range(len(ain)) if i not in (ix, iy)]
        if any(abs(A[r, i]) > 1e-12 for r in (ox, oy) for i in others):
            raise NgffError(f"[ngff] {what}: the affine mixes x / y with the axes {[ain[i] for i in others]}: not supported")
        for r, rr in ((0, ox), (1, oy)):
            M[r, 0], M[r, 1], M[r, 2] = A[rr, ix], A[rr, iy], A[rr, -1]
        return M
    if typ == "sequence":
        out = np.eye(3)
        for sub in t.get("transformations") or []:
            out = to_affine(sub, ain, aout, what) @ out
        return out
    raise NgffError(f"[ngff] {what}: a transformation of type {str(t.get('type') or '(none)')!r}. SODA reads identity, scale, "
                    f"translation, affine and sequence on the x / y axes; transform the element to one of those first")


def element_transforms(attrs: dict, what: str) -> dict[str, np.ndarray]:
    """{coordinate system name: 3x3 affine from the element's own (intrinsic) pixels / units} of a SpatialData
    element: the transformations that name an "output" coordinate system, from the multiscales of an image or
    the element's own coordinateTransformations."""
    ts = []
    for ms in attrs.get("multiscales") or []:
        ts += [t for t in ms.get("coordinateTransformations") or [] if isinstance(t.get("output"), dict)]
    ts += [t for t in attrs.get("coordinateTransformations") or [] if isinstance(t.get("output"), dict)]
    out = {}
    for t in ts:
        name = str(t["output"].get("name") or "")
        if name:
            out[name] = to_affine(t, what=what)
    return out


def apply(M: np.ndarray, x, y) -> tuple[np.ndarray, np.ndarray]:
    x = np.asarray(x, np.float64)
    y = np.asarray(y, np.float64)
    return M[0, 0] * x + M[0, 1] * y + M[0, 2], M[1, 0] * x + M[1, 1] * y + M[1, 2]


def scale_of(M: np.ndarray) -> float:
    """The mean length scale of an affine (sqrt of |det|)."""
    return float(np.sqrt(abs(np.linalg.det(M[:2, :2]))))


# ---------------------------------------------------------------------------------------
# images
# ---------------------------------------------------------------------------------------
def find_image_group(path, series=None, what: str = "omezarr") -> tuple[Path, dict, str]:
    """(folder of the multiscale image, its attrs, a note) for a path given by the user: the image itself,
    or the top of a bioformats2raw layout (series `series`, default 0)."""
    p = Path(path)
    if zarr_version(p) is None:
        raise NgffError(f"[{what}] not a zarr folder (no .zgroup / .zattrs / zarr.json): an OME-Zarr image is a folder")
    a = read_attrs(p)
    if "plate" in a or "well" in a:
        raise NgffError(f"[{what}] this is an OME-Zarr high-content screening {'plate' if 'plate' in a else 'well'}: "
                        f"SODA reads one image; point at one field's image folder (<well>/<field>)")
    if "multiscales" in a:
        if series is not None and str(series) not in ("0", ""):
            raise NgffError(f"[{what}] --series {series}: this store holds a single image")
        return p, a, ""
    if "bioformats2raw.layout" in a or (p / "OME").is_dir():
        names = (read_attrs(p / "OME").get("series") or
                 sorted((d.name for d in p.iterdir() if d.is_dir() and d.name.isdigit()), key=int))
        names = [str(n) for n in names]
        if not names:
            raise NgffError(f"[{what}] a bioformats2raw layout without a series folder")
        pick = str(series) if series is not None else names[0]
        if pick not in names:
            raise NgffError(f"[{what}] --series {series}: the series here are {', '.join(names)}")
        q = p / pick
        a2 = read_attrs(q)
        if "multiscales" not in a2:
            raise NgffError(f"[{what}] series {pick} has no multiscales metadata")
        note = f"series {pick} of {len(names)}" + (" (--series picks another)" if len(names) > 1 and series is None else "")
        return q, a2, note
    if "image-label" in a:
        raise NgffError(f"[{what}] this is a label image (a segmentation mask), not an image")
    raise NgffError(f"[{what}] no OME-NGFF image here: the .zattrs has no 'multiscales' (and it is not a "
                    f"bioformats2raw layout)")


class NgffImage:
    """A multiscale image: its axes, level 0 and the coarsest level, channel names / colours, pixel size."""

    def __init__(self, folder, attrs: dict, what: str = "omezarr"):
        self.folder, self.what = Path(folder), what
        ms = (attrs.get("multiscales") or [None])[0]
        if not ms or not ms.get("datasets"):
            raise NgffError(f"[{what}] multiscales without datasets")
        self.version = str(ms.get("version") or attrs.get("version") or "")
        ds = ms["datasets"]
        self.level_paths = [str(d["path"]) for d in ds]
        try:
            self.z0 = open_array(self.folder / self.level_paths[0])
        except SystemExit as e:                                 # zarrv3: outside the subset read, with the reason
            raise NgffError(f"[{what}] level 0: {str(e).replace('[zarr v3] ', 'zarr v3: ', 1)}")
        except Exception as e:
            raise NgffError(f"[{what}] level 0 ({self.level_paths[0]}) cannot be opened as a zarr array ({type(e).__name__})")
        try:
            self.zlow = open_array(self.folder / self.level_paths[-1])
        except (Exception, SystemExit):
            self.zlow = self.z0
        self.zarr_format = zarr_version(self.folder / self.level_paths[0])
        axes = ms.get("axes")
        if axes:
            self.axes = [_axis(a) for a in axes]
        elif self.z0.ndim == 5:                                 # NGFF 0.1 / 0.2: always t, c, z, y, x
            self.axes = [_axis(n) for n in "tczyx"]
        else:
            raise NgffError(f"[{what}] no axes in the multiscales metadata and a {self.z0.ndim}-D array: not understood")
        self.names = [a["name"] for a in self.axes]
        if len(self.names) != self.z0.ndim:
            raise NgffError(f"[{what}] {len(self.names)} axes for a {self.z0.ndim}-D array")
        if "x" not in self.names or "y" not in self.names:
            raise NgffError(f"[{what}] the axes {self.names} have no x and y")
        shape = dict(zip(self.names, self.z0.shape))
        self.shape = shape
        self.width, self.height = int(shape["x"]), int(shape["y"])
        for n, a in zip(self.names, self.axes):
            if n in ("x", "y", "c", "z", "t"):
                continue
            if a["type"] == "channel" and "c" not in self.names:
                continue
            if shape[n] > 1:
                raise NgffError(f"[{what}] the axis {n!r} ({a['type'] or 'no type'}) has {shape[n]} entries: not understood")
        if shape.get("t", 1) > 1:
            raise NgffError(f"[{what}] {shape['t']} time points: SODA shows one image; save the time point you want "
                            f"on its own")
        self.cax = "c" if "c" in self.names else next((n for n, a in zip(self.names, self.axes) if a["type"] == "channel"), None)
        self.nC = int(shape[self.cax]) if self.cax else 1
        self.nZ = int(shape.get("z", 1))
        # the pixel size: level 0's scale (dataset and multiscale level) in the x / y unit
        sx = sy = 1.0
        for t in (ds[0].get("coordinateTransformations") or []) + \
                [t for t in (ms.get("coordinateTransformations") or []) if not isinstance(t.get("output"), dict)]:
            if str(t.get("type")) == "scale":
                s = t.get("scale") or []
                if len(s) == len(self.names):
                    sx *= float(s[self.names.index("x")])
                    sy *= float(s[self.names.index("y")])
        ux = self.axes[self.names.index("x")]["unit"]
        uy = self.axes[self.names.index("y")]["unit"]
        f = UNIT_UM.get(str(ux).lower()) if ux else None
        self.pixel_size_um = None
        self.pixel_note = "no unit on the x / y axes" if not ux else ""
        if f and ux == uy:
            px, py = sx * f, sy * f
            if px > 0 and abs(px - py) / px < 0.01:
                self.pixel_size_um = px
            elif px > 0:
                self.pixel_note = f"pixels are not square ({px:g} x {py:g} µm)"
        elif ux and not f:
            self.pixel_note = f"the x axis unit {ux!r} is not a length"
        # channels: omero labels and colours
        om = (attrs.get("omero") or {}).get("channels") or []
        self.labels = [str(c.get("label") if c.get("label") is not None else "") for c in om][:self.nC]
        self.colors = [("#" + str(c.get("color")).lstrip("#").lower()[:6]) if c.get("color") else None for c in om][:self.nC]
        while len(self.labels) < self.nC:
            self.labels.append("")
        while len(self.colors) < self.nC:
            self.colors.append(None)
        self.dtype = np.dtype(self.z0.dtype)
        self.is_rgb = (self.nC == 3 and self.dtype == np.uint8
                       and ({c.lstrip("#") for c in self.colors if c} == _RGB_COLORS
                            or {lab.strip().lower() for lab in self.labels} in _RGB_NAMES
                            or not any(self.labels)))

    def channel_names(self) -> list[str]:
        if self.is_rgb:
            return ["Red", "Green", "Blue"]
        out = []
        for i, lab in enumerate(self.labels):
            out.append(lab.strip() if lab.strip() else f"ch{i}")
        return out

    def source(self, z=None, names=None, log=print) -> "OmeZarrSource":
        return OmeZarrSource(self, z=z, names=names, log=log)


class OmeZarrSource(Source):
    """Level 0 of an OME-Zarr image, one channel and one band of rows at a time (pyramid.Source).

    uint8 / uint16 are passed as they are. Other types are mapped to uint16 per channel: 0 .. the channel's
    99.99th percentile on the coarsest level (floats, int32 / uint32 / int64), negative values to 0; a bool
    mask to 0 / 255."""

    def __init__(self, img: NgffImage, z=None, names=None, log=print):
        self.img = img
        self.width, self.height = img.width, img.height
        self.names = list(names) if names else img.channel_names()
        if len(self.names) != img.nC:
            raise NgffError(f"[{img.what}] {len(self.names)} channel names for {img.nC} channels")
        self.is_rgb = img.is_rgb
        self.pixel_size_um = img.pixel_size_um
        if z is None:
            z = img.nZ // 2
        if not 0 <= int(z) < img.nZ:
            raise NgffError(f"[{img.what}] --z {z}: the image has {img.nZ} z-plane(s) (0 .. {img.nZ - 1})")
        self.z = int(z)
        src = img.dtype
        self._scale = None
        if src in (np.dtype(np.uint8), np.dtype(np.uint16)):
            self.dtype = src
        elif src == np.dtype(bool):
            self.dtype = np.dtype(np.uint8)
            self._scale = [255.0] * img.nC
        else:
            self.dtype = np.dtype(np.uint16)
            self._scale = []
            for c in range(img.nC):
                a = self._read(img.zlow, c, None).astype(np.float64)
                a = a[np.isfinite(a) & (a > 0)]
                top = float(np.percentile(a, 99.99)) if a.size else 1.0
                self._scale.append(65535.0 / top if top > 0 else 1.0)
            log(f"[{img.what}] {src.name} pixels mapped to 16 bits per channel (0 .. the 99.99th percentile)")
        self._lock = threading.Lock()

    def _read(self, arr, c, rows):
        idx = []
        for n, a in zip(self.img.names, self.img.axes):
            if n == "y":
                idx.append(rows if rows is not None else slice(None))
            elif n == "x":
                idx.append(slice(None))
            elif n == self.img.cax:
                idx.append(int(c))
            elif n == "z":
                idx.append(min(self.z, arr.shape[self.img.names.index("z")] - 1))
            else:
                idx.append(0)
        return np.asarray(arr[tuple(idx)])

    def read_band(self, c, y0, y1):
        with self._lock:
            band = self._read(self.img.z0, c, slice(y0, y1))
        if self._scale is None:
            return band
        b = band.astype(np.float64) * self._scale[c]
        b[~np.isfinite(b)] = 0
        return np.clip(b, 0, 255 if self.dtype == np.uint8 else 65535).astype(self.dtype)


def open_image(path, series=None, what: str = "omezarr") -> tuple[NgffImage, str]:
    folder, attrs, note = find_image_group(path, series, what)
    return NgffImage(folder, attrs, what), note
