"""A zipped instrument output, read where it is (xenium.py uses it for a zipped Xenium output).

    zf = zipinput.ZipInput(path, work=<folder beside the cache>, log=print)   # checks the whole file list first
    m = zf.member("experiment.xenium")          # paths are relative to the dataset folder inside the zip
    m.read_bytes() / m.open() / m.random_access() / m.local()

Reading straight from the zip is the default. Each member is reached the cheapest way that works:

    stored (not compressed)     read in place, with real random access (StoredFile: a window on the zip file)
    compressed, read in order   streamed (csv.gz, tsv.gz, matrix.mtx.gz, json)
    compressed, random access   (parquet, h5, a nested transcripts.zarr.zip) read into memory when <= MEMORY_MB
                                (SODA_ZIP_MEMORY_MB), else copied out
    copied out                  the stain images (morphology_focus/*.ome.tif) and a compressed transcripts.zarr.zip:
                                written to the work folder (never next to the zip), only those members, after a
                                free-space check, and said so in the log. `cleanup()` removes them again
                                (only the files this module wrote).

Safety, checked on the zip's file list before anything is read; a zip that fails is refused as a whole:

    zip-slip    a member name that is absolute (/x, C:x, \\\\server\\x), holds a .. part or a NUL is refused;
                every file written is also checked to land inside the work folder
    zip bomb    declared sizes: total uncompressed > SODA_ZIP_MAX_GB (default 2048 GB), a member whose
                declared ratio is above MAX_RATIO (deflate cannot pass 1032:1, so this only catches other
                methods and lies), members that overlap in the file (several names on the same bytes), a
                member that runs past the end of the file. Actual sizes: every member read is counted and must
                end at its declared size (zipfile checks the CRC too); more or less is refused
    disk        before copying out: the declared sizes must fit in the free space of the work folder's drive
                with SODA_UNZIP_RESERVE_GB (default 2) to spare, and stay under SODA_UNZIP_MAX_GB (default 1024)

Messages name members by their path inside the dataset folder, never the zip's own top folder (which is often
the sample's name).
"""
from __future__ import annotations

import collections.abc
import io
import json
import os
import re
import shutil
import struct
import time
import zipfile
from pathlib import Path, PurePosixPath

MEMORY_MB = 256                 # a compressed random-access member up to this size is read into memory (SODA_ZIP_MEMORY_MB)
MAX_RATIO = 1100                # declared uncompressed / compressed, per member above 1 MB
MAX_MEMBERS = 2_000_000
CHUNK = 8 << 20
MARKER = "_copied_from_zip.json"


def _gb_env(name: str, default: float) -> float:
    try:
        return float(os.environ.get(name, "") or default)
    except ValueError:
        return default


class ZipRefused(SystemExit):
    """The zip is not read: unsafe names, sizes that do not add up, not enough room. Nothing was written."""


def _short(e, full_name: str) -> str:
    """An error text with the member's full name (which may start with the sample's folder) cut to its file name."""
    return str(e).replace(full_name, PurePosixPath(full_name).name)


def _shown(name: str) -> str:
    n = name.replace("\\", "/")
    k = n.find("..")
    return n[k:] if k >= 0 else n


def _bad_name(name: str) -> str | None:
    n = name.replace("\\", "/")
    if "\x00" in n:
        return "a NUL character in a name"
    if n.startswith("/") or re.match(r"^[A-Za-z]:", n):
        return "an absolute path"
    if any(p == ".." for p in n.split("/")):
        return "a '..' part (it would write outside the folder)"
    return None


def memory_bytes() -> float:
    return _gb_env("SODA_ZIP_MEMORY_MB", MEMORY_MB) * 1e6


def _fmt_gb(n: float) -> str:
    return f"{n / 1e9:.2f} GB" if n >= 1e8 else f"{n / 1e6:.1f} MB"


class StoredFile(io.RawIOBase):
    """Random access to a stored (uncompressed) member: a read-only window on the zip file itself."""

    def __init__(self, zip_path: Path, info: zipfile.ZipInfo, archive_size: int):
        super().__init__()
        self._f = open(zip_path, "rb")
        self._f.seek(info.header_offset)
        h = self._f.read(zipfile.sizeFileHeader)
        if len(h) != zipfile.sizeFileHeader or h[:4] != zipfile.stringFileHeader:
            self._f.close()
            raise ZipRefused(f"[zip] {info.filename.rsplit('/', 1)[-1]}: no local header where the file list says (damaged zip)")
        fh = struct.unpack(zipfile.structFileHeader, h)
        self._start = info.header_offset + zipfile.sizeFileHeader + fh[zipfile._FH_FILENAME_LENGTH] + fh[zipfile._FH_EXTRA_FIELD_LENGTH]
        self._size = int(info.file_size)
        if self._start + self._size > archive_size:
            self._f.close()
            raise ZipRefused(f"[zip] {info.filename.rsplit('/', 1)[-1]}: declared {self._size:,} bytes, but the zip ends "
                             f"{max(0, archive_size - self._start):,} bytes after its start (truncated or a false size)")
        self._pos = 0
        self.name = info.filename

    def readable(self):
        return True

    def seekable(self):
        return True

    def tell(self):
        return self._pos

    def seek(self, off, whence=io.SEEK_SET):
        base = {io.SEEK_SET: 0, io.SEEK_CUR: self._pos, io.SEEK_END: self._size}[whence]
        self._pos = max(0, base + off)
        return self._pos

    def readinto(self, b):
        n = max(0, min(len(b), self._size - self._pos))
        if not n:
            return 0
        self._f.seek(self._start + self._pos)
        got = self._f.readinto(memoryview(b)[:n])
        self._pos += got
        return got

    def readall(self):
        return self.read(self._size - self._pos)

    def read(self, n=-1):
        if n is None or n < 0:
            n = self._size - self._pos
        buf = bytearray(max(0, min(n, self._size - self._pos)))
        got = self.readinto(buf)
        return bytes(buf[:got])

    def close(self):
        try:
            self._f.close()
        finally:
            super().close()


class _Counted(io.RawIOBase):
    """A member read in order, counted: it must end exactly at its declared size."""

    def __init__(self, raw, info: zipfile.ZipInfo, what: str):
        super().__init__()
        self._r, self._want, self._got, self._what, self._name = raw, int(info.file_size), 0, what, info.filename

    def readable(self):
        return True

    def readinto(self, b):
        try:
            n = self._r.readinto(b)
        except zipfile.BadZipFile as e:
            raise ZipRefused(f"[zip] {self._what}: the data does not match the file list ({_short(e, self._name)}): refused")
        self._got += n
        if self._got > self._want:
            raise ZipRefused(f"[zip] {self._what}: more data than the {self._want:,} bytes declared (refused)")
        if n == 0 and len(b) and self._got != self._want:
            raise ZipRefused(f"[zip] {self._what}: {self._got:,} bytes instead of the {self._want:,} declared (truncated)")
        return n

    def close(self):
        try:
            self._r.close()
        finally:
            super().close()


class Member:
    """One file (or folder) inside the dataset folder of the zip. Has .name and .suffix like a Path."""

    def __init__(self, zi: "ZipInput", rel: str):
        self.zi, self.rel = zi, rel
        self.info = zi._files.get(rel)
        self.name = PurePosixPath(rel).name
        self.suffix = PurePosixPath(rel).suffix
        self.size = int(self.info.file_size) if self.info else 0

    def __repr__(self):
        return f"<{self.rel} in zip>"

    def __str__(self):
        return self.rel

    def is_dir(self) -> bool:
        return self.info is None and self.zi.is_dir(self.rel)

    @property
    def stored(self) -> bool:
        return self.info is not None and self.info.compress_type == zipfile.ZIP_STORED

    def open(self):
        """A binary stream, read in order (counted against the declared size)."""
        if self.stored:
            return StoredFile(self.zi.path, self.info, self.zi.archive_size)
        return io.BufferedReader(_Counted(self.zi.zip.open(self.info), self.info, self.rel), CHUNK)

    def read_bytes(self) -> bytes:
        with self.open() as f:
            return f.read()

    def random_access(self):
        """Something pyarrow / h5py / zipfile can seek in: the member in place when stored, in memory when
        small, else a copy in the work folder (a Path)."""
        if self.stored:
            return StoredFile(self.zi.path, self.info, self.zi.archive_size)
        if self.size <= memory_bytes():
            return io.BytesIO(self.read_bytes())
        return self.local()

    def local(self) -> Path:
        """A copy of this member in the work folder (copied once; checked for room first)."""
        return self.zi.extract([self.rel], why=self.name)[0]


class ZipStoreMap(collections.abc.MutableMapping):
    """A read-only key -> bytes mapping over the members under a prefix of a ZipFile: a zarr v2 store
    (transcripts.zarr inside the zip, or the members of a nested transcripts.zarr.zip)."""

    def __init__(self, z: zipfile.ZipFile, prefix: str = "", closer=None):
        self._z, self._p, self._closer = z, (prefix.rstrip("/") + "/") if prefix else "", closer
        self._keys = {i.filename[len(self._p):]: i for i in z.infolist()
                      if i.filename.startswith(self._p) and not i.is_dir()}

    def __getitem__(self, k):
        i = self._keys.get(k)
        if i is None:
            raise KeyError(k)
        with self._z.open(i) as f:
            data = f.read()
        if len(data) != i.file_size:
            raise ZipRefused(f"[zip] {k}: {len(data):,} bytes instead of the {i.file_size:,} declared")
        return data

    def __contains__(self, k):
        return k in self._keys

    def __iter__(self):
        return iter(self._keys)

    def __len__(self):
        return len(self._keys)

    def keys(self):
        return self._keys.keys()

    def __setitem__(self, k, v):
        raise PermissionError("read-only")

    def __delitem__(self, k):
        raise PermissionError("read-only")

    def close(self):
        if self._closer:
            self._closer()


def zarr_store(member: Member, log=print):
    """A zarr store for transcripts.zarr(.zip) inside the zip, and a function that closes it.
    A folder: read member by member in place. A stored .zarr.zip: opened in place. A compressed one:
    copied out to the work folder first (random access into a compressed member is not practical)."""
    zi = member.zi
    if member.is_dir():
        return ZipStoreMap(zi.zip, zi._full(member.rel)), (lambda: None)
    if member.stored or member.size <= memory_bytes():
        raw = StoredFile(zi.path, member.info, zi.archive_size) if member.stored else io.BytesIO(member.read_bytes())
        try:
            inner = zipfile.ZipFile(raw)
        except zipfile.BadZipFile as e:
            raw.close()
            raise ZipRefused(f"[zip] {member.rel}: not a readable zip ({e})")
        check_list(inner.infolist(), member.size, member.name)
        m = ZipStoreMap(inner, "", closer=lambda: (inner.close(), raw.close()))
        log(f"[zip] {member.rel}: read in place ({'stored in the zip' if member.stored else 'in memory'}, {_fmt_gb(member.size)})")
        return m, m.close
    import zarr
    p = member.local()
    store = zarr.ZipStore(str(p), mode="r")
    return store, store.close


def check_list(infos, archive_size: int, what: str, max_total_gb: float | None = None):
    """The zip-slip / zip-bomb checks on a file list; raises ZipRefused with the first problem."""
    if len(infos) > MAX_MEMBERS:
        raise ZipRefused(f"[zip] {what}: {len(infos):,} files inside (more than {MAX_MEMBERS:,}): refused")
    cap = (max_total_gb if max_total_gb is not None else _gb_env("SODA_ZIP_MAX_GB", 2048)) * 1e9
    total = 0
    for i in infos:
        bad = _bad_name(i.filename)
        if bad:
            raise ZipRefused(f"[zip] {what}: refused, a file inside has {bad}: {_shown(i.filename)!r}. Nothing was read or written")
        if i.flag_bits & 0x1:
            raise ZipRefused(f"[zip] {what}: encrypted files inside; unzip it yourself (with the password) and use the folder")
        total += int(i.file_size)
        if i.file_size > 1e6 and i.file_size > MAX_RATIO * max(1, i.compress_size):
            raise ZipRefused(f"[zip] {what}: {PurePosixPath(i.filename).name} declares {i.file_size:,} bytes from "
                             f"{i.compress_size:,} (ratio above {MAX_RATIO}): refused as a possible zip bomb")
        if i.compress_type == zipfile.ZIP_STORED and i.file_size != i.compress_size:
            raise ZipRefused(f"[zip] {what}: {PurePosixPath(i.filename).name} is stored but declares two sizes "
                             f"({i.file_size:,} / {i.compress_size:,}): refused")
    if total > cap:
        raise ZipRefused(f"[zip] {what}: the files inside declare {_fmt_gb(total)} in all, above the limit of "
                         f"{_fmt_gb(cap)} (SODA_ZIP_MAX_GB): refused as a possible zip bomb")
    spans = sorted((int(i.header_offset), int(i.header_offset) + zipfile.sizeFileHeader + len(i.orig_filename.encode("utf-8", "replace"))
                    + int(i.compress_size), i.filename) for i in infos)
    for (a0, a1, an), (b0, _b1, bn) in zip(spans, spans[1:]):
        if b0 < a1:
            raise ZipRefused(f"[zip] {what}: two files inside share the same bytes ({PurePosixPath(an).name}, "
                             f"{PurePosixPath(bn).name}): refused as a possible zip bomb")
    if spans and spans[-1][1] > archive_size:
        raise ZipRefused(f"[zip] {what}: {PurePosixPath(spans[-1][2]).name} runs past the end of the zip "
                         f"(truncated, or a false size): refused")
    return total


class ZipInput:
    """A zip holding one dataset folder (at the top of the zip or a few folders down). `root_test(names)`
    says whether a folder's file / folder names make it the dataset folder; without it the zip's top is used."""

    def __init__(self, path, work, log=print, root_test=None, tag="zip"):
        self.path = Path(path)
        self.work = Path(work)
        self.log, self.tag = log, tag
        try:
            self.archive_size = self.path.stat().st_size
            self.zip = zipfile.ZipFile(str(self.path))
        except (zipfile.BadZipFile, OSError) as e:
            raise ZipRefused(f"[{tag}] {self.path.name}: not a readable zip ({e})")
        infos = self.zip.infolist()
        self.total = check_list(infos, self.archive_size, self.path.name)
        self._all = {i.filename.replace("\\", "/").strip("/"): i for i in infos
                     if not i.is_dir() and not i.filename.startswith("__MACOSX/")}
        dirs: dict[str, set] = {"": set()}
        for n in self._all:
            parts = n.split("/")
            for k in range(len(parts)):
                d = "/".join(parts[:k])
                dirs.setdefault(d, set()).add(parts[k])
        self._dirs = dirs
        roots = [d for d, names in dirs.items() if root_test(names)] if root_test else [""]
        if not roots:
            raise ZipRefused(f"[{tag}] {self.path.name}: no dataset folder found inside")
        top = min(len(r.split("/")) if r else 0 for r in roots)
        shallow = [r for r in roots if (len(r.split("/")) if r else 0) == top]
        if len(shallow) > 1:
            raise ZipRefused(f"[{tag}] {self.path.name}: {len(shallow)} dataset folders inside ({', '.join(sorted(shallow)[:4])}): "
                             f"unzip it and point at the one you want")
        self.root = shallow[0]
        self.depth = top
        self._files = {n[len(self.root) + 1:] if self.root else n: i for n, i in self._all.items()
                       if not self.root or n.startswith(self.root + "/")}
        self._extracted: list[Path] = []
        self._made_work = self.work.parent.name == "_unzipped"     # the default place (<out>/_unzipped/<id>) is ours
        self.copied_bytes = 0

    # ---- the dataset folder's view
    def _full(self, rel: str) -> str:
        return f"{self.root}/{rel}" if self.root else rel

    def where(self) -> str:
        return "at the top of the zip" if not self.depth else (
            "in " + self.root.rsplit("/", 1)[-1] + "/" if self.root.rsplit("/", 1)[-1].lower() in ("data", "outs", "output")
            else f"{self.depth} folder{'s' if self.depth > 1 else ''} down")

    def names(self, rel: str = "") -> set:
        return set(self._dirs.get(self._full(rel) if rel else self.root, set()))

    def exists(self, rel: str) -> bool:
        return rel in self._files or self.is_dir(rel)

    def is_dir(self, rel: str) -> bool:
        return self._full(rel) in self._dirs and rel not in self._files

    def member(self, rel: str) -> Member | None:
        return Member(self, rel) if self.exists(rel) else None

    def glob(self, folder: str, pattern: str) -> list[Member]:
        rx = re.compile(pattern, re.I)
        return [Member(self, f"{folder}/{n}") for n in sorted(self.names(folder)) if f"{folder}/{n}" in self._files and rx.search(n)]

    # ---- copying out
    def _marker(self) -> dict:
        try:
            return json.loads((self.work / MARKER).read_text(encoding="utf-8"))
        except (OSError, ValueError):
            return {}

    def extract(self, rels: list[str], why: str = "") -> list[Path]:
        """Copy these members into the work folder (once: a copy from an earlier run of the same zip is
        reused). The free space is checked first; sizes and CRCs while copying."""
        st = self.path.stat()
        stamp = {"zipSize": st.st_size, "zipMtime": int(st.st_mtime)}
        mk = self._marker()
        if {k: mk.get(k) for k in stamp} != stamp:
            mk = dict(stamp, files={})
        out, todo = [], []
        root = self.work.resolve()
        for rel in rels:
            info = self._files[rel]
            dst = (self.work / rel).resolve()
            if root != dst and root not in dst.parents:
                raise ZipRefused(f"[{self.tag}] {rel}: would be written outside the work folder: refused")
            out.append(dst)
            have = mk.get("files", {}).get(rel)
            if dst.exists() and have == {"size": int(info.file_size), "crc": int(info.CRC)} and dst.stat().st_size == info.file_size:
                if dst not in self._extracted:
                    self._extracted.append(dst)
                continue
            todo.append((rel, info, dst))
        if not todo:
            if rels:
                self.log(f"[{self.tag}] {', '.join(PurePosixPath(r).name for r in rels[:4])}{' ...' if len(rels) > 4 else ''}: "
                         f"the copy from an earlier run is used ({self.work})")
            return out
        need = sum(int(i.file_size) for _, i, _ in todo)
        cap = _gb_env("SODA_UNZIP_MAX_GB", 1024) * 1e9
        if need > cap:
            raise ZipRefused(f"[{self.tag}] copying {_fmt_gb(need)} out of {self.path.name} is above the limit of {_fmt_gb(cap)} "
                             f"(SODA_UNZIP_MAX_GB): refused. Nothing was written")
        probe = self.work
        while not probe.exists() and probe.parent != probe:
            probe = probe.parent
        free = shutil.disk_usage(str(probe)).free
        reserve = _gb_env("SODA_UNZIP_RESERVE_GB", 2) * 1e9
        if need + reserve > free:
            raise ZipRefused(f"[{self.tag}] not enough free space to copy {_fmt_gb(need)} out of {self.path.name} into "
                             f"{self.work}: {_fmt_gb(free)} free, {_fmt_gb(need + reserve)} needed with "
                             f"{_fmt_gb(reserve)} to spare (SODA_UNZIP_RESERVE_GB). Free some space or give --unzip-to <folder "
                             f"on another drive>. Nothing was written")
        if not self.work.exists():
            self._made_work = True
            self.work.mkdir(parents=True)
        t0 = time.time()
        self.log(f"[{self.tag}] copying {len(todo)} file(s), {_fmt_gb(need)}, out of the zip into the work folder {self.work} "
                 f"({why or 'needed whole'}; {_fmt_gb(free)} free there). Your zip is not changed")
        for rel, info, dst in todo:
            dst.parent.mkdir(parents=True, exist_ok=True)
            part = dst.with_name(dst.name + ".part")
            got = 0
            with self.zip.open(info) as src, open(part, "wb") as f:
                while True:
                    try:
                        b = src.read(CHUNK)
                    except zipfile.BadZipFile as e:
                        f.close()
                        part.unlink(missing_ok=True)
                        raise ZipRefused(f"[{self.tag}] {rel}: the data does not match the file list ({_short(e, info.filename)}): refused")
                    if not b:
                        break
                    got += len(b)
                    if got > info.file_size:
                        f.close()
                        part.unlink(missing_ok=True)
                        raise ZipRefused(f"[{self.tag}] {rel}: more data than the {info.file_size:,} bytes declared: refused")
                    f.write(b)
            if got != info.file_size:
                part.unlink(missing_ok=True)
                raise ZipRefused(f"[{self.tag}] {rel}: {got:,} bytes instead of the {info.file_size:,} declared (truncated): refused")
            os.replace(part, dst)
            self._extracted.append(dst)
            mk.setdefault("files", {})[rel] = {"size": int(info.file_size), "crc": int(info.CRC)}
            (self.work / MARKER).write_text(json.dumps(mk, indent=1), encoding="utf-8")
            self.copied_bytes += got
        dt = time.time() - t0
        self.log(f"[{self.tag}]   copied {_fmt_gb(need)} in {dt:.0f}s")
        return out

    def cleanup(self, keep: bool = False):
        """Remove the copies this run made or reused (never anything else), then the work folder if empty."""
        if keep or not self._extracted:
            if self._extracted:
                self.log(f"[{self.tag}] the files copied out of the zip are kept in {self.work} (--keep-unzipped)")
            return
        n = 0
        for p in self._extracted:
            try:
                p.unlink()
                n += 1
            except OSError:
                pass
        (self.work / MARKER).unlink(missing_ok=True)
        w = self.work.resolve()
        for d in sorted({q for p in self._extracted for q in p.parents if w in q.parents}, key=lambda q: -len(q.parts)):
            try:
                d.rmdir()                                  # sub-folders this module made, only when empty
            except OSError:
                pass
        for d in ([w] if self._made_work else []) + ([w.parent] if self._made_work and w.parent.name == "_unzipped" else []):
            try:
                d.rmdir()                                  # the work folder (and <out>/_unzipped) if this run made it and it is empty
            except OSError:
                pass
        self.log(f"[{self.tag}] removed the {n} file(s) copied out of the zip (the work folder's copies only; "
                 f"--keep-unzipped keeps them)")
        self._extracted = []

    def close(self):
        self.zip.close()
