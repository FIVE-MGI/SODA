r"""SODA desktop app: one program that prepares a folder of datasets and opens the viewer.

    SODA                                    double-click: the SODA home page opens in your browser
    SODA prepare --src <data folder> --work <output folder> [any prepare.py batch option]
                                            every dataset under <data folder>, exported to
                                            <output folder>/export (cache/ and logs/ beside it),
                                            then datasets.json and the sample ID check (always a gate)
    SODA view <export folder> [--port N] [--no-browser]
                                            serve the page and that folder on 127.0.0.1 and open the viewer
    SODA check <folder or file> ...         the sample ID check (deid.py check); exit 1 if it finds any
    SODA selftest [--quick]                 import and exercise the packaged libraries; JSON on stdout
    SODA version
    SODA <script>.py [args]                 run one of SODA's scripts as python would, e.g.
                                            SODA prepare.py detect <folder>   or   SODA deid.py check <folder>
                                            (a bare name is looked up among SODA's own scripts)

From source the same commands are  python app/soda_app.py ...  with the repository's venv.

How it is built (docs/DESKTOP_APP.md): PyInstaller one-folder, the same recipe as FIVE Atlas.
SODA's own .py files ship as plain source next to the page (web/SpatialViewer is copied whole,
minus private files), so prepare.py, batch.py and deid.py run exactly as they do from a checkout.
batch.py starts every dataset in its own process with sys.executable; frozen, that is this
program, which is why it also behaves like python for `SODA <script>.py`.

Nothing is uploaded anywhere: the server listens on 127.0.0.1 only, the actions on the home page
need a per-launch token, and the request log is off (paths can carry dataset names).
Environment: SODA_PORT (default 8770), SODA_NO_BROWSER=1, SODA_WORKDIR (logs and caches; default
%LOCALAPPDATA%\SODA or ~/Library/Application Support/SODA), SODA_NO_UPDATE_CHECK=1.
"""
from __future__ import annotations

import json
import os
import re
import secrets
import socket
import subprocess
import sys
import threading
import time
import urllib.parse
import urllib.request
import webbrowser
from http.server import ThreadingHTTPServer
from pathlib import Path

FROZEN = bool(getattr(sys, "frozen", False))
if FROZEN:
    BASE = Path(getattr(sys, "_MEIPASS", Path(sys.executable).parent))
    SODA_DIR = BASE / "soda"
    HOME_HTML = BASE / "app" / "home.html"
else:
    _HERE = Path(__file__).resolve().parent
    SODA_DIR = _HERE.parent / "web" / "SpatialViewer"
    HOME_HTML = _HERE / "home.html"

APP_NAME = "SODA"
DEFAULT_PORT = 8770
IDLE_EXIT_S = 30 * 60          # macOS only: no Dock icon, so an idle server stops itself
TEXT_TYPES = {
    ".html": "text/html; charset=utf-8", ".js": "text/javascript; charset=utf-8",
    ".mjs": "text/javascript; charset=utf-8", ".css": "text/css; charset=utf-8",
    ".json": "application/json; charset=utf-8", ".geojson": "application/json; charset=utf-8",
    ".svg": "image/svg+xml", ".png": "image/png", ".jpg": "image/jpeg", ".jpeg": "image/jpeg",
    ".webp": "image/webp", ".ico": "image/x-icon", ".bin": "application/octet-stream",
    ".gz": "application/gzip", ".wasm": "application/wasm", ".csv": "text/csv; charset=utf-8",
    ".txt": "text/plain; charset=utf-8", ".woff2": "font/woff2",
}
# The viewer's config.js in the app (no site config is shipped): the page's defaults.
APP_CONFIG_JS = b"window.SpatialViewerConfig = window.SpatialViewerConfig || {};\n"
# Files next to the page that are never served: the tools, not the page.
_BLOCKED_SUFFIXES = (".py", ".pyc", ".ps1", ".sh", ".txt", ".md", ".password", ".csv")
_BLOCKED_PARTS = ("__pycache__", ".venv", "soda-env", ".password")


# ---------------------------------------------------------------------------------------------
# small helpers
# ---------------------------------------------------------------------------------------------
def app_version() -> str:
    if FROZEN:
        try:
            from _buildinfo import VERSION    # written by SODA.spec at build time
            return VERSION
        except Exception:
            return "0.0.0-dev"
    # from source: the checkout's own version, never a leftover app/_buildinfo.py from a local build
    try:
        v = subprocess.run(["git", "-C", str(SODA_DIR), "describe", "--tags", "--dirty", "--always"],
                           capture_output=True, text=True, timeout=5).stdout.strip()
        if v:
            return v.lstrip("v") + "-source"
    except Exception:
        pass
    return "0.0.0-dev"


def update_repo() -> str:
    if FROZEN:
        try:
            from _buildinfo import UPDATE_REPO
            return os.environ.get("SODA_UPDATE_REPO", UPDATE_REPO)
        except Exception:
            pass
    return os.environ.get("SODA_UPDATE_REPO", "")


def workdir() -> Path:
    if os.environ.get("SODA_WORKDIR"):
        d = Path(os.environ["SODA_WORKDIR"])
    elif sys.platform == "win32":
        d = Path(os.environ.get("LOCALAPPDATA") or Path.home() / "AppData" / "Local") / APP_NAME
    elif sys.platform == "darwin":
        d = Path.home() / "Library" / "Application Support" / APP_NAME
    else:
        d = Path(os.environ.get("XDG_DATA_HOME") or Path.home() / ".local" / "share") / APP_NAME
    d.mkdir(parents=True, exist_ok=True)
    return d


def soda_import(name: str):
    """One of SODA's own modules (deid, batch ...), from the page folder."""
    if str(SODA_DIR) not in sys.path:
        sys.path.insert(0, str(SODA_DIR))
    import importlib
    return importlib.import_module(name)


def scrub(text: str) -> str:
    """Anything deid.py would call a sample ID -> [sample ID] (batch.scrub). Nothing this app
    shows or writes skips it."""
    try:
        return soda_import("batch").scrub(text)
    except Exception:
        return text


def python_cmd() -> list[str]:
    """How to start a SODA script: this program when frozen (it acts as python), else python."""
    return [sys.executable]


def self_cmd() -> list[str]:
    """How to start this app again (for the folder picker)."""
    return [sys.executable] if FROZEN else [sys.executable, str(Path(__file__).resolve())]


def setup_streams():
    """UTF-8 out, line by line: batch.py reads a child's output as UTF-8, live. A macOS .app has
    no console at all (sys.stdout is None); the launcher's own text then goes to the log file."""
    for name in ("stdout", "stderr"):
        s = getattr(sys, name)
        if s is None:
            continue
        try:
            s.reconfigure(encoding="utf-8", errors="replace", line_buffering=True)
        except Exception:
            pass


class _Tee:
    def __init__(self, stream, handle):
        self._s, self._h = stream, handle

    def write(self, text):
        for t in (self._s, self._h):
            try:
                if t is not None:
                    t.write(text)
                    t.flush()
            except Exception:
                pass
        return len(text)

    def flush(self):
        for t in (self._s, self._h):
            try:
                if t is not None:
                    t.flush()
            except Exception:
                pass

    def isatty(self):
        return False


def start_log(wd: Path):
    """The launcher's own messages also go to <workdir>/SODA.log (the only place they exist in a
    macOS .app). Only this app's banner and errors are written there, never a request log."""
    try:
        h = open(wd / "SODA.log", "a", encoding="utf-8", buffering=1)
        h.write(f"\n==== {time.strftime('%Y-%m-%d %H:%M:%S')} {APP_NAME} {app_version()} ====\n")
        sys.stdout = _Tee(sys.stdout, h)
        sys.stderr = _Tee(sys.stderr, h)
    except Exception:
        pass


def free_port(preferred: int) -> int:
    for port in (preferred, 0):
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
            try:
                s.bind(("127.0.0.1", port))
                return s.getsockname()[1]
            except OSError:
                continue
    raise RuntimeError("no free port")


def ping(port: int) -> dict | None:
    try:
        with urllib.request.urlopen(f"http://127.0.0.1:{port}/app/api/ping", timeout=1.5) as r:
            d = json.loads(r.read().decode())
            return d if d.get("app") == APP_NAME else None
    except Exception:
        return None


def open_browser(url: str):
    if os.environ.get("SODA_NO_BROWSER") or os.environ.get("CI"):
        return
    try:
        webbrowser.open(url)
    except Exception:
        pass


def _version_key(v: str):
    out = []
    for part in str(v).lstrip("v").split("-")[0].split("."):
        m = re.match(r"\d+", part)
        out.append(int(m.group(0)) if m else 0)
    return tuple((out + [0, 0, 0])[:3])


def check_update(current: str) -> dict | None:
    """Is a newer SODA on the Releases page of the repository it was built from? One 4-second
    request at start-up; offline or a private repository is simply 'no news'. Only the version is
    sent (as the User-Agent)."""
    repo = update_repo()
    if not repo or os.environ.get("SODA_NO_UPDATE_CHECK") or current.endswith("-dev"):
        return None
    try:
        req = urllib.request.Request(f"https://api.github.com/repos/{repo}/releases/latest",
                                     headers={"Accept": "application/vnd.github+json",
                                              "User-Agent": f"SODA/{current}"})
        with urllib.request.urlopen(req, timeout=4) as r:
            d = json.loads(r.read().decode())
        latest = str(d.get("tag_name") or "").lstrip("v")
        if latest and _version_key(latest) > _version_key(current):
            return {"latest": latest, "url": d.get("html_url") or f"https://github.com/{repo}/releases/latest"}
    except Exception:
        pass
    return None


# ---------------------------------------------------------------------------------------------
# prepare = prepare.py batch, with the sample ID check always a gate
# ---------------------------------------------------------------------------------------------
def batch_argv(src: Path, work: Path, extra: list[str]) -> list[str]:
    """prepare.py batch arguments for <work>/export, <work>/cache and <work>/logs.
    --check-ids is always added: the app never leaves the sample ID check as a mere warning."""
    src, work = Path(src).expanduser().resolve(), Path(work).expanduser().resolve()
    if not src.exists():
        raise ValueError(f"the data folder does not exist: {src}")
    if work == src or work.is_relative_to(src):
        raise ValueError("the output folder must not be inside the data folder (it would be read as a dataset)")
    if src.is_relative_to(work):
        raise ValueError("the data folder must not be inside the output folder")
    argv = ["batch", "--src", str(src), "--out", str(work / "cache"), "--dest", str(work / "export"),
            "--logs", str(work / "logs"), *extra]
    if "--check-ids" not in argv:
        argv.append("--check-ids")
    return argv


def read_report(logs: Path, since: float) -> dict | None:
    """The newest batch_*.json that batch.py wrote in <logs> after `since` (already scrubbed)."""
    try:
        cands = [p for p in logs.glob("batch_*.json") if p.stat().st_mtime >= since - 1]
        if not cands:
            return None
        p = max(cands, key=lambda q: q.stat().st_mtime)
        d = json.loads(p.read_text(encoding="utf-8"))
        rows = [{"id": r.get("id"), "kind": r.get("kind"), "status": r.get("status"), "reason": r.get("reason")}
                for r in d.get("datasets", [])]
        return {"report": str(p.with_suffix(".md")), "counts": d.get("counts"), "exitCode": d.get("exitCode"),
                "deid": d.get("deid"), "catalog": (d.get("catalog") or {}).get("text"),
                "export": d.get("export"), "datasets": rows}
    except Exception:
        return None


class Job:
    """One prepare.py batch run in its own process; its output (scrubbed) is kept for the page."""

    MAX_LINES = 4000

    def __init__(self, argv: list[str], work: Path, log_path: Path):
        self.argv, self.work, self.started = argv, Path(work), time.time()
        self.lines: list[str] = []
        self.offset = 0                        # index of self.lines[0] in the whole output
        self.exit: int | None = None
        self.report: dict | None = None
        self.stopped = False
        self._lock = threading.Lock()
        self._log = open(log_path, "a", encoding="utf-8", buffering=1)
        self._log.write(f"\n==== {time.strftime('%Y-%m-%d %H:%M:%S')} prepare (batch) ====\n")
        kw = {}
        if sys.platform == "win32":
            kw["creationflags"] = subprocess.CREATE_NEW_PROCESS_GROUP
        else:
            kw["start_new_session"] = True
        self.proc = subprocess.Popen(
            [*python_cmd(), str(SODA_DIR / "prepare.py"), *argv], stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL, text=True, encoding="utf-8", errors="replace",
            cwd=str(self.work), env={**os.environ, "PYTHONIOENCODING": "utf-8", "PYTHONUNBUFFERED": "1"}, **kw)
        threading.Thread(target=self._pump, daemon=True).start()

    def _pump(self):
        for line in self.proc.stdout:
            line = scrub(line.rstrip("\n"))
            with self._lock:
                self.lines.append(line)
                if len(self.lines) > self.MAX_LINES:
                    drop = len(self.lines) - self.MAX_LINES
                    del self.lines[:drop]
                    self.offset += drop
            try:
                self._log.write(line + "\n")
            except Exception:
                pass
        rc = self.proc.wait()
        self.report = read_report(self.work / "logs", self.started)
        self.exit = rc
        try:
            self._log.write(f"==== exit {rc}\n")
            self._log.close()
        except Exception:
            pass

    @property
    def running(self) -> bool:
        return self.exit is None

    def stop(self):
        self.stopped = True
        if self.proc.poll() is not None:
            return
        if sys.platform == "win32":
            subprocess.run(["taskkill", "/T", "/F", "/PID", str(self.proc.pid)], capture_output=True)
        else:
            import signal
            try:
                os.killpg(self.proc.pid, signal.SIGINT)      # batch.py stops its children on Ctrl+C
            except OSError:
                pass

    def view(self, since: int = 0) -> dict:
        with self._lock:
            start = max(since - self.offset, 0)
            new = self.lines[start:]
            nxt = self.offset + len(self.lines)
        return {"running": self.running, "exit": self.exit, "stopped": self.stopped, "next": nxt,
                "lines": new, "report": self.report, "work": str(self.work),
                "seconds": round(time.time() - self.started)}


# ---------------------------------------------------------------------------------------------
# the local server: home page, the viewer page, export folders, the home page's actions
# ---------------------------------------------------------------------------------------------
class App:
    def __init__(self, port: int, wd: Path):
        self.port, self.wd = port, wd
        self.token = secrets.token_urlsafe(24)
        self.version = app_version()
        self.exports: dict[str, Path] = {}
        self.job: Job | None = None
        self.update: dict | None = None
        self.last_seen = time.time()
        self.httpd: ThreadingHTTPServer | None = None
        self._lock = threading.Lock()

    def base(self) -> str:
        return f"http://127.0.0.1:{self.port}"

    def add_export(self, folder: Path) -> str:
        folder = Path(folder).expanduser().resolve()
        with self._lock:
            for k, v in self.exports.items():
                if v == folder:
                    return k
            k = str(len(self.exports) + 1)
            self.exports[k] = folder
            return k

    def viewer_url(self, folder: str | Path) -> str:
        """The page's address for an export folder: ?data= when it has datasets.json (every dataset
        in the list), the folder of one dataset as ?data= too (the page reads its manifest.json),
        else ?static=&ds= for the only dataset in it."""
        if not str(folder).strip():
            raise ValueError("pick the export folder")
        folder = Path(str(folder).strip()).expanduser().resolve()
        if not folder.is_dir():
            raise ValueError(f"not a folder: {folder}")
        sets = sorted(p.name for p in folder.iterdir() if (p / "manifest.json").is_file())
        if not ((folder / "datasets.json").is_file() or (folder / "manifest.json").is_file() or sets):
            raise ValueError("no SODA export in this folder (no datasets.json and no dataset folder with a "
                             "manifest.json). Pick the 'export' folder that Prepare made.")
        k = self.add_export(folder)
        data = f"{self.base()}/data/{k}/"
        if (folder / "datasets.json").is_file() or (folder / "manifest.json").is_file() or len(sets) != 1:
            q = "data=" + urllib.parse.quote(data, safe=":/")
        else:
            q = "static=" + urllib.parse.quote(data.rstrip("/"), safe=":/") + "&ds=" + urllib.parse.quote(sets[0])
        return f"{self.base()}/viewer/index.html?{q}"

    def status(self) -> dict:
        key = Path(os.environ.get("SPATIALVIZ_KEY") or Path.home() / "private" / "soda" / "sample_id_key.csv")
        return {"app": APP_NAME, "version": self.version, "frozen": FROZEN, "port": self.port,
                "keyFile": str(key), "keyFound": key.exists(), "update": self.update,
                "jobRunning": bool(self.job and self.job.running), "workdir": str(self.wd)}

    def serve(self):
        app = self

        # Files are served by static_serve.py's handler (one implementation of Range / 206 /
        # Content-Range for the gene packs, shared with `python static_serve.py`). This subclass
        # only routes /viewer/ and /data/<n>/ to their folders and tightens it for an app: no
        # CORS (the page is same-origin here; an open CORS header would let any website read
        # the export while SODA runs), no directory listings, no request log.
        Base = soda_import("static_serve").CORSHandler

        class Handler(Base):
            server_version = "SODA"
            sys_version = ""

            def log_message(self, fmt, *args):       # no request log: paths can name datasets
                pass

            def send_header(self, keyword, value):
                if keyword.lower().startswith("access-control-"):
                    return
                super().send_header(keyword, value)

            def end_headers(self):
                self.send_header("X-Content-Type-Options", "nosniff")
                super().end_headers()

            def guess_type(self, path):              # not the Windows registry's idea of .js
                return TEXT_TYPES.get(Path(str(path)).suffix.lower()) or super().guess_type(path)

            def list_directory(self, path):
                self.send_error(404, "not found")
                return None

            def do_OPTIONS(self):                    # same-origin only: no preflight to answer
                self._send(405, b"", "text/plain")

            # -- plumbing ------------------------------------------------------------------
            def _host_ok(self) -> bool:
                host = (self.headers.get("Host") or "").lower()
                return host in (f"127.0.0.1:{app.port}", f"localhost:{app.port}")

            def _send(self, code: int, body: bytes, ctype: str, extra: dict | None = None, head=False):
                self.send_response(code)
                self.send_header("Content-Type", ctype)
                self.send_header("Content-Length", str(len(body)))
                for k, v in (extra or {}).items():
                    self.send_header(k, v)
                self.end_headers()
                if not head:
                    self.wfile.write(body)

            def _json(self, obj, code=200):
                self._send(code, json.dumps(obj).encode(), "application/json; charset=utf-8")

            def _file(self, root: Path, rel: str, head: bool, blocked=False):
                """static_serve's GET / HEAD (Range included) on <root>/<rel>. Its translate_path
                drops '..' and drive parts, so nothing outside <root> is reachable."""
                low = "/" + urllib.parse.unquote(rel).lower()
                if blocked and (low.endswith(_BLOCKED_SUFFIXES) or any(p in low for p in _BLOCKED_PARTS)):
                    return self._send(404, b"not found", "text/plain", head=head)
                if rel and not rel.endswith("/") and (Path(root) / urllib.parse.unquote(rel)).is_dir():
                    return self._send(404, b"not found", "text/plain", head=head)   # no redirect: wrong prefix
                self.directory = str(root)
                self.path = "/" + rel
                return super().do_HEAD() if head else super().do_GET()

            def _body(self) -> dict:
                n = int(self.headers.get("Content-Length") or 0)
                if n > 1 << 20:
                    raise ValueError("request too large")
                raw = self.rfile.read(n) if n else b"{}"
                return json.loads(raw.decode("utf-8") or "{}")

            # -- GET -----------------------------------------------------------------------
            def do_HEAD(self):
                self.do_GET(head=True)

            def do_GET(self, head=False):
                if not self._host_ok():
                    return self._send(403, b"forbidden", "text/plain", head=head)
                app.last_seen = time.time()
                url = urllib.parse.urlsplit(self.path)
                path = url.path
                try:
                    if path == "/":
                        return self._send(302, b"", "text/plain", {"Location": "/app/"}, head=head)
                    if path in ("/app", "/app/", "/app/index.html"):
                        html = HOME_HTML.read_text(encoding="utf-8")
                        html = html.replace("{{SODA_TOKEN}}", app.token).replace("{{SODA_VERSION}}", app.version)
                        return self._send(200, html.encode("utf-8"), TEXT_TYPES[".html"], head=head)
                    if path == "/app/api/ping":
                        return self._json({"app": APP_NAME, "version": app.version})
                    if path == "/app/api/status":
                        return self._json(app.status())
                    if path == "/app/api/job":
                        if not app.job:
                            return self._json({"none": True})
                        since = int(urllib.parse.parse_qs(url.query).get("since", ["0"])[0] or 0)
                        return self._json(app.job.view(since))
                    if path == "/viewer/config.js" and not (SODA_DIR / "config.js").is_file():
                        # config.js is each site's own file and is not shipped; local viewing needs only the defaults
                        return self._send(200, APP_CONFIG_JS, TEXT_TYPES[".js"], head=head)
                    if path.startswith("/viewer/"):
                        return self._file(SODA_DIR, path[len("/viewer/"):], head, blocked=True)
                    m = re.match(r"^/data/(\d+)/(.*)$", path)
                    if m and m.group(1) in app.exports:
                        return self._file(app.exports[m.group(1)], m.group(2), head)
                    if path == "/favicon.ico":
                        return self._send(204, b"", "image/x-icon", head=head)
                    return self._send(404, b"not found", "text/plain", head=head)
                except (BrokenPipeError, ConnectionResetError, ConnectionAbortedError):
                    pass

            # -- POST: the home page's actions -----------------------------------------------
            def do_POST(self):
                if not self._host_ok() or not secrets.compare_digest(
                        self.headers.get("X-SODA-Token", ""), app.token):
                    return self._json({"error": "forbidden"}, 403)
                app.last_seen = time.time()
                path = urllib.parse.urlsplit(self.path).path
                try:
                    b = self._body()
                    if path == "/app/api/prepare":
                        return self._json(app.start_prepare(b))
                    if path == "/app/api/stop":
                        if app.job:
                            app.job.stop()
                        return self._json({"ok": True})
                    if path == "/app/api/view":
                        return self._json({"url": app.viewer_url(str(b.get("folder", "")))})
                    if path == "/app/api/check":
                        return self._json(check_folder(str(b.get("folder", ""))))
                    if path == "/app/api/pick":
                        return self._json({"path": pick(str(b.get("kind", "folder")))})
                    if path == "/app/api/quit":
                        self._json({"ok": True})
                        threading.Thread(target=app.shutdown, daemon=True).start()
                        return
                    return self._json({"error": "unknown action"}, 404)
                except ValueError as e:
                    return self._json({"error": scrub(str(e))}, 400)
                except Exception as e:
                    return self._json({"error": scrub(f"{type(e).__name__}: {e}")}, 500)

        self.httpd = ThreadingHTTPServer(("127.0.0.1", self.port), Handler)
        self.httpd.daemon_threads = True
        return self.httpd

    def start_prepare(self, b: dict) -> dict:
        with self._lock:
            if self.job and self.job.running:
                raise ValueError("a prepare run is still going; stop it first")
            src, work = str(b.get("src", "")).strip(), str(b.get("work", "")).strip()
            if not src or not work:
                raise ValueError("pick both the data folder and the output folder")
            extra: list[str] = []
            if str(b.get("options", "")).strip():
                extra += ["--options", str(Path(str(b["options"]).strip()).expanduser().resolve())]
            if str(b.get("title", "")).strip():
                extra += ["--title", str(b["title"]).strip()]
            jobs = int(b.get("jobs") or 1)
            if jobs > 1:
                extra += ["--jobs", str(min(jobs, 8))]
            if b.get("rebuild"):
                extra.append("--rebuild")
            if b.get("ids") is False:                  # the barcode / cell id lists are written unless unticked
                extra.append("--no-ids")
            argv = batch_argv(Path(src), Path(work), extra)
            Path(work).expanduser().resolve().mkdir(parents=True, exist_ok=True)
            self.job = Job(argv, Path(work).expanduser().resolve(), self.wd / "prepare.log")
            return {"ok": True}

    def shutdown(self):
        if self.job and self.job.running:
            self.job.stop()
        time.sleep(0.3)
        if self.httpd:
            self.httpd.shutdown()


def check_folder(folder: str | Path) -> dict:
    """The sample ID check on one folder (deid.check_paths). What it finds is shown with the IDs
    themselves replaced by [sample ID]: the place is what you need, not the ID."""
    if not str(folder).strip():
        raise ValueError("pick a folder to check")
    folder = Path(str(folder).strip()).expanduser()
    if not folder.exists():
        raise ValueError(f"does not exist: {folder}")
    probs = soda_import("deid").check_paths([folder])
    return {"clean": not probs, "count": len(probs), "lines": [scrub(p) for p in probs[:50]]}


def pick(kind: str) -> str:
    """A native folder / file picker, in its own process (a modal dialog must not run inside the
    server). macOS: AppleScript's chooser. Elsewhere: Tk. Empty string when cancelled."""
    if sys.platform == "darwin":
        what = "choose file" if kind == "file" else "choose folder"
        r = subprocess.run(["osascript", "-e", f'POSIX path of ({what} with prompt "SODA")'],
                           capture_output=True, text=True, timeout=600)
        return r.stdout.strip()
    r = subprocess.run([*self_cmd(), "--pick", kind], capture_output=True, text=True, timeout=600,
                       encoding="utf-8", errors="replace")
    return (r.stdout.strip().splitlines() or [""])[-1]


def run_picker(kind: str) -> int:
    import tkinter
    from tkinter import filedialog
    root = tkinter.Tk()
    root.withdraw()
    root.attributes("-topmost", True)
    if kind == "file":
        p = filedialog.askopenfilename(title="SODA: pick a file")
    else:
        p = filedialog.askdirectory(title="SODA: pick a folder", mustexist=True)
    root.destroy()
    print(p or "")
    return 0


# ---------------------------------------------------------------------------------------------
# the commands
# ---------------------------------------------------------------------------------------------
def cmd_home() -> int:
    wd = workdir()
    if sys.platform == "darwin" or sys.stdout is None:
        start_log(wd)
    version = app_version()
    preferred = int(os.environ.get("SODA_PORT", DEFAULT_PORT))
    running = ping(preferred)
    if running and running.get("version") == version:
        print(f"SODA {version} is already running: opening http://127.0.0.1:{preferred}/app/")
        open_browser(f"http://127.0.0.1:{preferred}/app/")
        return 0
    port = free_port(preferred)
    app = App(port, wd)
    httpd = app.serve()
    url = f"{app.base()}/app/"
    if not os.environ.get("SODA_NO_UPDATE_CHECK"):
        threading.Thread(target=lambda: setattr(app, "update", check_update(version)), daemon=True).start()
    print(f"SODA {version}")
    print(f"  running at {url}")
    print("  your browser opens on it. Keep this window open while you work; closing it stops SODA.")
    print(f"  logs: {wd}")
    threading.Timer(0.5, open_browser, args=(url,)).start()
    if sys.platform == "darwin":
        def idle():
            while True:
                time.sleep(60)
                if time.time() - app.last_seen > IDLE_EXIT_S and not (app.job and app.job.running):
                    print("no page open for 30 minutes: SODA stops")
                    httpd.shutdown()
                    return
        threading.Thread(target=idle, daemon=True).start()
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        if app.job and app.job.running:
            app.job.stop()
    return 0


def cmd_view(args: list[str]) -> int:
    import argparse
    ap = argparse.ArgumentParser(prog="SODA view", description="serve the viewer and an export folder")
    ap.add_argument("folder", help="the export folder (the one with datasets.json)")
    ap.add_argument("--port", type=int, default=int(os.environ.get("SODA_PORT", DEFAULT_PORT)))
    ap.add_argument("--no-browser", action="store_true")
    a = ap.parse_args(args)
    app = App(free_port(a.port), workdir())
    httpd = app.serve()
    try:
        url = app.viewer_url(a.folder)
    except ValueError as e:
        print(f"SODA view: {scrub(str(e))}", file=sys.stderr)
        return 2
    print(f"SODA {app.version}: viewer at {url}", flush=True)
    print("  Ctrl+C (or close this window) stops it.", flush=True)
    if not a.no_browser:
        threading.Timer(0.5, open_browser, args=(url,)).start()
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        pass
    return 0


def cmd_prepare(args: list[str]) -> int:
    import argparse
    ap = argparse.ArgumentParser(
        prog="SODA prepare",
        description="prepare.py batch for a folder of datasets, into <work>/export, <work>/cache and <work>/logs. "
                    "Any other prepare.py batch option (--options, --jobs, --key, --title, --only, --rebuild, "
                    "--dry-run, --no-ids ...) is passed on. The sample ID check is always a gate (--check-ids). Each export holds "
                    "the barcode / cell id lists by default (checked for sample IDs; --no-ids leaves them out).")
    ap.add_argument("--src", required=True, help="the folder of dataset folders")
    ap.add_argument("--work", required=True, help="the output folder")
    a, extra = ap.parse_known_args(args)
    try:
        argv = batch_argv(Path(a.src), Path(a.work), extra)
    except ValueError as e:
        print(f"SODA prepare: {e}", file=sys.stderr)
        return 2
    Path(a.work).expanduser().resolve().mkdir(parents=True, exist_ok=True)
    sys.argv = [str(SODA_DIR / "prepare.py"), *argv]
    rc = soda_import("batch").main(argv[1:])
    rep = read_report(Path(a.work).expanduser().resolve() / "logs", 0)
    if rep and rep.get("deid"):
        d = rep["deid"]
        print("")
        print("SODA: sample ID check " + ("CLEAN: the export folder can be shared" if d.get("clean") else
                                          "NOT CLEAN: do not upload or share the export folder yet"))
    return rc


def cmd_check(args: list[str]) -> int:
    if not args:
        print("SODA check <folder or file> ...", file=sys.stderr)
        return 2
    total = 0
    for a in args:
        r = check_folder(a)
        total += r["count"]
        for line in r["lines"]:
            print(line)
    print(f"[deid] {total} sample ID(s) found" if total else "[deid] clean")
    return 1 if total else 0


def cmd_selftest(args: list[str]) -> int:
    """Every library the tools use, imported and made to do one small real thing, plus every SODA
    module imported: a missing hidden import or DLL in the packaged app shows up here."""
    import tempfile
    quick = "--quick" in args
    t0 = time.time()
    checks: dict[str, str] = {}

    def t(name, fn):
        s = time.time()
        try:
            r = fn()
            checks[name] = f"ok ({time.time() - s:.1f}s)" + (f" {r}" if r else "")
        except BaseException as e:            # noqa: BLE001  (SystemExit from a module counts too)
            checks[name] = f"FAIL {type(e).__name__}: {e}"

    with tempfile.TemporaryDirectory(prefix="soda_selftest_") as td:
        tmp = Path(td)

        def numpy_():
            import numpy as np
            a = np.random.default_rng(0).random((64, 8))
            np.linalg.svd(a)
            return np.__version__

        def scipy_():
            import scipy.ndimage
            import scipy.sparse
            m = scipy.sparse.random(50, 50, density=0.1, format="csr", random_state=0)
            scipy.ndimage.gaussian_filter(m.toarray(), 1.0)
            import scipy
            return scipy.__version__

        def tiff_():
            import numpy as np
            import tifffile
            a = (np.arange(256 * 256, dtype=np.uint16) % 4096).reshape(256, 256)
            for comp in ("zlib", "jpeg2000", "lzw"):
                tifffile.imwrite(tmp / f"t_{comp}.tif", a if comp != "jpeg2000" else a.astype(np.uint8),
                                 compression=comp)
                tifffile.imread(tmp / f"t_{comp}.tif")
            return tifffile.__version__

        def codecs_():
            import imagecodecs
            import numpy as np
            a = np.zeros((32, 32), np.uint8)
            imagecodecs.jpeg8_decode(imagecodecs.jpeg8_encode(a))
            imagecodecs.png_decode(imagecodecs.png_encode(a))
            return imagecodecs.__version__

        def zarr_():
            import numcodecs
            import numpy as np
            import zarr
            z = zarr.open(str(tmp / "z.zarr"), mode="w", shape=(64, 64), chunks=(32, 32), dtype="u2",
                          compressor=numcodecs.Blosc(cname="zstd"))
            z[:] = 7
            assert int(zarr.open(str(tmp / "z.zarr"), mode="r")[3, 3]) == 7
            return f"zarr {zarr.__version__}, numcodecs {numcodecs.__version__}"

        def h5_():
            import h5py
            with h5py.File(tmp / "a.h5", "w") as f:
                f["x"] = [1, 2, 3]
            with h5py.File(tmp / "a.h5", "r") as f:
                assert list(f["x"][:]) == [1, 2, 3]
            return h5py.__version__

        def arrow_():
            import pyarrow as pa
            import pyarrow.csv
            import pyarrow.parquet as pq
            pq.write_table(pa.table({"a": [1, 2]}), tmp / "a.parquet")
            assert pq.read_table(tmp / "a.parquet").num_rows == 2
            return pa.__version__

        def pil_():
            from PIL import Image
            Image.new("RGB", (8, 8)).save(tmp / "a.png")
            Image.new("L", (8, 8)).save(tmp / "a.jpg")
            import PIL
            return PIL.__version__

        def crypto_():
            from cryptography.hazmat.primitives.ciphers.aead import AESGCM
            k = AESGCM.generate_key(bit_length=256)
            n = os.urandom(12)
            assert AESGCM(k).decrypt(n, AESGCM(k).encrypt(n, b"x", None), None) == b"x"

        def sklearn_():
            import numpy as np
            from sklearn.neighbors import NearestNeighbors
            NearestNeighbors(n_neighbors=3).fit(np.random.default_rng(0).random((50, 4)))
            import sklearn
            return sklearn.__version__

        def numba_():
            import numba
            import numpy as np

            @numba.njit(cache=False)
            def s(a):
                t = 0.0
                for x in a:
                    t += x
                return t
            assert s(np.ones(10)) == 10.0
            return numba.__version__

        def leiden_():
            import igraph as ig
            import leidenalg as la
            g = ig.Graph.Famous("Zachary")
            p = la.find_partition(g, la.RBConfigurationVertexPartition, seed=0)
            return f"{len(p)} communities"

        def umap_():
            import numpy as np
            sys.modules.setdefault("tensorflow", None)
            import umap
            e = umap.UMAP(n_neighbors=10, random_state=0, n_epochs=50).fit_transform(
                np.random.default_rng(0).random((200, 8)))
            assert e.shape == (200, 2)
            return umap.__version__ if hasattr(umap, "__version__") else ""

        def soda_modules():
            names = sorted(p.stem for p in SODA_DIR.glob("*.py"))
            bad = []
            for n in names:
                try:
                    soda_import(n)
                except BaseException as e:      # noqa: BLE001
                    bad.append(f"{n}: {type(e).__name__}: {e}")
            if bad:
                raise RuntimeError("; ".join(bad))
            return f"{len(names)} modules"

        for name, fn in (("numpy", numpy_), ("scipy", scipy_), ("tifffile", tiff_), ("imagecodecs", codecs_),
                         ("zarr", zarr_), ("h5py", h5_), ("pyarrow", arrow_), ("pillow", pil_),
                         ("cryptography", crypto_), ("scikit-learn", sklearn_), ("numba", numba_),
                         ("igraph+leidenalg", leiden_),
                         ("soda modules", soda_modules)):
            t(name, fn)
        if not quick:
            t("umap", umap_)

        def page_():
            missing = [n for n in ("index.html", "viewer.js", "styles.css") if not (SODA_DIR / n).is_file()]
            if missing:
                raise FileNotFoundError(", ".join(missing))
        t("page files", page_)
    ok = all(v.startswith("ok") for v in checks.values())
    print(json.dumps({"ok": ok, "app": APP_NAME, "version": app_version(), "frozen": FROZEN,
                      "python": sys.version.split()[0], "platform": sys.platform,
                      "seconds": round(time.time() - t0, 1), "checks": checks}, indent=1))
    return 0 if ok else 1


def resolve_script(name: str) -> Path | None:
    p = Path(name)
    if p.is_file():
        return p.resolve()
    q = SODA_DIR / p.name
    return q if q.is_file() else None


def run_script(script: Path, args: list[str]) -> int:
    """`SODA <script>.py args` = `python <script>.py args`: how batch.py's per-dataset processes,
    catalog.py and anything else started with sys.executable run inside the packaged app."""
    import runpy
    sys.argv = [str(script), *args]
    d = str(script.parent)
    if d not in sys.path:
        sys.path.insert(0, d)
    if str(SODA_DIR) not in sys.path:
        sys.path.insert(1, str(SODA_DIR))
    runpy.run_path(str(script), run_name="__main__")
    return 0


def setup_frozen_env():
    """Frozen only: numba's on-disk cache (umap-learn, pynndescent) needs a writable folder, and
    the install folder (Program Files, /Applications) may not be."""
    if FROZEN:
        os.environ.setdefault("NUMBA_CACHE_DIR", str(workdir() / "numba-cache"))


def main(argv: list[str] | None = None) -> int:
    if FROZEN:
        import multiprocessing
        multiprocessing.freeze_support()
    argv = list(sys.argv[1:] if argv is None else argv)
    argv = [a for a in argv if not a.startswith("-psn_")]      # old macOS Finder launch argument
    setup_streams()
    setup_frozen_env()
    if not argv:
        return cmd_home()
    cmd, rest = argv[0], argv[1:]
    if cmd == "--pick":
        return run_picker(rest[0] if rest else "folder")
    if cmd in ("version", "--version"):
        print(f"{APP_NAME} {app_version()}")
        return 0
    if cmd in ("-h", "--help", "help"):
        print(__doc__)
        return 0
    if cmd == "prepare":
        return cmd_prepare(rest)
    if cmd == "view":
        return cmd_view(rest)
    if cmd == "check":
        return cmd_check(rest)
    if cmd == "selftest":
        return cmd_selftest(rest)
    if cmd.lower().endswith(".py"):
        script = resolve_script(cmd)
        if script is None:
            print(f"SODA: no such script: {cmd}", file=sys.stderr)
            return 2
        return run_script(script, rest)
    print(f"SODA: unknown command {cmd!r}. Commands: prepare, view, check, selftest, version, <script>.py "
          "(SODA --help for more)", file=sys.stderr)
    return 2


if __name__ == "__main__":
    sys.exit(main())
