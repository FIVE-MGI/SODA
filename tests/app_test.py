#!/usr/bin/env python
r"""Test the SODA desktop app (app/soda_app.py) without opening any window.

    python tests/app_test.py                       the app from source (app/soda_app.py)
    python tests/app_test.py --exe <SODA.exe | SODA.app/Contents/MacOS/SODA>
                                                   a built app (app/build_release.ps1 / _mac.sh)
           [--work <folder>] [--full-selftest] [--browser] [--keep]

What it checks, through the app's command line and its local server only:
  1. `SODA version` and `SODA selftest` (every packaged library loads and does one real thing).
  2. The python stand-in batch.py relies on: `SODA deid.py text ...`, `SODA prepare.py detect ...`.
  3. `SODA prepare`: the README's demo data (examples/make_demo.py: a Visium HD folder and an
     OME-TIFF) as a folder of datasets -> <out>/export with datasets.json and both datasets, the
     batch log, and the sample ID check run as a gate (--check-ids) and clean. An output folder
     inside the data folder is refused.
  4. `SODA check <export>`: clean.
  5. `SODA view <export> --no-browser`: the viewer page, its scripts as JavaScript, an empty
     config.js, the export's files; the tools (.py) and paths outside the folder are not served,
     and a request with another Host header is refused.
  6. The home page (`SODA` with SODA_NO_BROWSER=1): served with its token; an action without the
     token is refused; Open viewer, Check and Prepare through its API, the Prepare result with the
     sample ID check; a second launch of the same version hands over to the first; Quit stops it.
  7. --browser: the viewer and the home page in headless Chromium (tests/app_browser_test.js;
     needs node and Playwright on NODE_PATH, and internet for deck.gl; skipped otherwise).

Every dataset name here is a made-up demo name. Exit code 0 = every check passed.
"""
from __future__ import annotations

import argparse
import http.client
import json
import os
import re
import shutil
import socket
import subprocess
import sys
import tempfile
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
DATASETS = ("demo_visiumhd", "demo_multiplex")


class Checker:
    def __init__(self):
        self.passed, self.failed = 0, []

    def ok(self, cond, what: str) -> bool:
        print(f"  {'ok  ' if cond else 'FAIL'} {what}", flush=True)
        if cond:
            self.passed += 1
        else:
            self.failed.append(what)
        return bool(cond)


def free_port() -> int:
    with socket.socket() as s:
        s.bind(("127.0.0.1", 0))
        return s.getsockname()[1]


def get(url: str, host: str | None = None, timeout=20, headers: dict | None = None, full=False):
    """(status, content-type, body) of a GET, with an optional Host header; full=True adds the
    response headers."""
    u = urllib.parse.urlsplit(url)
    c = http.client.HTTPConnection(u.hostname, u.port, timeout=timeout)
    path = u.path + (f"?{u.query}" if u.query else "")
    c.putrequest("GET", path, skip_host=True)
    c.putheader("Host", host or f"{u.hostname}:{u.port}")
    for k, v in (headers or {}).items():
        c.putheader(k, v)
    c.endheaders()
    r = c.getresponse()
    body = r.read()
    hdrs = {k.lower(): v for k, v in r.getheaders()}
    c.close()
    return (r.status, r.getheader("Content-Type") or "", body, hdrs) if full else (r.status, r.getheader("Content-Type") or "", body)


def post(base: str, path: str, body: dict, token: str | None):
    req = urllib.request.Request(base + path, data=json.dumps(body).encode(), method="POST",
                                 headers={"Content-Type": "application/json", **({"X-SODA-Token": token} if token else {})})
    try:
        with urllib.request.urlopen(req, timeout=60) as r:
            return r.status, json.loads(r.read().decode())
    except urllib.error.HTTPError as e:
        return e.code, json.loads(e.read().decode() or "{}")


def wait_ping(base: str, proc: subprocess.Popen, secs=120) -> bool:
    for _ in range(secs * 4):
        if proc.poll() is not None:
            return False
        try:
            with urllib.request.urlopen(base + "/app/api/ping", timeout=1) as r:
                return json.loads(r.read().decode()).get("app") == "SODA"
        except Exception:
            time.sleep(0.25)
    return False


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--exe", help="a built SODA (default: python app/soda_app.py)")
    ap.add_argument("--work", help="scratch folder (default: a new temporary folder)")
    ap.add_argument("--python", default=sys.executable, help="python with SODA's packages, for make_demo.py")
    ap.add_argument("--full-selftest", action="store_true", help="selftest with umap (slower)")
    ap.add_argument("--browser", action="store_true", help="also open the pages in headless Chromium")
    ap.add_argument("--keep", action="store_true", help="keep the scratch folder")
    a = ap.parse_args()

    soda = [str(Path(a.exe).resolve())] if a.exe else [sys.executable, str(REPO / "app" / "soda_app.py")]
    work = Path(a.work).resolve() if a.work else Path(tempfile.mkdtemp(prefix="soda_app_test_"))
    work.mkdir(parents=True, exist_ok=True)
    env = {**os.environ, "SODA_NO_BROWSER": "1", "SODA_NO_UPDATE_CHECK": "1", "SODA_WORKDIR": str(work / "appdata"),
           "PYTHONIOENCODING": "utf-8"}
    c = Checker()
    t_all = time.time()
    procs: list[subprocess.Popen] = []

    def run(args, timeout=1800):
        t0 = time.time()
        p = subprocess.run([*soda, *args], capture_output=True, text=True, encoding="utf-8", errors="replace",
                           env=env, timeout=timeout, cwd=str(work))
        with open(work / "app_test.log", "a", encoding="utf-8") as lg:
            lg.write(f"\n$ SODA {' '.join(args)}\n(exit {p.returncode}, {time.time() - t0:.1f}s)\n{p.stdout}\n{p.stderr}\n")
        return p, time.time() - t0

    print(f"SODA app test: {' '.join(soda)}\n  scratch: {work}", flush=True)
    try:
        # 1. version, selftest ----------------------------------------------------------------
        print("[1] version and selftest", flush=True)
        p, _ = run(["version"])
        version = p.stdout.strip().removeprefix("SODA ").strip()
        c.ok(p.returncode == 0 and p.stdout.startswith("SODA "), f"SODA version -> {p.stdout.strip()!r}")
        p, dt = run(["selftest"] + ([] if a.full_selftest else ["--quick"]))
        try:
            st = json.loads(p.stdout[p.stdout.index("{"):])
        except Exception:
            st = {}
        bad = [f"{k}: {v}" for k, v in (st.get("checks") or {}).items() if not v.startswith("ok")]
        c.ok(p.returncode == 0 and st.get("ok"), f"selftest ok in {dt:.1f}s ({len(st.get('checks') or {})} checks)"
             + (f"; failed: {'; '.join(bad)}" if bad else ""))

        # 2. python stand-in -------------------------------------------------------------------
        print("[2] SODA <script>.py (what batch.py's per-dataset processes use)", flush=True)
        p, _ = run(["deid.py", "text", "plain words only"])
        c.ok(p.returncode == 0 and p.stdout.strip() == "plain words only", "SODA deid.py text ... runs deid.py")
        p, _ = run(["no_such_script.py"])
        c.ok(p.returncode == 2, "an unknown script is refused (exit 2)")

        # 3. prepare ---------------------------------------------------------------------------
        print("[3] SODA prepare on the demo data", flush=True)
        raw, inputs = work / "raw", work / "inputs"
        if not (raw / "visiumhd" / "outs").is_dir():
            q = subprocess.run([a.python, str(REPO / "examples" / "make_demo.py"), "--out", str(raw)],
                               capture_output=True, text=True, encoding="utf-8", errors="replace", cwd=str(REPO))
            if not c.ok(q.returncode == 0, "examples/make_demo.py made the demo data"):
                print(q.stdout[-2000:], q.stderr[-2000:])
                return 1
        if inputs.exists():
            shutil.rmtree(inputs)
        shutil.copytree(raw / "visiumhd", inputs / "demo_visiumhd")
        (inputs / "demo_multiplex").mkdir()
        shutil.copy2(raw / "multiplex.ome.tif", inputs / "demo_multiplex" / "multiplex.ome.tif")
        p, _ = run(["prepare.py", "detect", str(inputs)])
        c.ok(p.returncode == 0 and p.stdout.count("CERTAIN") >= 2, "SODA prepare.py detect: both demo datasets CERTAIN")

        out = work / "out"
        if out.exists():
            shutil.rmtree(out)
        p, dt = run(["prepare", "--src", str(inputs), "--work", str(out)])
        export = out / "export"
        c.ok(p.returncode == 0, f"SODA prepare exit 0 ({dt:.0f}s)" + ("" if p.returncode == 0 else f": exit {p.returncode}, "
                                                                     + p.stdout[-800:].replace("\n", " | ")))
        cat = json.loads((export / "datasets.json").read_text(encoding="utf-8")) if (export / "datasets.json").exists() else {}
        ids = sorted(d.get("id") for d in cat.get("datasets", []))
        c.ok(ids == sorted(DATASETS), f"datasets.json lists {ids}")
        c.ok(all((export / d / "manifest.json").is_file() for d in DATASETS), "each dataset has its manifest.json")
        reps = sorted((out / "logs").glob("batch_*.json"))
        rep = json.loads(reps[-1].read_text(encoding="utf-8")) if reps else {}
        de = rep.get("deid") or {}
        c.ok(bool(reps) and (out / "logs" / reps[-1].with_suffix(".md").name).is_file(), "the batch log (md, json) is in <out>/logs")
        c.ok(de.get("gate") is True and de.get("clean") is True, f"sample ID check ran as a gate and is clean ({de.get('text')})")
        c.ok("--check-ids" in (rep.get("command") or ""), "the app passed --check-ids")
        c.ok("CLEAN: the export folder can be shared" in p.stdout, "SODA prepare says the check is clean")
        p, _ = run(["prepare", "--src", str(inputs), "--work", str(inputs / "out")])
        c.ok(p.returncode == 2 and not (inputs / "out").exists(), "an output folder inside the data folder is refused")

        # 4. check -----------------------------------------------------------------------------
        print("[4] SODA check", flush=True)
        p, _ = run(["check", str(export)])
        c.ok(p.returncode == 0 and "[deid] clean" in p.stdout, "SODA check <export>: clean")

        # 5. view ------------------------------------------------------------------------------
        print("[5] SODA view", flush=True)
        port = free_port()
        vp = subprocess.Popen([*soda, "view", str(export), "--no-browser", "--port", str(port)], stdout=subprocess.PIPE,
                              stderr=subprocess.STDOUT, text=True, encoding="utf-8", errors="replace", env=env)
        procs.append(vp)
        base = f"http://127.0.0.1:{port}"
        up = wait_ping(base, vp)
        line = vp.stdout.readline() if up else ""
        m = re.search(r"(http://127\.0\.0\.1:\d+/viewer/index\.html\?\S+)", line)
        if c.ok(up and m, f"SODA view serves on {base}"):
            url = m.group(1)
            c.ok("data=" in url and "/data/1/" in url, "the viewer address opens the folder's datasets.json (?data=)")
            s, ct, body = get(url)
            c.ok(s == 200 and ct.startswith("text/html") and b"viewer.js" in body, "the viewer page")
            s, ct, _ = get(base + "/viewer/viewer.js")
            c.ok(s == 200 and "javascript" in ct, f"viewer.js as JavaScript ({ct})")
            s, ct, body = get(base + "/viewer/config.js")
            c.ok(s == 200 and b"SpatialViewerConfig" in body, "an empty config.js")
            s, ct, body = get(base + "/data/1/datasets.json")
            c.ok(s == 200 and json.loads(body).get("datasets"), "datasets.json from the export")
            s, _, body, h = get(base + f"/data/1/{DATASETS[0]}/manifest.json", full=True)
            c.ok(s == 200, "a dataset's manifest.json")
            c.ok("access-control-allow-origin" not in h, "no CORS header (only the page on this address reads the export)")
            # HTTP Range, as blob storage: the viewer reads one gene out of a gene pack (genepack.py)
            size = len(body)
            s, _, part, h = get(base + f"/data/1/{DATASETS[0]}/manifest.json", headers={"Range": "bytes=5-14"}, full=True)
            c.ok(s == 206 and part == body[5:15] and h.get("content-range") == f"bytes 5-14/{size}",
                 f"Range request: 206 with Content-Range ({s}, {h.get('content-range')})")
            s, _, part, h = get(base + f"/data/1/{DATASETS[0]}/manifest.json", headers={"Range": "bytes=-4"}, full=True)
            c.ok(s == 206 and part == body[-4:], "suffix Range (bytes=-4)")
            s, _, _, h = get(base + f"/data/1/{DATASETS[0]}/manifest.json", headers={"Range": f"bytes={size + 10}-"}, full=True)
            c.ok(s == 416 and h.get("content-range") == f"bytes */{size}", "an unsatisfiable Range: 416")
            packs = sorted(q for q in export.rglob("*") if q.is_file() and q.suffix == ".pack" and q.stat().st_size > 64)
            if packs:
                rel = packs[0].relative_to(export).as_posix()
                raw_ = packs[0].read_bytes()
                s, _, part, h = get(base + f"/data/1/{rel}", headers={"Range": "bytes=16-47"}, full=True)
                c.ok(s == 206 and part == raw_[16:48], f"Range on a gene pack ({packs[0].name}): the right 32 bytes")
            s, _, _ = get(base + f"/data/1/{DATASETS[0]}")
            c.ok(s == 404, "a folder is not listed")
            for path in ("/viewer/prepare.py", "/viewer/deid.py", "/viewer/requirements-lock.txt",
                         "/data/1/../../out/logs/", "/data/1/%2e%2e/%2e%2e/raw/multiplex.ome.tif", "/data/2/datasets.json"):
                s, _, _ = get(base + path)
                c.ok(s == 404, f"not served: {path} ({s})")
            s, _, _ = get(base + "/viewer/index.html", host="evil.example:80")
            c.ok(s == 403, "another Host header is refused (DNS rebinding)")
            if a.browser:
                browser_test(c, url, base + "/app/", work, home=False)
        vp.terminate()

        # 6. the home page ---------------------------------------------------------------------
        print("[6] the home page (SODA with no arguments)", flush=True)
        port = free_port()
        henv = {**env, "SODA_PORT": str(port)}
        hp = subprocess.Popen(soda, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, encoding="utf-8",
                              errors="replace", env=henv)
        procs.append(hp)
        base = f"http://127.0.0.1:{port}"
        if c.ok(wait_ping(base, hp), f"SODA starts its home page on {base}"):
            s, ct, body = get(base + "/")
            c.ok(s == 302, "/ goes to /app/")
            s, ct, body = get(base + "/app/")
            m = re.search(rb'name="soda-token" content="([^"]+)"', body)
            token = m.group(1).decode() if m else None
            c.ok(s == 200 and token and b"{{" not in body and version.encode() in body,
                 "the home page, with its token and the version")
            s, d = post(base, "/app/api/check", {"folder": str(export)}, None)
            c.ok(s == 403, "an action without the token is refused")
            s, d = post(base, "/app/api/check", {"folder": str(export)}, "wrong")
            c.ok(s == 403, "an action with a wrong token is refused")
            s, d = post(base, "/app/api/view", {"folder": str(export)}, token)
            c.ok(s == 200 and "/viewer/index.html?data=" in d.get("url", ""), "Open viewer gives the viewer address")
            s, d = post(base, "/app/api/view", {"folder": str(inputs)}, token)
            c.ok(s == 400 and "no SODA export" in d.get("error", ""), "Open viewer on a folder that is not an export says so")
            s, d = post(base, "/app/api/check", {"folder": str(export)}, token)
            c.ok(s == 200 and d.get("clean") is True, "Check: the export is clean")
            s, d = post(base, "/app/api/prepare", {"src": str(inputs), "work": str(inputs / "x")}, token)
            c.ok(s == 400, "Prepare refuses an output folder inside the data folder")

            out2 = work / "out2"
            if out2.exists():
                shutil.rmtree(out2)
            t0 = time.time()
            s, d = post(base, "/app/api/prepare", {"src": str(inputs), "work": str(out2), "title": "App test"}, token)
            c.ok(s == 200 and d.get("ok"), "Prepare starts")
            s, d = post(base, "/app/api/prepare", {"src": str(inputs), "work": str(out2)}, token)
            c.ok(s == 400, "a second Prepare while one runs is refused")
            j, since, lines = {}, 0, []
            while time.time() - t0 < 1800:
                with urllib.request.urlopen(f"{base}/app/api/job?since={since}", timeout=10) as r:
                    j = json.loads(r.read().decode())
                lines += j.get("lines", [])
                since = j.get("next", since)
                if not j.get("running"):
                    break
                time.sleep(1)
            rep = j.get("report") or {}
            c.ok(j.get("exit") == 0, f"Prepare through the home page: exit {j.get('exit')} in {time.time() - t0:.0f}s, "
                                     f"{len(lines)} log lines streamed")
            c.ok((rep.get("deid") or {}).get("clean") is True and (rep.get("deid") or {}).get("gate") is True,
                 "the page gets the sample ID check result (gate, clean)")
            c.ok(sorted(r_["id"] for r_ in rep.get("datasets", []) if r_.get("status") == "ok") == sorted(DATASETS),
                 "the page gets the dataset table")
            c.ok(json.loads((out2 / "export" / "datasets.json").read_text(encoding="utf-8")).get("title") == "App test",
                 "the collection title reached datasets.json")
            plog = work / "appdata" / "prepare.log"
            c.ok(plog.is_file() and "[batch]" in plog.read_text(encoding="utf-8"), "the app's prepare.log is written")

            # a second launch of the same version hands over to the running one
            p2 = subprocess.run(soda, capture_output=True, text=True, encoding="utf-8", errors="replace", env=henv,
                                timeout=120)
            c.ok(p2.returncode == 0 and "already running" in p2.stdout, "a second launch hands over to the first")
            if a.browser:
                browser_test(c, "", base + "/app/", work, home=True)
            s, d = post(base, "/app/api/quit", {}, token)
            try:
                hp.wait(timeout=20)
            except subprocess.TimeoutExpired:
                pass
            c.ok(s == 200 and hp.poll() is not None, "Quit stops SODA")
        # the privacy rule for what the app writes: its logs and the batch logs. (Not appdata/numba-cache:
        # numba names its compiled-code files after the jitted function and a line number, and
        # pynndescent's tau_rand_int-<line> reads like an ID to deid.py; nothing there is ever shared.)
        logs_ = [str(q) for q in (work / "appdata").glob("*.log")] + [str(out / "logs"), str(out2 / "logs")]
        deid_hits = run(["check", *logs_])[0]
        c.ok(deid_hits.returncode == 0, "no sample ID in the app's own logs or the batch logs")
    finally:
        for pr in procs:
            if pr.poll() is None:
                pr.kill()
        print(f"\napp test: {c.passed} passed, {len(c.failed)} failed, in {time.time() - t_all:.0f}s", flush=True)
        for f in c.failed:
            print(f"  FAILED: {f}")
        if not a.keep and not a.work and not c.failed:
            shutil.rmtree(work, ignore_errors=True)
    return 0 if not c.failed else 1


def browser_test(c: Checker, viewer_url: str, home_url: str, work: Path, home: bool):
    node = shutil.which("node")
    if not node:
        print("  (node not found: browser check skipped)")
        return
    shots = work / "shots"
    shots.mkdir(exist_ok=True)
    args = [node, str(REPO / "tests" / "app_browser_test.js"), "--shots", str(shots)]
    args += ["--home", home_url] if home else ["--viewer", viewer_url, "--datasets", ",".join(DATASETS)]
    p = subprocess.run(args, capture_output=True, text=True, encoding="utf-8", errors="replace", timeout=300)
    print("    " + "\n    ".join(p.stdout.strip().splitlines()[-12:]))
    if "SKIPPED" in p.stdout:
        return
    c.ok(p.returncode == 0, f"headless Chromium: {'home page' if home else 'viewer'} (tests/app_browser_test.js)")


if __name__ == "__main__":
    sys.exit(main())
