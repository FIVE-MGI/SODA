#!/usr/bin/env python
r"""Serve an export folder over HTTP with CORS, for testing before upload.

    python static_serve.py <export root> 8765

then open the page with ?static=http://localhost:8765 - that path reads the files
directly and skips the password and the blob token entirely.

It answers HTTP Range requests (206 Partial Content), as blob storage does: the viewer reads
one gene out of an export's gene packs (genepack.py) with a Range request. Python's plain
`python -m http.server` ignores Range and would send the whole pack for every gene.

The README used to point at a copy of this that only existed in a scratch folder;
it lives here now so the local workflow keeps working.
"""
import os
import re
import sys
from functools import partial
from http import HTTPStatus
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer

_RANGE = re.compile(r"bytes=(\d*)-(\d*)")


class _Slice:
    """A file object that reads at most `n` bytes from where it stands (the body of a 206)."""

    def __init__(self, f, n):
        self.f, self.left = f, n

    def read(self, size=-1):
        if self.left <= 0:
            return b""
        size = self.left if size is None or size < 0 else min(size, self.left)
        b = self.f.read(size)
        self.left -= len(b)
        return b

    def close(self):
        self.f.close()


class CORSHandler(SimpleHTTPRequestHandler):
    """Static files with open CORS, no caching, and single-range Range requests."""

    def end_headers(self):
        self.send_header("Access-Control-Allow-Origin", "*")
        # the viewer checks a ranged answer against the shard size in Content-Range
        self.send_header("Access-Control-Expose-Headers", "Content-Range, Content-Length, Accept-Ranges")
        self.send_header("Accept-Ranges", "bytes")
        self.send_header("Cache-Control", "no-store")   # re-exports show up without a hard reload
        super().end_headers()

    def do_OPTIONS(self):
        # CORS preflight: a browser that does not count "Range: bytes=a-b" as a simple header asks first
        self.send_response(HTTPStatus.NO_CONTENT)
        self.send_header("Access-Control-Allow-Methods", "GET, HEAD, OPTIONS")
        self.send_header("Access-Control-Allow-Headers", self.headers.get("Access-Control-Request-Headers") or "Range")
        self.send_header("Access-Control-Max-Age", "3600")
        self.send_header("Content-Length", "0")
        self.end_headers()

    def send_head(self):
        rng = (self.headers.get("Range") or "").strip()
        m = _RANGE.fullmatch(rng)
        path = self.translate_path(self.path)
        if not m or not (m.group(1) or m.group(2)) or os.path.isdir(path):
            return super().send_head()     # no Range, or one this does not do (several ranges): the whole file
        try:
            f = open(path, "rb")
        except OSError:
            self.send_error(HTTPStatus.NOT_FOUND, "File not found")
            return None
        try:
            size = os.fstat(f.fileno()).st_size
            if m.group(1):
                start = int(m.group(1))
                end = min(int(m.group(2)) if m.group(2) else size - 1, size - 1)
            else:                                       # "bytes=-N": the last N bytes
                start, end = max(0, size - int(m.group(2))), size - 1
            if start >= size or start > end:
                f.close()
                self.send_response(HTTPStatus.REQUESTED_RANGE_NOT_SATISFIABLE)
                self.send_header("Content-Range", f"bytes */{size}")
                self.send_header("Content-Length", "0")
                self.end_headers()
                return None
            self.send_response(HTTPStatus.PARTIAL_CONTENT)
            self.send_header("Content-Type", self.guess_type(path))
            self.send_header("Content-Range", f"bytes {start}-{end}/{size}")
            self.send_header("Content-Length", str(end - start + 1))
            self.end_headers()
            f.seek(start)
            return _Slice(f, end - start + 1)
        except Exception:
            f.close()
            raise

    def log_message(self, fmt, *args):
        if "404" in (fmt % args):
            super().log_message(fmt, *args)


def serve(root, port, host="127.0.0.1"):
    """A threading server for `root` (port 0 = any free port); call .serve_forever() on it."""
    return ThreadingHTTPServer((host, port), partial(CORSHandler, directory=str(root)))


if __name__ == "__main__":
    root = sys.argv[1] if len(sys.argv) > 1 else "."
    port = int(sys.argv[2]) if len(sys.argv) > 2 else 8765
    print(f"serving {root} on http://localhost:{port}  (Ctrl+C to stop)")
    serve(root, port).serve_forever()
