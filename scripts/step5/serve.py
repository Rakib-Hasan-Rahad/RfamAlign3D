"""Serve a results app locally with caching disabled (so a rebuild is always picked up).

  python3 scripts/step5/serve.py <app_dir> [port]
"""
import functools
import http.server
import sys


class NoCache(http.server.SimpleHTTPRequestHandler):
    def end_headers(self):
        self.send_header("Cache-Control", "no-store, max-age=0")
        self.send_header("Pragma", "no-cache")
        super().end_headers()


if __name__ == "__main__":
    d, port = sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 8770
    http.server.ThreadingHTTPServer(("127.0.0.1", port), functools.partial(NoCache, directory=d)).serve_forever()
