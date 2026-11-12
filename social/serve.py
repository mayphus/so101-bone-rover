"""Local render studio server: serves the repo and saves rendered PNGs to docs/social/.
usage: python3 social/serve.py   then open http://localhost:8791/social/studio.html"""
import http.server, os, re

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "docs", "social")


class Handler(http.server.SimpleHTTPRequestHandler):
    def __init__(self, *a, **k):
        super().__init__(*a, directory=ROOT, **k)

    def do_POST(self):
        m = re.fullmatch(r"/save/([a-z0-9_-]{1,40})\.png", self.path)
        if not m:
            self.send_error(400)
            return
        n = int(self.headers.get("Content-Length", 0))
        if n <= 0 or n > 40_000_000:
            self.send_error(413)
            return
        data = self.rfile.read(n)
        if not data.startswith(b"\x89PNG"):
            self.send_error(415)
            return
        os.makedirs(OUT, exist_ok=True)
        with open(os.path.join(OUT, m.group(1) + ".png"), "wb") as f:
            f.write(data)
        self.send_response(204)
        self.end_headers()

    def log_message(self, *a):
        pass


http.server.ThreadingHTTPServer(("127.0.0.1", 8791), Handler).serve_forever()
