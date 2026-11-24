"""Launch the local classroom interface for the existing alignment implementation."""
import argparse
import base64
import json
import subprocess
import sys
import tempfile
import threading
import webbrowser
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

from alignment import smith_waterman

ROOT = Path(__file__).resolve().parent
PILOT_LOCK = threading.Lock()
DATASETS = {'BRCA1': ROOT / 'data', 'TP53': ROOT / 'data/TP53', 'BRCA2': ROOT / 'data/BRCA2'}


class DemoHandler(BaseHTTPRequestHandler):
    def send_json(self, value, status=200):
        body = json.dumps(value).encode()
        self.send_response(status)
        self.send_header('Content-Type', 'application/json')
        self.send_header('Content-Length', str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def do_GET(self):
        if self.path == '/api/datasets':
            self.send_json([name for name, folder in DATASETS.items() if (folder / 'manifest.json').is_file()])
            return
        files = {'/': ('web/index.html', 'text/html; charset=utf-8'),
                 '/app.js': ('web/app.js', 'text/javascript; charset=utf-8'),
                 '/style.css': ('web/style.css', 'text/css; charset=utf-8')}
        if self.path not in files:
            self.send_error(404)
            return
        path, content_type = files[self.path]
        body = (ROOT / path).read_bytes()
        self.send_response(200)
        self.send_header('Content-Type', content_type)
        self.send_header('Content-Length', str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def do_POST(self):
        # Only the page served by this local process can invoke its actions.
        if self.headers.get('Origin') != f'http://{self.headers.get("Host")}':
            self.send_json({'error': 'Open the demo using the local address printed in Terminal.'}, 403)
            return
        if self.path == '/api/toy':
            result = smith_waterman('TTACGTAA', 'GGACGTCC')
            self.send_json(result)
            return
        if self.path != '/api/pilot':
            self.send_error(404)
            return
        try:
            size = int(self.headers.get('Content-Length', '0'))
            if size > 4096:
                raise ValueError('Request too large')
            options = json.loads(self.rfile.read(size)) if size else {}
            gene = options.get('gene', 'BRCA1')
            folder = DATASETS[gene]
            manifest = json.loads((folder / 'manifest.json').read_text())
            window = min(600, min(s['cds_length'] for s in manifest['sequences']) // 3)
        except (ValueError, KeyError, OSError):
            self.send_json({'error': 'Choose an available gene dataset.'}, 400)
            return
        if not PILOT_LOCK.acquire(blocking=False):
            self.send_json({'error': 'A real-data run is already in progress. Please wait for it to finish.'}, 409)
            return
        try:
            # Each fresh process preserves the benchmark's memory measurement scope.
            with tempfile.TemporaryDirectory(prefix='brca1-classroom-') as folder:
                completed = subprocess.run(
                    [sys.executable, str(ROOT / 'demo.py'), '--output', folder, '--data-dir', str(DATASETS[gene]), '--window', str(window)],
                    cwd=ROOT, capture_output=True, text=True, timeout=180)
                if completed.returncode:
                    self.send_json({'error': 'The real-data run failed. In your activated environment, run python -m pip install -r requirements.txt, then restart the demo. Details are in Terminal.'}, 500)
                    print(completed.stderr, file=sys.stderr)
                    return
                output = Path(folder)
                payload = {name: json.loads((output / filename).read_text()) for name, filename in
                           [('alignments', 'alignments.json'), ('benchmark', 'benchmark.json')]}
                payload['figure'] = base64.b64encode((output / 'conservation.png').read_bytes()).decode()
                payload['manifest'] = manifest
                self.send_json(payload)
        except subprocess.TimeoutExpired:
            self.send_json({'error': 'This run exceeded three minutes. See Terminal, stop with Ctrl+C, and retry.'}, 504)
        finally:
            PILOT_LOCK.release()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--port', type=int, default=8000)
    parser.add_argument('--no-browser', action='store_true')
    args = parser.parse_args()
    try:
        server = ThreadingHTTPServer(('127.0.0.1', args.port), DemoHandler)
    except OSError as error:
        parser.exit(1, f'Could not start the local demo: {error}\nTry python classroom_demo.py --port 8001\n')
    address = f'http://127.0.0.1:{server.server_port}'
    print(f'\nBRCA1 classroom demo\nOpen {address}\nKeep this Terminal window open. Press Ctrl+C to stop.\n', flush=True)
    if not args.no_browser:
        webbrowser.open(address)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print('\nDemo stopped.')
    finally:
        server.server_close()


if __name__ == '__main__':
    main()
