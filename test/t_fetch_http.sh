#!/bin/bash
## A failed download says why. The mirror here is a scripted server rather
## than http.server, so it can answer 404, 429 with a Retry-After, or nothing
## at all, and the test reads the failure line for the status, the reason
## phrase and the seconds the server asked for. The transient case is the one
## that matters: a 429 that clears after two tries must end in a fetched,
## verified file with no user involvement, and a Retry-After longer than a
## fetch will wait must be reported rather than slept.
##
## Same fixture and manifest as t_fetch.sh; skips without python3 or libcurl.
set -euo pipefail
YAME=${YAME:?export YAME=/path/to/yame}
here=$(cd "$(dirname "$0")" && pwd)
root=$(cd "$here/.." && pwd)

command -v python3 >/dev/null || { echo "skip: no python3 for the mirror" >&2; exit 0; }
help=$("$YAME" fetch -h </dev/null 2>&1) || true
printf '%s\n' "$help" | grep -q 'fetch available' ||
  { echo "skip: built without libcurl" >&2; exit 0; }

d=$(mktemp -d); trap 'kill $srv 2>/dev/null || true; rm -rf "$d"' EXIT
cd "$d"

scope=EPIC/KYCG; tag=v8.1
sums="$root/tools/registry/sums/InfiniumAnnotation/$tag/$scope/SHA256SUMS"
tree="mirror/zhou-lab/InfiniumAnnotation/raw/$tag/$scope"
mkdir -p "$tree"
cp "$sums" "$tree/SHA256SUMS"
cp "$here/fixtures/Blacklist.20220304.cm" "$tree/"
asset="$scope/Blacklist.20220304.cm"

## The server's behaviour for the fixture is read from a file on every
## request, so one server covers every case: "ok", "404", "429:<n>:<secs>"
## (429 for the first n requests, then the bytes), "429:always:<secs>".
cat > srv.py <<'PY'
import os, sys
from http.server import BaseHTTPRequestHandler, HTTPServer
port, root = int(sys.argv[1]), sys.argv[2]
hits = {}
class H(BaseHTTPRequestHandler):
    def do_GET(s):
        mode = open("mode").read().strip()
        if s.path.endswith("Blacklist.20220304.cm"):
            n = hits[mode] = hits.get(mode, 0) + 1
            if mode == "404":
                s.send_response(404); s.end_headers(); return
            if mode.startswith("429:"):
                _, k, secs = mode.split(":")
                if k == "always" or n <= int(k):
                    s.send_response(429); s.send_header("Retry-After", secs)
                    s.end_headers(); return
        f = root + s.path
        if not os.path.isfile(f):
            s.send_response(404); s.end_headers(); return
        b = open(f, "rb").read()
        s.send_response(200); s.send_header("Content-Length", str(len(b)))
        s.end_headers(); s.wfile.write(b)
    def log_message(s, *a): pass
HTTPServer(("127.0.0.1", port), H).serve_forever()
PY
port=$(python3 -c 'import socket; s=socket.socket(); s.bind(("127.0.0.1",0)); print(s.getsockname()[1])')
echo ok > mode
python3 srv.py "$port" "$d/mirror" >server.log 2>&1 &
srv=$!
for i in $(seq 1 50); do
  python3 -c "import socket; socket.create_connection(('127.0.0.1', $port), timeout=1)" 2>/dev/null && break
  sleep 0.1
done
export YAME_ASSETS_MIRROR="http://127.0.0.1:$port"
export YAME_DATA_HOME="$d/store"
mkdir -p "$YAME_DATA_HOME"

## ---- 1. a 404 is named as one, with the URL that was asked ----------------
echo 404 > mode
rc=0; "$YAME" fetch -y "$asset" </dev/null > f404.log 2>&1 || rc=$?
[ "$rc" -ne 0 ] || { echo "a 404 fetch succeeded"; cat f404.log; exit 1; }
grep -q 'download failed: HTTP 404 Not Found: http://127.0.0.1:'"$port"'/.*Blacklist.20220304.cm' f404.log ||
  { echo "404 not reported as such:"; cat f404.log; exit 1; }
[ ! -e "$YAME_DATA_HOME/$asset" ] || { echo "a file landed after a 404"; exit 1; }

## ---- 2. a 429 that clears is retried, honouring Retry-After, and succeeds --
echo '429:2:1' > mode
t0=$(date +%s)
"$YAME" fetch -y "$asset" </dev/null > f429.log 2>&1 || { echo "fetch did not recover from a transient 429"; cat f429.log; exit 1; }
t1=$(date +%s)
cmp -s "$YAME_DATA_HOME/$asset" "$here/fixtures/Blacklist.20220304.cm" || { echo "recovered file differs from the fixture"; exit 1; }
grep -q 'HTTP 429 Too Many Requests (Retry-After: 1 s); retrying in 1 s (try 2 of 3)' f429.log ||
  { echo "first retry not announced with the server's Retry-After:"; cat f429.log; exit 1; }
grep -q 'retrying in 1 s (try 3 of 3)' f429.log || { echo "second retry not announced:"; cat f429.log; exit 1; }
## two 1 s waits: well under the 30 s default backoff, so Retry-After was honoured
[ $((t1 - t0)) -lt 20 ] || { echo "retry waited $((t1 - t0)) s: Retry-After ignored"; exit 1; }

## ---- 3. a 429 that never clears gives up after three tries, and says so ----
rm -f "$YAME_DATA_HOME/$asset"
echo '429:always:1' > mode
rc=0; "$YAME" fetch -y "$asset" </dev/null > f429b.log 2>&1 || rc=$?
[ "$rc" -ne 0 ] || { echo "a persistent 429 fetch succeeded"; cat f429b.log; exit 1; }
grep -q 'download failed after 3 tries: HTTP 429 Too Many Requests, the server is rate-limiting this client; wait a few minutes and re-run: http://' f429b.log ||
  { echo "persistent 429 not explained:"; cat f429b.log; exit 1; }
[ "$(grep -c 'retrying in' f429b.log)" -eq 2 ] || { echo "expected exactly 2 retry lines:"; cat f429b.log; exit 1; }

## ---- 4. a Retry-After longer than a fetch waits is reported, not slept -----
echo '429:always:600' > mode
t0=$(date +%s)
rc=0; "$YAME" fetch -y "$asset" </dev/null > f429c.log 2>&1 || rc=$?
t1=$(date +%s)
[ "$rc" -ne 0 ] || { echo "a 600 s Retry-After fetch succeeded"; exit 1; }
[ $((t1 - t0)) -lt 10 ] || { echo "slept on a 600 s Retry-After"; exit 1; }
grep -q 'download failed: HTTP 429 Too Many Requests, the server asks for a 600 s wait, longer than a fetch waits on its own (120 s); re-run after that: http://' f429c.log ||
  { echo "long Retry-After not reported:"; cat f429c.log; exit 1; }
! grep -q 'retrying in' f429c.log || { echo "retried despite a 600 s Retry-After"; cat f429c.log; exit 1; }

## ---- 5. below HTTP: a refused connection names the curl failure ------------
## The manifest is the first request, so that is the message that comes back.
kill $srv; wait $srv 2>/dev/null || true
rc=0; "$YAME" fetch -y "$asset" </dev/null > frefused.log 2>&1 || rc=$?
[ "$rc" -ne 0 ] || { echo "fetch succeeded with no server"; exit 1; }
grep -q "cannot fetch the manifest: curl: .*: http://127.0.0.1:$port/" frefused.log ||
  { echo "refused connection not explained:"; cat frefused.log; exit 1; }

echo "ok: t_fetch_http"
