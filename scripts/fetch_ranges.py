"""Bounded parallel HTTP byte ranges; validates each range before concatenation.

Completed chunks are retained for audit/resume. This is not upstream checksum verification.
"""

import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import time
import urllib.request


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("url")
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--workers", type=int, default=8)
    args = p.parse_args()
    if not 1 <= args.workers <= 8:
        raise ValueError("Use 1..8 workers")
    if args.out.exists():
        raise FileExistsError(args.out)
    request = urllib.request.Request(args.url, method="HEAD")
    with urllib.request.urlopen(request, timeout=45) as response:
        total = int(response.headers["Content-Length"])
        etag = response.headers.get("ETag")
    if total > 400_000_000:
        raise ValueError("Only bounded small resource downloads (<=400 MB)")
    folder = args.out.with_name(args.out.name+".ranges")
    folder.mkdir(parents=True, exist_ok=True)
    metadata = {"url": args.url, "bytes": total, "etag": etag, "chunk_size": 4*1024*1024}
    meta_path = folder/"source.json"
    if meta_path.exists() and json.loads(meta_path.read_text()) != metadata:
        raise ValueError("Remote resource changed; use a new destination")
    meta_path.write_text(json.dumps(metadata, indent=2)+"\n")
    start_time = time.monotonic()

    def fetch(start):
        end = min(start+metadata["chunk_size"], total)-1
        path = folder/f"{start:012}.bin"
        if path.exists() and path.stat().st_size == end-start+1:
            return path
        headers = {"Range": f"bytes={start}-{end}", "User-Agent": "at-pheno-research/0.1"}
        if etag:
            headers["If-Match"] = etag
        for attempt in range(3):
            try:
                req = urllib.request.Request(args.url, headers=headers)
                with urllib.request.urlopen(req, timeout=45) as response:
                    if response.status != 206 or response.headers.get("Content-Range") != f"bytes {start}-{end}/{total}":
                        raise ValueError("Server did not honor exact requested range")
                    payload = response.read(end-start+2)
                    if len(payload) != end-start+1:
                        raise ValueError("Incomplete or oversized range")
                path.write_bytes(payload)
                return path
            except Exception:
                if attempt == 2:
                    raise
        raise RuntimeError("Unreachable")

    positions = list(range(0, total, metadata["chunk_size"]))
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        paths = []
        for index, path in enumerate(executor.map(fetch, positions), 1):
            paths.append(path)
            print(f"Verified ranges {index}/{len(positions)}", flush=True)
    partial = args.out.with_name(args.out.name+".joining")
    digest = hashlib.sha256()
    with partial.open("xb") as stream:
        for path in paths:
            payload = path.read_bytes()
            digest.update(payload)
            stream.write(payload)
    partial.rename(args.out)
    metadata.update(sha256=digest.hexdigest(), completed_at=datetime.now(timezone.utc).isoformat(),
                    seconds=time.monotonic()-start_time, checksum_kind="locally computed", status="ok")
    args.out.with_name(args.out.name+".manifest.json").write_text(json.dumps(metadata, indent=2)+"\n")
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
