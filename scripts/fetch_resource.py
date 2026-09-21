"""Fetch one public resource without overwriting, and record checksum provenance."""

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import time
import urllib.request


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("url")
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--max-bytes", type=int, default=400_000_000)
    args = parser.parse_args()
    partial = args.out.with_name(args.out.name+".part")
    manifest = args.out.with_name(args.out.name+".manifest.json")
    if any(p.exists() for p in (args.out, partial, manifest)):
        raise SystemExit("Destination/partial/manifest already exists; use a new output path")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    record = {"url": args.url, "retrieved_at": datetime.now(timezone.utc).isoformat(),
              "checksum_kind": "locally computed; no upstream checksum asserted"}
    digest, total, start = hashlib.sha256(), 0, time.monotonic()
    try:
        req = urllib.request.Request(args.url, headers={"User-Agent": "at-pheno-research/0.1"})
        with urllib.request.urlopen(req, timeout=45) as response, partial.open("xb") as stream:
            expected = response.headers.get("Content-Length")
            if expected and int(expected) > args.max_bytes:
                raise ValueError("Server content exceeds limit")
            record.update(resolved_url=response.url, content_length=expected,
                          last_modified=response.headers.get("Last-Modified"))
            while chunk := response.read(1024*1024):
                total += len(chunk)
                if total > args.max_bytes:
                    raise ValueError("Download exceeds limit")
                digest.update(chunk)
                stream.write(chunk)
            if expected and total != int(expected):
                raise ValueError("Incomplete response")
        partial.rename(args.out)
        record.update(status="ok", bytes=total, sha256=digest.hexdigest(), seconds=time.monotonic()-start)
    except (Exception, KeyboardInterrupt) as exc:
        record.update(status="failed", bytes=total, error=f"{type(exc).__name__}: {exc}")
        raise
    finally:
        manifest.write_text(json.dumps(record, indent=2)+"\n")
    print(json.dumps(record, indent=2))


if __name__ == "__main__":
    main()
