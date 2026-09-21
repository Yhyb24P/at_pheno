"""Extract only regular, expected HDF5/README members into a new directory."""

import argparse
import json
from pathlib import Path
import shutil
import tarfile


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("archive", type=Path)
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    records = []
    with tarfile.open(args.archive, "r:gz") as archive:
        for member in archive:
            record = {"name": member.name, "size": member.size, "extracted": False}
            records.append(record)
            name = Path(member.name).name
            if not member.isfile() or not (name.endswith(".hdf5") or name.lower().startswith("readme")):
                continue
            if member.size > 20_000_000_000:
                raise ValueError("Unexpectedly large archive member")
            # Flatten only explicit regular files; never trust archive directories/links.
            target = args.out/name
            with archive.extractfile(member) as src, target.open("xb") as dest:
                shutil.copyfileobj(src, dest, length=1024*1024)
            if target.stat().st_size != member.size:
                raise ValueError("Truncated archive member")
            record["extracted"] = True
            print(f"Extracted {name}: {member.size} bytes", flush=True)
    (args.out/"archive_members.json").write_text(json.dumps(records, indent=2)+"\n")


if __name__ == "__main__":
    main()
