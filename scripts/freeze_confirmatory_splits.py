"""Generate G0 split manifests without reading phenotype numeric values."""

import argparse
import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path


CANONICAL_SALT = "at_pheno_split_v1"
ALTERNATE_SALTS = (CANONICAL_SALT, "at_pheno_split_v1_s2", "at_pheno_split_v1_s3")
ROLE_USE = {"PRIMARY", "INDEPENDENT_REPLICATION", "SECONDARY_SENSITIVITY"}


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def read_tsv(path):
    with Path(path).open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def objective(folds, n_total, group_total, k):
    sizes = np_counts(folds.values())
    value = sum(((sizes.get(fold, 0) - n_total / k) / (n_total / k)) ** 2 for fold in range(k))
    group_term = 0.0
    for group, total in group_total.items():
        for fold in range(k):
            actual = sum(1 for _, target in folds.items() if target == fold and _[1] == group)
            group_term += ((actual - total / k) / max(1, total / k)) ** 2
    return value + group_term / max(1, len(group_total))


def np_counts(values):
    return Counter(values)


def block_order(blocks, salt):
    def key(members):
        payload = salt + "\0" + "\0".join(sorted(member[0] for member in members))
        return (-len(members), hashlib.sha256(payload.encode()).hexdigest())
    return sorted(blocks, key=key)


def assign_blocks(blocks, k, salt):
    """Greedy deterministic assignment using only IDs and published groups."""
    total = sum(len(block) for block in blocks)
    group_total = Counter(group for block in blocks for _, group, _ in block)
    assignments = {}
    for block in block_order(blocks, salt):
        candidates = []
        for fold in range(k):
            trial = dict(assignments)
            for accession, group, block_id in block:
                trial[(accession, group, block_id)] = fold
            candidates.append((objective(trial, total, group_total, k), fold, trial))
        _, _, assignments = min(candidates, key=lambda item: (item[0], item[1]))
    return {accession: fold for (accession, _, _), fold in assignments.items()}


def blocks_for_trait(ids, global_blocks, groups):
    present = set(ids)
    blocks = []
    for block_id, members in global_blocks.items():
        subset = [(accession, groups[accession], block_id) for accession in members if accession in present]
        if subset:
            blocks.append(subset)
    if sum(map(len, blocks)) != len(present):
        raise ValueError("trait IDs do not map exactly to global blocks")
    return blocks


def write_assignment(path, assignment, groups, block_lookup):
    rows = [{"accession_id": accession, "fold": assignment[accession], "genetic_group": groups[accession],
             "nearclone_block": block_lookup[accession]} for accession in sorted(assignment, key=lambda item: (assignment[item], item))]
    with Path(path).open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def read_trait_ids(path):
    """Read only accession IDs; numeric phenotype fields are never accessed."""
    records = json.loads(Path(path).read_text())
    ids = [str(record["accession_id"]) for record in records]
    if len(ids) != len(set(ids)):
        raise ValueError(f"{path}: duplicate accession IDs require resolution")
    return ids


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--trait-panel", type=Path, required=True)
    parser.add_argument("--blocks", type=Path, required=True)
    parser.add_argument("--samples", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    protected = [args.out / "blocks.tsv", args.out / "split_summary.json", args.out / "group_logo_manifest.tsv"]
    if any(path.exists() for path in protected):
        raise ValueError(f"refusing to overwrite existing split artifacts in {args.out}")
    sample_rows = read_tsv(args.samples)
    groups = {row["accession_id"]: row["genetic_group"] for row in sample_rows}
    if len(groups) != len(sample_rows) or any(not value for value in groups.values()):
        raise ValueError("sample manifest IDs/groups invalid")
    global_blocks = defaultdict(list)
    block_lookup = {}
    for row in read_tsv(args.blocks):
        accession, block = row["accession_id"], int(row["nearclone_block"])
        if accession not in groups or accession in block_lookup:
            raise ValueError("block manifest differs from samples or repeats accession")
        global_blocks[block].append(accession)
        block_lookup[accession] = block
    if set(block_lookup) != set(groups):
        raise ValueError("block manifest does not cover samples exactly")
    panel = [row for row in read_tsv(args.trait_panel) if row["role"] in ROLE_USE]
    if not panel:
        raise ValueError("no usable traits in frozen panel")
    args.out.mkdir(parents=True, exist_ok=True)
    # Canonical global block table is copied under its design-facing name.
    with (args.out / "blocks.tsv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=["accession_id", "nearclone_block"], delimiter="\t")
        writer.writeheader()
        writer.writerows({"accession_id": accession, "nearclone_block": block_lookup[accession]} for accession in sorted(groups))
    summaries, logo_rows = [], []
    for trait in panel:
        label = trait["trait_name"].replace(" ", "_")
        ids = read_trait_ids(trait["raw_values_path"])
        if set(ids) - set(groups):
            raise ValueError(f"{trait['trait_id']}: raw IDs absent from genotype samples")
        if len(ids) != int(trait["finite_panel_n"]):
            raise ValueError(f"{trait['trait_id']}: ID presence {len(ids)} != frozen finite_panel_n {trait['finite_panel_n']}")
        blocks = blocks_for_trait(ids, global_blocks, groups)
        for salt in ALTERNATE_SALTS:
            outer = assign_blocks(blocks, 5, salt + "\0" + trait["trait_id"])
            suffix = "" if salt == CANONICAL_SALT else "_" + salt.rsplit("_", 1)[-1]
            write_assignment(args.out / f"{label}_outer_blocked{suffix}.tsv", outer, groups, block_lookup)
            if salt == CANONICAL_SALT:
                for held_out in range(5):
                    train_ids = [accession for accession, fold in outer.items() if fold != held_out]
                    inner_blocks = blocks_for_trait(train_ids, global_blocks, groups)
                    inner = assign_blocks(inner_blocks, 4, salt + "\0" + trait["trait_id"] + f"\0outer{held_out}")
                    write_assignment(args.out / f"{label}_inner_blocked_outer{held_out}.tsv", inner, groups, block_lookup)
        singletons = [[(accession, groups[accession], block_lookup[accession])] for accession in ids]
        iid = assign_blocks(singletons, 5, CANONICAL_SALT + "\0iid\0" + trait["trait_id"])
        write_assignment(args.out / f"{label}_outer_iid.tsv", iid, groups, block_lookup)
        for held_out in range(5):
            train = [accession for accession, fold in iid.items() if fold != held_out]
            inner = assign_blocks([[(accession, groups[accession], block_lookup[accession])] for accession in train], 4,
                                  CANONICAL_SALT + f"\0iid\0{trait['trait_id']}\0outer{held_out}")
            write_assignment(args.out / f"{label}_inner_iid_outer{held_out}.tsv", inner, groups, block_lookup)
        for group in sorted({groups[accession] for accession in ids}):
            logo_rows.extend({"trait_id": trait["trait_id"], "trait_name": trait["trait_name"], "accession_id": accession,
                              "held_out_group": group, "logo_fold": group, "is_test": groups[accession] == group}
                             for accession in ids)
        summaries.append({"trait_id": trait["trait_id"], "trait_name": trait["trait_name"], "n": len(ids),
                          "groups": dict(sorted(Counter(groups[accession] for accession in ids).items())),
                          "global_blocks_present": len(blocks), "max_present_block": max(map(len, blocks))})
    with (args.out / "group_logo_manifest.tsv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(logo_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(logo_rows)
    summary = {"schema": "confirmatory_splits_v1", "status": "CONFIRMATORY_SPLITS_FROZEN", "canonical_salt": CANONICAL_SALT,
               "alternate_outer_salts": list(ALTERNATE_SALTS[1:]), "design": "IDs/presence/published-groups/exact-nearclone-blocks only; no phenotype numeric values",
               "inputs_sha256": {"trait_panel": digest(args.trait_panel), "blocks": digest(args.blocks), "samples": digest(args.samples)}, "traits": summaries}
    (args.out / "split_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({"status": summary["status"], "traits": [(item["trait_id"], item["n"]) for item in summaries]}))


if __name__ == "__main__":
    main()
