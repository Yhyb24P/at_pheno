"""Small executable pilot. Large-scale ingestion and deep models are future stages."""

import argparse
import csv
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import time
import tomllib

import numpy as np

from .core import (additive_kernel, fit_qc, grouped_folds, marker_order,
                   predict, random_folds, scores, select)


def write_json(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False)+"\n")


def read_table(path, required):
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        if not set(required).issubset(reader.fieldnames or []):
            raise ValueError(f"{path}: required columns {required}")
        return list(reader)


def write_table(path, rows, fields):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024*1024), b""):
            h.update(chunk)
    return h.hexdigest()


def load_dataset(directory, trait):
    samples = read_table(directory/"samples.tsv", ["accession_id", "genetic_group"])
    variants = read_table(directory/"variants.tsv", ["chromosome", "position", "ref", "alt"])
    phenotypes = read_table(directory/"phenotypes.tsv", ["accession_id", "trait_id", "value"])
    ids = [r["accession_id"] for r in samples]
    if not ids or len(set(ids)) != len(ids) or any(not s or s != s.strip() for s in ids):
        raise ValueError("Sample IDs must be unique nonempty canonical strings")
    marker_ids, previous = [], None
    for row in variants:
        chrom, pos, ref, alt = (row[k] for k in ("chromosome", "position", "ref", "alt"))
        if chrom not in {"1", "2", "3", "4", "5"} or int(pos) < 1:
            raise ValueError("Require nuclear chromosome 1..5 and 1-based positions")
        if ref not in "ACGT" or alt not in "ACGT" or len(ref) != 1 or len(alt) != 1 or ref == alt:
            raise ValueError("Require biallelic SNP REF/ALT")
        coordinate = (int(chrom), int(pos), ref, alt)
        if previous is not None and coordinate <= previous:
            raise ValueError("variants.tsv must be strictly ordered by chromosome, position, REF, ALT")
        previous = coordinate
        marker_ids.append(f"{chrom}:{pos}:{ref}:{alt}")
    if not marker_ids or len(set(marker_ids)) != len(marker_ids):
        raise ValueError("Empty or duplicate variant manifest")
    x = np.load(directory/"genotypes.npy", mmap_mode="r", allow_pickle=False)
    if x.ndim != 2 or x.shape != (len(ids), len(marker_ids)) or x.dtype.kind not in "fiu":
        raise ValueError("Genotype shape/type does not match manifests")
    # Validate the entire matrix in bounded memory, including future test calls.
    for start in range(0, x.shape[1], 4096):
        g = x[:, start:start+4096]
        if np.any(~(np.isnan(g) | (g == 0) | (g == 1) | (g == 2))):
            raise ValueError("Genotypes must be diploid hard calls 0/1/2 or NaN")
    values, missing = {}, 0
    all_trait_ids = set()
    for row in phenotypes:
        if row["trait_id"] != trait:
            continue
        accession = row["accession_id"]
        if accession in all_trait_ids:
            raise ValueError("Repeated phenotype accession: aggregate under a documented protocol first")
        all_trait_ids.add(accession)
        if row["value"] in {"", "NA", "NaN", "nan"}:
            missing += 1
            continue
        value = float(row["value"])
        if not np.isfinite(value):
            raise ValueError("Nonfinite phenotype")
        values[accession] = value
    rows = np.array([i for i, accession in enumerate(ids) if accession in values], dtype=int)
    aligned_ids = [ids[i] for i in rows]
    y = np.array([values[s] for s in aligned_ids])
    groups = [samples[i]["genetic_group"] for i in rows]
    audit = {"trait_id": trait, "genotyped_n": len(ids), "marker_n": len(marker_ids),
             "phenotype_nonmissing_n": len(values), "phenotype_missing_n": missing,
             "intersection_n": len(rows), "phenotype_unmatched_ids": sorted(set(values)-set(ids)),
             "groups": {g: groups.count(g) for g in sorted(set(groups))}}
    if not len(rows):
        raise ValueError("No matched observations for this trait")
    return x, rows, aligned_ids, marker_ids, y, groups, audit


def formal_provenance(directory, trait):
    path = directory / "provenance.json"
    if not path.is_file():
        raise ValueError("Formal mode requires data/provenance.json")
    try:
        record = json.loads(path.read_text())
    except json.JSONDecodeError as exc:
        raise ValueError("Formal provenance is not valid JSON") from exc
    required = {"dataset_id", "genotype_representation", "reference_assembly",
                "source_urls", "input_sha256", "source_vcf_sha256", "conversion_script", "conversion_commit",
                "allele_encoding", "phenotype_registry"}
    missing = sorted(required-set(record))
    if missing:
        raise ValueError(f"Formal provenance missing fields: {', '.join(missing)}")
    if record["genotype_representation"] != "vcf_alt_dosage":
        raise ValueError("Formal CLI only accepts genotype_representation='vcf_alt_dosage'")
    if not isinstance(record["source_urls"], list) or not record["source_urls"]:
        raise ValueError("Formal provenance requires nonempty source_urls")
    if record["reference_assembly"] != "TAIR10":
        raise ValueError("Formal provenance requires reference_assembly='TAIR10'")
    if record["allele_encoding"] != "ALT_dosage_0_1_2":
        raise ValueError("Formal provenance requires ALT_dosage_0_1_2 encoding")
    if not all(isinstance(record[key], str) and len(record[key]) == 64
               for key in ("input_sha256", "source_vcf_sha256")):
        raise ValueError("Formal provenance requires SHA256 hashes")
    if record["input_sha256"] != sha256(directory / "genotypes.npy"):
        raise ValueError("Formal provenance input_sha256 does not bind genotypes.npy")
    registry = record["phenotype_registry"]
    if not isinstance(registry, dict) or registry.get("trait_id") != trait or registry.get("resolved") is not True:
        raise ValueError("Formal provenance requires resolved phenotype_registry for requested trait")
    return path, record


def save_qc(path, qc):
    np.savez_compressed(path, **vars(qc))


def run(args):
    config = tomllib.loads(args.config.read_text())
    mode = getattr(args, "mode", "pilot")
    formal_path, formal_record = (formal_provenance(args.data, args.trait) if mode == "formal" else (None, None))
    x, rows, ids, variants, y, groups, audit = load_dataset(args.data, args.trait)
    if len(y) < 12:
        raise ValueError("Pilot requires at least 12 matched accessions")
    if args.protocol == "group":
        fold = grouped_folds(groups)
    else:
        fold = random_folds(ids, config["outer_folds"], config["seed"])
    if args.splits:
        supplied = read_table(args.splits, ["accession_id", "fold"])
        mapping = {r["accession_id"]: int(r["fold"]) for r in supplied}
        if len(mapping) != len(supplied) or set(mapping) != set(ids):
            raise ValueError("Split manifest must cover matched accessions exactly once")
        fold = np.array([mapping[s] for s in ids])
        if len(set(fold)) < 2 or np.min(fold) < 0:
            raise ValueError("At least two nonnegative outer fold labels required")
        if args.protocol == "group":
            for group in set(groups):
                if len(set(fold[np.array(groups) == group])) != 1:
                    raise ValueError("A genetic group crosses held-out folds")
    args.out.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    order = marker_order(variants, config["hash_salt"])
    qc_args = {k: config[k] for k in ("min_call_rate", "min_maf", "min_mac", "block_size")}
    block_size = config["block_size"]
    alphas = config["alphas"]
    if not alphas or any(not np.isfinite(a) or a <= 0 for a in alphas):
        raise ValueError("Invalid alpha grid")
    predictions, tuning, omitted = [], [], []
    split_rows = [{"accession_id": s, "fold": int(f), "genetic_group": g}
                  for s, f, g in zip(ids, fold, groups)]
    write_table(args.out/"splits.tsv", split_rows, list(split_rows[0]))
    write_json(args.out/"audit.json", audit)
    write_json(args.out/"config.json", config)
    for outer in sorted(set(fold)):
        train, test = np.flatnonzero(fold != outer), np.flatnonzero(fold == outer)
        train_ids = [ids[i] for i in train]
        inner = (grouped_folds([groups[i] for i in train]) if args.protocol == "group"
                 else random_folds(train_ids, config["inner_folds"], config["seed"]+int(outer)+1))
        outer_qc = fit_qc(x, rows[train], **qc_args)
        save_qc(args.out/f"qc_outer_{outer}.npz", outer_qc)
        inner_parts = []
        for part in sorted(set(inner)):
            itrain, valid = train[inner != part], train[inner == part]
            iqc = fit_qc(x, rows[itrain], **qc_args)
            save_qc(args.out/f"qc_outer_{outer}_inner_{part}.npz", iqc)
            inner_parts.append((int(part), itrain, valid, iqc))
        write_table(args.out/f"inner_splits_{outer}.tsv",
                    [{"accession_id": s, "fold": int(f)} for s, f in zip(train_ids, inner)],
                    ["accession_id", "fold"])
        for density in config["densities"]:
            if density != "all" and len(outer_qc.columns) < density:
                omitted.append({"outer_fold": int(outer), "density": density,
                                "reason": "outer training universe smaller than requested density"})
                continue
            error = np.zeros(len(alphas))
            inner_counts = []
            for part, itrain, valid, iqc in inner_parts:
                chosen = select(iqc, order, density)
                inner_counts.append(len(chosen.columns))
                k, cross = additive_kernel(x, rows[itrain], rows[valid], chosen, block_size)
                for index, alpha in enumerate(alphas):
                    estimate = predict(k, cross, y[itrain], alpha)
                    error[index] += np.sum((y[valid]-estimate)**2)
            alpha = alphas[int(np.argmin(error))]
            chosen = select(outer_qc, order, density)
            save_qc(args.out/f"selected_outer_{outer}_{density}.npz", chosen)
            k, cross = additive_kernel(x, rows[train], rows[test], chosen, block_size)
            estimate = predict(k, cross, y[train], alpha)
            baseline = float(np.mean(y[train]))
            tuning.append({"outer_fold": int(outer), "density": density, "alpha": alpha,
                           "inner_mse": (error/len(train)).tolist(), "inner_marker_counts": inner_counts,
                           "marker_count": len(chosen.columns),
                           "inner_density_capped": density != "all" and any(c < density for c in inner_counts)})
            for i, p in zip(test, estimate):
                predictions.append({"accession_id": ids[i], "trait_id": args.trait,
                    "protocol": args.protocol, "fold": int(outer), "genetic_group": groups[i],
                    "density": str(density), "marker_count": len(chosen.columns),
                    "model": "additive_krr", "y": float(y[i]), "prediction": float(p),
                    "train_mean": baseline, "alpha": alpha})
    metrics = []
    for density in config["densities"]:
        records = [r for r in predictions if r["density"] == str(density)]
        if not records:
            continue
        for f in ["pooled"] + sorted(set(r["fold"] for r in records)):
            rs = records if f == "pooled" else [r for r in records if r["fold"] == f]
            truth, pred, base = [np.array([r[k] for r in rs]) for k in ("y", "prediction", "train_mean")]
            for model, pp in [("additive_krr", pred), ("train_mean", base)]:
                metrics.append({"density": str(density), "fold": f, "model": model,
                                "complete_oof": len(records) == len(ids), **scores(truth, pp, base)})
    if not predictions:
        write_json(args.out/"omitted.json", omitted)
        raise ValueError("No requested density was feasible; see omitted.json")
    write_table(args.out/"predictions.tsv", predictions, list(predictions[0]))
    write_json(args.out/"metrics.json", metrics)
    write_json(args.out/"tuning.json", tuning)
    write_json(args.out/"omitted.json", omitted)
    inputs = [args.data/name for name in ("genotypes.npy", "samples.tsv", "variants.tsv", "phenotypes.tsv")]
    inputs += [args.config] + ([args.splits] if args.splits else []) + ([formal_path] if formal_path else [])
    try:
        commit = subprocess.check_output(["git", "rev-parse", "HEAD"], stderr=subprocess.DEVNULL, text=True).strip()
    except subprocess.CalledProcessError:
        commit = None
    source_paths = sorted(Path(__file__).parent.glob("*.py"))
    write_json(args.out/"provenance.json", {"status": "formal_data_contract_passed_not_confirmatory" if mode == "formal" else "pilot_not_confirmatory", "mode": mode, "code_commit": commit,
        "python": platform.python_version(), "numpy": np.__version__, "elapsed_seconds": time.monotonic()-started,
        "source_sha256": {str(p): sha256(p) for p in source_paths},
        "input_sha256": {str(p): sha256(p) for p in inputs}, "formal_dataset_provenance": formal_record,
        "notes": ["Inner QC is refit, including frequency and eligibility.",
                  "IID protocol is not kinship controlled.",
                  "Pooled group PCC can be driven by between-group means; inspect fold metrics.",
                  "No significance tests or biological claims are produced."]})
    print(json.dumps({"output": str(args.out), "matched_n": len(ids), "prediction_rows": len(predictions)}))


def demo(directory):
    directory.mkdir(parents=True, exist_ok=False)
    rng = np.random.default_rng(20260922)
    n, p = 60, 160
    x = (2*rng.binomial(1, rng.uniform(0.08, 0.5, p), size=(n, p))).astype("float32")
    y = x[:, 5]*1.4 - x[:, 35]*0.9 + rng.normal(size=n)
    x[rng.random(x.shape) < 0.02] = np.nan
    np.save(directory/"genotypes.npy", x)
    write_table(directory/"samples.tsv", [{"accession_id": f"synthetic_{i:03}", "genetic_group": f"G{i%3}"}
                for i in range(n)], ["accession_id", "genetic_group"])
    write_table(directory/"variants.tsv", [{"chromosome": str(j//32+1), "position": str((j%32+1)*1000),
                "ref": "A", "alt": "G"} for j in range(p)], ["chromosome", "position", "ref", "alt"])
    write_table(directory/"phenotypes.tsv", [{"accession_id": f"synthetic_{i:03}", "trait_id": "synthetic",
                "value": float(y[i])} for i in range(n)], ["accession_id", "trait_id", "value"])
    write_json(directory/"provenance.json", {"synthetic": True, "seed": 20260922,
               "purpose": "Software smoke test; not Arabidopsis experimental evidence"})
    print(f"Synthetic fixture: {directory}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    d = sub.add_parser("demo")
    d.add_argument("--out", type=Path, required=True)
    a = sub.add_parser("audit")
    a.add_argument("--data", type=Path, required=True)
    a.add_argument("--trait", required=True)
    r = sub.add_parser("run")
    r.add_argument("--data", type=Path, required=True)
    r.add_argument("--trait", required=True)
    r.add_argument("--config", type=Path, default=Path("configs/pilot.toml"))
    r.add_argument("--out", type=Path, required=True)
    r.add_argument("--protocol", choices=["iid", "group"], default="iid")
    r.add_argument("--mode", choices=["pilot", "formal"], default="pilot")
    r.add_argument("--splits", type=Path, help="Frozen TSV accession_id/fold manifest")
    args = parser.parse_args()
    if args.command == "demo":
        demo(args.out)
    elif args.command == "audit":
        print(json.dumps(load_dataset(args.data, args.trait)[-1], ensure_ascii=False, indent=2))
    else:
        run(args)


if __name__ == "__main__":
    main()
