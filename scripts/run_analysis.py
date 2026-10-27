#!/usr/bin/env python3
"""Logged, non-overwriting entry point for existing Glide results (no docking submission)."""
import argparse
import datetime
import json
import subprocess
import sys
from pathlib import Path
from prepare_top_hits import file_sha256


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--n-clusters", type=int, default=100)
    parser.add_argument("--score-top", type=int, default=2000)
    parser.add_argument("--le-top", type=int, default=500)
    parser.add_argument("--max-per-cluster", type=int, default=1)
    parser.add_argument("--exclude-review-flags", action="store_true")
    parser.add_argument("--skip-density", action="store_true")
    parser.add_argument("--density-id-property", help="Reliable parent-ID property for the spatial analysis; default uses titles")
    parser.add_argument("--density-keep-all-poses", action="store_true", help="Use pose frequency, without title/ID deduplication")
    parser.add_argument("--number", type=int, default=10000, help="Number of score-sorted poses, not parents")
    parser.add_argument("--receptor-entry", type=int, help="Explicit unscored receptor entry in the main input")
    parser.add_argument("--with-ef", action="store_true", help="Also compute regional EF from the COMPLETE supplied EF cohort")
    parser.add_argument("--ef-input", type=Path, help="Full pose file for EF; default is --input, never the exported Top10k")
    parser.add_argument("--receptor", type=Path, help="Original-coordinate prepared receptor for EF if absent from that input")
    parser.add_argument("--chain", help="Protein chain for EF anchors")
    parser.add_argument("--region-definitions", type=Path, help="Reviewed four-region center/radius JSON")
    parser.add_argument("--ef-cohort-label", default="supplied input cohort; full-library coverage unverified")
    args = parser.parse_args()
    source, out = args.input.resolve(), args.out.resolve()
    if not source.is_file() or min(args.number, args.n_clusters, args.score_top, args.le_top, args.max_per_cluster) < 1:
        parser.error("Input must exist and counts must be positive")
    if not args.with_ef and any([args.ef_input, args.receptor, args.chain, args.region_definitions]):
        parser.error("EF options require --with-ef")
    ef_source = (args.ef_input or source).resolve()
    if args.with_ef and not ef_source.is_file():
        parser.error("EF input must exist")
    out.mkdir(parents=True, exist_ok=False)
    scripts = Path(__file__).resolve().parent
    prepared = out / "prepared"
    jobs = [
        ("prepare_top_hits.py", ["--input", str(source), "--out", str(prepared),
          "--number", str(args.number),
          *(["--receptor-entry", str(args.receptor_entry)] if args.receptor_entry else [])]),
        ("analyze_global_glide_1w.py", ["--input", str(prepared / "top1w_SP_pv.maegz"),
          "--output", str(out / "global"), *(["--skip-density"] if args.skip_density else []),
          *(["--id-property", args.density_id_property] if args.density_id_property else []),
          *(["--keep-all-poses"] if args.density_keep_all_poses else [])]),
        ("cluster_glide_top1w.py", ["--input", str(prepared / "top1w_SP_pv.maegz"),
          "--output-dir", str(out / "cluster"), "--n-clusters", str(args.n_clusters),
          "--score-top", str(args.score_top), "--le-top", str(args.le_top),
          "--max-molecules", str(args.number),
          "--max-per-cluster", str(args.max_per_cluster),
          *(["--exclude-review-flags"] if args.exclude_review_flags else [])]),
        ("draw_top100_2d.py", ["--input", str(prepared / "top100_SP_pv.maegz"),
          "--output", str(out / "top100_2D.png")]),
        ("plot_le_analysis.py", ["--csv", str(out / "cluster/cluster_assignments.csv"),
          "--out", str(out / "le")]),
    ]
    if args.with_ef:
        jobs.append(("calculate_pocket_ef.py", ["--input", str(ef_source),
            "--output-dir", str(out / "regional_ef"), "--cohort-label", args.ef_cohort_label,
            *(["--receptor", str(args.receptor.resolve())] if args.receptor else []),
            *(["--chain", args.chain] if args.chain else []),
            *(["--definitions", str(args.region_definitions.resolve())] if args.region_definitions else [])]))
    report = {"input_sha256": file_sha256(source), "python": sys.version,
              "started_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
              "analysis_parameters": {k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()},
              "ef_input_sha256": file_sha256(ef_source) if args.with_ef else None,
              "selection": "score -> LE -> cluster size; no weighted score", "steps": [], "status": "running"}
    manifest = out / "run_manifest.json"
    def save():
        manifest.write_text(json.dumps(report, indent=2), encoding="utf-8")
    save()
    try:
        for index, (name, flags) in enumerate(jobs, 1):
            command = [sys.executable, str(scripts / name), *flags]
            entry = {"script": name, "arguments": flags, "status": "running"}
            report["steps"].append(entry)
            save()
            print(f"[{index}/{len(jobs)}] {name}", flush=True)
            with (out / f"{index:02d}_{Path(name).stem}.log").open("w", encoding="utf-8") as log:
                result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT)
            entry["returncode"] = result.returncode
            entry["status"] = "complete" if result.returncode == 0 else "failed"
            save()
            if result.returncode:
                raise RuntimeError(f"{name} failed: see log in {out}")
        if report["input_sha256"] != file_sha256(source):
            raise RuntimeError("Source file changed during analysis")
        if args.with_ef and report["ef_input_sha256"] != file_sha256(ef_source):
            raise RuntimeError("EF source changed during analysis")
        report["status"] = "complete"
        report["artifacts"] = {
            str(p.relative_to(out)): file_sha256(p)
            for p in sorted(out.rglob("*")) if p.is_file() and p != manifest
        }
    except Exception as exc:
        report.update(status="failed", error=str(exc))
        raise
    finally:
        save()
    print("Review Top60 in Maestro. Historical manual ranks must use their ORIGINAL saved Top60.")


if __name__ == "__main__":
    main()
