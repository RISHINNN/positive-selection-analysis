#!/usr/bin/env python3
"""
Parallel PAML branch-site analysis with the correct mixture-distribution LRT
and Benjamini-Hochberg FDR correction.

The script is designed as a drop-in replacement for a typical step-3 runner:
  1. discover codon alignments (*.phy)
  2. create an isolated alt/null working directory for every gene
  3. run codeml
  4. parse log-likelihoods
  5. calculate branch-site LRT p-values using
         0.5 * chi-square(df=0) + 0.5 * chi-square(df=1)
  6. correct all valid gene-level p-values by BH-FDR

Python compatibility: 3.6+
External dependency: codeml (PAML)
"""

from __future__ import print_function

import argparse
import csv
import math
import os
import re
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


LNL_PATTERN = re.compile(
    r"lnL\s*\(\s*ntime\s*:\s*\d+\s+np\s*:\s*\d+\s*\)\s*:\s*"
    r"([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?)"
)

CONTROL_KEYS = ("seqfile", "treefile", "outfile")
SUMMARY_ALIASES = (
    "PAML_LRT_Results.tsv",
    "lrt_results.tsv",
    "paml_results_summary.tsv",
    "positive_selection_results.tsv",
)


def parse_lnl_text(text):
    """Return the last finite CODEML lnL value found in text."""
    values = []
    for match in LNL_PATTERN.finditer(text):
        try:
            value = float(match.group(1))
        except (TypeError, ValueError):
            continue
        if math.isfinite(value):
            values.append(value)
    if not values:
        raise ValueError("No finite CODEML lnL line was found")
    return values[-1]


def parse_lnl_file(path):
    path = Path(path)
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError("Missing or empty CODEML result: {0}".format(path))
    return parse_lnl_text(path.read_text(encoding="utf-8", errors="replace"))


def branch_site_lrt(lnl_alt, lnl_null, tolerance=1e-6):
    """
    Calculate the branch-site Model A likelihood-ratio test.

    For LRT > 0, the null distribution is a 50:50 mixture of a point mass
    at zero and chi-square(df=1).  Because chi-square(df=1) survival is
    erfc(sqrt(LRT / 2)), the mixture p-value is:

        p = 0.5 * erfc(sqrt(LRT / 2))

    Returns
    -------
    (lrt, p_value, status)
    """
    try:
        lnl_alt = float(lnl_alt)
        lnl_null = float(lnl_null)
    except (TypeError, ValueError):
        return float("nan"), float("nan"), "non_numeric_lnl"

    if not math.isfinite(lnl_alt) or not math.isfinite(lnl_null):
        return float("nan"), float("nan"), "non_finite_lnl"

    delta = lnl_alt - lnl_null

    # Tiny negative values are ordinary numerical noise.
    if -tolerance <= delta < 0.0:
        delta = 0.0

    # The alternative model contains the null model. A clearly lower alt lnL
    # therefore indicates failed optimisation, a bad parse, or stale output.
    if delta < -tolerance:
        return float("nan"), float("nan"), "alt_lnl_lower_than_null"

    lrt = 2.0 * delta
    if lrt <= 0.0:
        return 0.0, 1.0, "ok"

    p_value = 0.5 * math.erfc(math.sqrt(lrt / 2.0))
    p_value = min(1.0, max(0.0, p_value))
    return lrt, p_value, "ok"


def benjamini_hochberg(p_values):
    """
    Benjamini-Hochberg FDR adjustment preserving input order.

    Non-finite/invalid values are returned as NaN and are excluded from the
    number of tests. This is intentional: failed CODEML jobs are not tests.
    """
    adjusted = [float("nan")] * len(p_values)
    valid = []

    for index, value in enumerate(p_values):
        try:
            p = float(value)
        except (TypeError, ValueError):
            continue
        if math.isfinite(p) and 0.0 <= p <= 1.0:
            valid.append((index, p))

    m = len(valid)
    if m == 0:
        return adjusted

    valid.sort(key=lambda item: item[1])
    running_min = 1.0

    for rank_index in range(m - 1, -1, -1):
        original_index, p = valid[rank_index]
        rank = rank_index + 1
        candidate = p * m / float(rank)
        running_min = min(running_min, candidate, 1.0)
        adjusted[original_index] = running_min

    return adjusted


def render_ctl(template_path, output_path, seqfile, treefile, outfile):
    """Copy a CODEML control template and replace its three file paths."""
    template_path = Path(template_path)
    output_path = Path(output_path)

    if not template_path.is_file():
        raise FileNotFoundError("Control template not found: {0}".format(template_path))

    replacements = {
        "seqfile": str(seqfile),
        "treefile": str(treefile),
        "outfile": str(outfile),
    }
    found = set()
    rendered = []

    for raw_line in template_path.read_text(
        encoding="utf-8", errors="replace"
    ).splitlines():
        match = re.match(r"^(\s*)(seqfile|treefile|outfile)(\s*)=.*$", raw_line)
        if match:
            key = match.group(2)
            rendered.append("{0}{1} = {2}".format(match.group(1), key, replacements[key]))
            found.add(key)
        else:
            rendered.append(raw_line)

    for key in CONTROL_KEYS:
        if key not in found:
            rendered.append("{0} = {1}".format(key, replacements[key]))

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(rendered) + "\n", encoding="utf-8")


def _is_integer_string(value):
    try:
        int(value)
        return True
    except (TypeError, ValueError):
        return False


def normalize_legacy_argv(argv):
    """
    Accept both common historical uses of '-t':
      -t 20             -> workers
      -t labeled.nwk    -> tree
    """
    normalized = []
    index = 0
    while index < len(argv):
        token = argv[index]
        if token == "-t" and index + 1 < len(argv):
            value = argv[index + 1]
            if _is_integer_string(value):
                normalized.extend(["--workers", value])
            else:
                normalized.extend(["--tree", value])
            index += 2
            continue
        normalized.append(token)
        index += 1
    return normalized


def first_existing(candidates, want_directory=False):
    for candidate in candidates:
        path = Path(candidate)
        if want_directory and path.is_dir():
            return path
        if not want_directory and path.is_file():
            return path
    return None


def resolve_executable(command):
    command_path = Path(command)
    if command_path.parent != Path(".") or os.path.sep in command:
        if command_path.is_file() and os.access(str(command_path), os.X_OK):
            return str(command_path.resolve())
        raise FileNotFoundError("CODEML executable is not executable: {0}".format(command))

    resolved = shutil.which(command)
    if resolved is None:
        raise FileNotFoundError(
            "Cannot find '{0}' in PATH. Install PAML or pass --codeml PATH.".format(command)
        )
    return resolved


def safe_gene_name(alignment_path, alignment_root):
    relative = alignment_path.relative_to(alignment_root)
    without_suffix = relative.with_suffix("")
    return "__".join(without_suffix.parts)


def discover_alignments(alignment_dir, pattern):
    alignment_dir = Path(alignment_dir)
    paths = sorted(path for path in alignment_dir.rglob(pattern) if path.is_file())
    if not paths:
        raise FileNotFoundError(
            "No alignment files matching '{0}' found under {1}".format(
                pattern, alignment_dir
            )
        )

    names = {}
    for path in paths:
        gene = safe_gene_name(path, alignment_dir)
        if gene in names:
            raise ValueError(
                "Duplicate output gene name '{0}' from {1} and {2}".format(
                    gene, names[gene], path
                )
            )
        names[gene] = path
    return [(gene, path) for gene, path in sorted(names.items())]


def result_is_reusable(result_path, dependencies):
    result_path = Path(result_path)
    if not result_path.is_file() or result_path.stat().st_size == 0:
        return False
    result_mtime = result_path.stat().st_mtime
    for dependency in dependencies:
        if Path(dependency).stat().st_mtime > result_mtime:
            return False
    try:
        parse_lnl_file(result_path)
    except ValueError:
        return False
    return True


def run_one_model(
    model_name,
    work_dir,
    template_ctl,
    alignment_source,
    tree_source,
    codeml_executable,
    force=False,
):
    work_dir = Path(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)

    alignment_local = work_dir / "alignment.phy"
    tree_local = work_dir / "tree.nwk"
    ctl_local = work_dir / "codeml.ctl"
    result_local = work_dir / "result.txt"
    log_local = work_dir / "codeml.log"

    shutil.copy2(str(alignment_source), str(alignment_local))
    shutil.copy2(str(tree_source), str(tree_local))
    render_ctl(template_ctl, ctl_local, alignment_local.name, tree_local.name, result_local.name)

    dependencies = (alignment_local, tree_local, ctl_local)
    if not force and result_is_reusable(result_local, dependencies):
        return {
            "model": model_name,
            "lnl": parse_lnl_file(result_local),
            "status": "reused",
            "error": "",
        }

    if result_local.exists():
        result_local.unlink()

    with log_local.open("w", encoding="utf-8") as log_handle:
        completed = subprocess.run(
            [codeml_executable, ctl_local.name],
            cwd=str(work_dir),
            stdout=log_handle,
            stderr=subprocess.STDOUT,
            universal_newlines=True,
        )

    if completed.returncode != 0:
        return {
            "model": model_name,
            "lnl": float("nan"),
            "status": "codeml_failed",
            "error": "codeml exit code {0}; see {1}".format(
                completed.returncode, log_local
            ),
        }

    try:
        lnl = parse_lnl_file(result_local)
    except ValueError as exc:
        return {
            "model": model_name,
            "lnl": float("nan"),
            "status": "parse_failed",
            "error": "{0}; see {1}".format(exc, log_local),
        }

    return {"model": model_name, "lnl": lnl, "status": "completed", "error": ""}


def run_gene_job(job):
    gene = job["gene"]
    alignment = Path(job["alignment"])
    gene_root = Path(job["output_dir"]) / gene

    try:
        alt = run_one_model(
            "alt",
            gene_root / "alt",
            job["alt_ctl"],
            alignment,
            job["tree"],
            job["codeml"],
            force=job["force"],
        )
        null = run_one_model(
            "null",
            gene_root / "null",
            job["null_ctl"],
            alignment,
            job["tree"],
            job["codeml"],
            force=job["force"],
        )
    except Exception as exc:
        return {
            "Gene": gene,
            "Alignment": str(alignment),
            "lnL_Alt": float("nan"),
            "lnL_Null": float("nan"),
            "LRT": float("nan"),
            "P_value": float("nan"),
            "Status": "job_failed",
            "Error": str(exc),
        }

    model_errors = [item["error"] for item in (alt, null) if item["error"]]
    if alt["status"] in ("codeml_failed", "parse_failed") or null["status"] in (
        "codeml_failed",
        "parse_failed",
    ):
        return {
            "Gene": gene,
            "Alignment": str(alignment),
            "lnL_Alt": alt["lnl"],
            "lnL_Null": null["lnl"],
            "LRT": float("nan"),
            "P_value": float("nan"),
            "Status": "model_failed",
            "Error": " | ".join(model_errors),
        }

    lrt, p_value, lrt_status = branch_site_lrt(alt["lnl"], null["lnl"])
    status = lrt_status
    if lrt_status == "ok":
        if alt["status"] == "reused" and null["status"] == "reused":
            status = "ok_reused"
        else:
            status = "ok"

    return {
        "Gene": gene,
        "Alignment": str(alignment),
        "lnL_Alt": alt["lnl"],
        "lnL_Null": null["lnl"],
        "LRT": lrt,
        "P_value": p_value,
        "Status": status,
        "Error": " | ".join(model_errors),
    }


def format_number(value):
    try:
        number = float(value)
    except (TypeError, ValueError):
        return "NA"
    if not math.isfinite(number):
        return "NA"
    return "{0:.12g}".format(number)


def write_summary(rows, output_path, alpha):
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = [
        "Gene",
        "lnL_Alt",
        "lnL_Null",
        "LRT",
        "P_value",
        # Kept for compatibility; now intentionally based on BH-FDR.
        "Significant",
        "BH_FDR",
        "q_value",
        "Significant_BH_FDR",
        "Significant_raw_P",
        "Status",
        "Error",
        "Alignment",
    ]

    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            q_value = row.get("BH_FDR", float("nan"))
            p_value = row.get("P_value", float("nan"))
            fdr_significant = (
                math.isfinite(q_value) and q_value <= alpha and row["Status"].startswith("ok")
            )
            raw_significant = (
                math.isfinite(p_value) and p_value <= alpha and row["Status"].startswith("ok")
            )
            writer.writerow(
                {
                    "Gene": row["Gene"],
                    "lnL_Alt": format_number(row["lnL_Alt"]),
                    "lnL_Null": format_number(row["lnL_Null"]),
                    "LRT": format_number(row["LRT"]),
                    "P_value": format_number(row["P_value"]),
                    "Significant": "Yes" if fdr_significant else "No",
                    "BH_FDR": format_number(q_value),
                    "q_value": format_number(q_value),
                    "Significant_BH_FDR": "Yes" if fdr_significant else "No",
                    "Significant_raw_P": "Yes" if raw_significant else "No",
                    "Status": row["Status"],
                    "Error": row.get("Error", ""),
                    "Alignment": row["Alignment"],
                }
            )


def create_summary_aliases(primary_path):
    primary_path = Path(primary_path)
    for alias_name in SUMMARY_ALIASES:
        alias_path = primary_path.parent / alias_name
        if alias_path.resolve() == primary_path.resolve():
            continue
        shutil.copy2(str(primary_path), str(alias_path))


def build_parser():
    parser = argparse.ArgumentParser(
        description=(
            "Run PAML branch-site alternative/null models in parallel, calculate "
            "the correct 50:50 mixture LRT p-value, and add BH-FDR."
        )
    )
    parser.add_argument(
        "-i",
        "--input",
        "--input-dir",
        "--alignment-dir",
        "--phy-dir",
        dest="alignment_dir",
        help="Directory containing codon alignments (*.phy)",
    )
    parser.add_argument("--tree", help="Foreground-labelled Newick tree")
    parser.add_argument(
        "-a", "--alt-ctl", "--alt_ctl", dest="alt_ctl", help="Alternative-model CTL"
    )
    parser.add_argument(
        "-n", "--null-ctl", "--null_ctl", dest="null_ctl", help="Null-model CTL"
    )
    parser.add_argument(
        "-o",
        "--output",
        "--output-dir",
        dest="output_dir",
        default="PAML_results",
        help="Output directory [default: PAML_results]",
    )
    parser.add_argument(
        "-j",
        "-p",
        "--threads",
        "--workers",
        dest="workers",
        type=int,
        default=max(1, min(8, os.cpu_count() or 1)),
        help="Concurrent gene jobs",
    )
    parser.add_argument("--codeml", default="codeml", help="CODEML executable")
    parser.add_argument("--pattern", default="*.phy", help="Alignment glob pattern")
    parser.add_argument(
        "--summary",
        default=None,
        help="Primary summary TSV [default: OUTPUT/paml_lrt_results.tsv]",
    )
    parser.add_argument(
        "--alpha", type=float, default=0.05, help="FDR significance cutoff [0.05]"
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Rerun CODEML even when a fresh parseable result exists",
    )
    parser.add_argument(
        "legacy_args",
        nargs="*",
        help=argparse.SUPPRESS,
    )
    return parser


def apply_positional_compatibility(args, parser):
    values = list(args.legacy_args)
    targets = ("alignment_dir", "tree", "alt_ctl", "null_ctl", "output_dir")
    for target in targets:
        if not values:
            break
        current = getattr(args, target)
        # output_dir has a default, but a fifth positional argument should override it.
        if current is None or target == "output_dir":
            setattr(args, target, values.pop(0))
    if values:
        if len(values) == 1 and _is_integer_string(values[0]):
            args.workers = int(values[0])
        else:
            parser.error("Unrecognized positional arguments: {0}".format(" ".join(values)))
    return args


def resolve_default_inputs(args, parser):
    if args.alignment_dir is None:
        found = first_existing(
            ("Codon_Alignments", "codon_alignments", "CodonAlignments", "alignments"),
            want_directory=True,
        )
        if found is not None:
            args.alignment_dir = str(found)

    if args.tree is None:
        found = first_existing(
            (
                "labeled_tree.nwk",
                "foreground_labeled_tree.nwk",
                "SpeciesTree_rooted_labeled.txt",
                "SpeciesTree_rooted.txt",
            )
        )
        if found is not None:
            args.tree = str(found)

    if args.alt_ctl is None:
        found = first_existing(("alt_model.ctl", "alternative_model.ctl", "alt.ctl"))
        if found is not None:
            args.alt_ctl = str(found)

    if args.null_ctl is None:
        found = first_existing(("null_model.ctl", "null.ctl"))
        if found is not None:
            args.null_ctl = str(found)

    missing = []
    for name in ("alignment_dir", "tree", "alt_ctl", "null_ctl"):
        if getattr(args, name) is None:
            missing.append(name)
    if missing:
        parser.error(
            "Missing required inputs: {0}. Run with --help for accepted options.".format(
                ", ".join(missing)
            )
        )

    args.alignment_dir = str(Path(args.alignment_dir).resolve())
    args.tree = str(Path(args.tree).resolve())
    args.alt_ctl = str(Path(args.alt_ctl).resolve())
    args.null_ctl = str(Path(args.null_ctl).resolve())
    args.output_dir = str(Path(args.output_dir).resolve())

    if not Path(args.alignment_dir).is_dir():
        parser.error("Alignment directory does not exist: {0}".format(args.alignment_dir))
    for label, value in (
        ("Tree", args.tree),
        ("Alternative CTL", args.alt_ctl),
        ("Null CTL", args.null_ctl),
    ):
        if not Path(value).is_file():
            parser.error("{0} does not exist: {1}".format(label, value))

    if args.workers < 1:
        parser.error("--workers must be >= 1")
    if not 0.0 < args.alpha < 1.0:
        parser.error("--alpha must be between 0 and 1")

    return args


def main(argv=None):
    raw_argv = sys.argv[1:] if argv is None else list(argv)
    normalized_argv = normalize_legacy_argv(raw_argv)
    parser = build_parser()
    args = parser.parse_args(normalized_argv)
    args = apply_positional_compatibility(args, parser)
    args = resolve_default_inputs(args, parser)

    codeml_executable = resolve_executable(args.codeml)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    alignments = discover_alignments(args.alignment_dir, args.pattern)
    print("[INFO] Found {0} codon alignments".format(len(alignments)))
    print("[INFO] Running up to {0} genes concurrently".format(args.workers))

    jobs = []
    for gene, alignment in alignments:
        jobs.append(
            {
                "gene": gene,
                "alignment": str(alignment),
                "output_dir": str(output_dir),
                "tree": args.tree,
                "alt_ctl": args.alt_ctl,
                "null_ctl": args.null_ctl,
                "codeml": codeml_executable,
                "force": args.force,
            }
        )

    rows = []
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        future_to_gene = {executor.submit(run_gene_job, job): job["gene"] for job in jobs}
        completed_count = 0
        for future in as_completed(future_to_gene):
            gene = future_to_gene[future]
            try:
                row = future.result()
            except Exception as exc:
                row = {
                    "Gene": gene,
                    "Alignment": "",
                    "lnL_Alt": float("nan"),
                    "lnL_Null": float("nan"),
                    "LRT": float("nan"),
                    "P_value": float("nan"),
                    "Status": "worker_failed",
                    "Error": str(exc),
                }
            rows.append(row)
            completed_count += 1
            print(
                "[INFO] {0}/{1} {2}: {3}".format(
                    completed_count, len(jobs), gene, row["Status"]
                )
            )

    rows.sort(key=lambda row: row["Gene"])
    q_values = benjamini_hochberg([row["P_value"] for row in rows])
    for row, q_value in zip(rows, q_values):
        row["BH_FDR"] = q_value

    primary_summary = (
        Path(args.summary).resolve()
        if args.summary
        else output_dir / "paml_lrt_results.tsv"
    )
    write_summary(rows, primary_summary, args.alpha)
    create_summary_aliases(primary_summary)

    valid_count = sum(
        1
        for row in rows
        if row["Status"].startswith("ok") and math.isfinite(row["P_value"])
    )
    fdr_count = sum(
        1
        for row in rows
        if row["Status"].startswith("ok")
        and math.isfinite(row["BH_FDR"])
        and row["BH_FDR"] <= args.alpha
    )
    failed_count = len(rows) - valid_count

    print("[INFO] Summary: {0}".format(primary_summary))
    print("[INFO] Valid LRT tests: {0}".format(valid_count))
    print("[INFO] BH-FDR <= {0}: {1}".format(args.alpha, fdr_count))
    print("[INFO] Failed/invalid tests: {0}".format(failed_count))

    # Return non-zero only if every test failed. Partial failure remains visible in TSV.
    return 0 if valid_count > 0 else 2


if __name__ == "__main__":
    sys.exit(main())
