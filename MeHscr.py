#!/usr/bin/env python3
"""Profile methylation heterogeneity from bisulfite-sequencing BAM files.

This Python 3 rewrite keeps the original CLI while supporting current pandas
versions and processing BAM files without creating destructive temporary copies.
"""

import argparse
import logging
import math
import os
import random
import sys
from collections import deque
from pathlib import Path

import numpy as np
import pandas as pd
import pysam
from joblib import Parallel, delayed

VERSION = "2.0.0"
CONTEXTS = ("CG", "CHG", "CHH")
RESULT_COLUMNS = ["chrom", "pos", "MeH", "dis", "strand"]
ML_COLUMNS = ["chrom", "pos", "ML", "strand", "depth"]
LOGGER = logging.getLogger("MeH")


def open_log(fname):
    """Configure console and file logging (legacy public helper)."""
    LOGGER.handlers.clear()
    LOGGER.setLevel(logging.INFO)
    formatter = logging.Formatter("[%(asctime)s] %(message)s", "%Y-%m-%d %H:%M:%S")
    for handler in (logging.StreamHandler(sys.stdout), logging.FileHandler(fname, mode="w")):
        handler.setFormatter(formatter)
        LOGGER.addHandler(handler)


def logm(message):
    LOGGER.info(message)


def close_log():
    for handler in list(LOGGER.handlers):
        handler.close()
        LOGGER.removeHandler(handler)


def enough_reads(window, w, complete):
    values = np.asarray(window, dtype=float)
    full = np.isnan(values).sum(axis=1) == 0
    if complete:
        return bool(full.sum() >= 2**w)
    one_missing = np.isnan(values).sum(axis=1) == 1
    return bool(full.sum() >= 2 ** (w - 2) and one_missing.any())


def coverage(methbin, complete, w):
    meth = methbin.drop(columns="Qname", errors="ignore")
    total = len(meth.columns) - w + 1
    if total <= 0:
        return 0.0
    eligible = sum(enough_reads(meth.iloc[:, i:i + w].to_numpy(), w, complete)
                   for i in range(total))
    return eligible / total * 100.0


def impute(window, w, rng=None):
    """Impute rows with one missing value using the original BSImp strategy."""
    values = np.asarray(window, dtype=float).copy()
    chooser = rng or random
    full_indices = np.flatnonzero(np.isnan(values).sum(axis=1) == 0)
    for row_index in np.flatnonzero(np.isnan(values).sum(axis=1) == 1):
        missing_col = int(np.flatnonzero(np.isnan(values[row_index]))[0])
        observed = values[~np.isnan(values[:, missing_col]), missing_col]
        if np.unique(observed).size == 1:
            replacement = observed[0]
        else:
            matches = np.flatnonzero(
                np.sum(values[row_index] == values[full_indices], axis=1) == w - 1
            )
            if matches.size:
                replacement = values[full_indices[chooser.choice(matches.tolist())], missing_col]
            else:
                replacement = chooser.choice(observed.tolist())
        values[row_index, missing_col] = replacement
    return values


def getcomplete(window, w):
    values = np.asarray(window, dtype=float)
    return values[np.isnan(values).sum(axis=1) == 0, :]


def Ham_d(pat1, pat2):
    return int(np.sum(np.asarray(pat1) != np.asarray(pat2)))


def WDK_d(pat1, pat2):
    first, second = np.asarray(pat1), np.asarray(pat2)
    width = first.shape[0]
    distance = 0
    for length_minus_one in range(width):
        weight = width - length_minus_one - 1
        for start in range(width - length_minus_one):
            stop = start + length_minus_one + 1
            distance += weight * int(not np.array_equal(first[start:stop], second[start:stop]))
    return int(distance)


def _all_patterns(width):
    numbers = np.arange(2**width, dtype=np.uint64)
    return ((numbers[:, None] >> np.arange(width, dtype=np.uint64)) & 1).astype(float)


def PattoDis(mat, dist=1):
    values = np.asarray(mat)
    if dist == 1:
        return np.sum(values[:, None, :] != values[None, :, :], axis=2).astype(float)
    size = values.shape[0]
    output = np.zeros((size, size), dtype=float)
    for i in range(size):
        for j in range(i):
            output[i, j] = output[j, i] = WDK_d(values[i], values[j])
    return output


def window_summ(pat, start, dis, chrom):
    values = pat.dropna().to_numpy(dtype=int)
    width = pat.shape[1]
    indices = np.sum(values * (1 << np.arange(width)), axis=1) if len(values) else []
    counts = np.bincount(indices, minlength=2**width) if len(values) else np.zeros(2**width, int)
    row = {"chrom": chrom, "pos": start}
    row.update({"p{:02d}".format(i + 1): int(value) for i, value in enumerate(counts)})
    row["dis"] = dis
    return pd.DataFrame([row])


def _phylogeny_score(counts, width, dist):
    count = np.concatenate(([counts[0]], counts))
    if width == 3:
        leaves = 0.5 if dist == 1 else 1.5
        internal = 0.25 if dist == 1 else 0.75
        tree = np.r_[[0], np.repeat(leaves, 8), np.repeat(internal, 3), leaves, internal]
        nodes = np.zeros(14)
        nodes[1:9] = count[1:9]
        nodes[9] = nodes[1] + nodes[2]
        nodes[10] = nodes[5] + nodes[6]
        nodes[11] = nodes[3] + nodes[4]
        nodes[12] = nodes[9] + nodes[10]
        nodes[13] = nodes[11] + nodes[12]
    elif width == 4:
        nodes = np.zeros(30)
        nodes[1:17] = count[[1, 9, 5, 3, 2, 13, 11, 10, 7, 6, 4, 15, 14, 12, 8, 16]]
        if dist == 1:
            tree = np.r_[[0], np.repeat(0.5, 16), np.repeat(0.25, 6), 0.5,
                         np.repeat(0.25, 6)]
            for target, pair in zip(range(17, 21), ((4, 7), (9, 12), (1, 2), (3, 6))):
                nodes[target] = nodes[pair[0]] + nodes[pair[1]]
            nodes[21] = nodes[17] + nodes[18]
            nodes[22] = nodes[19] + nodes[20]
            nodes[23] = nodes[21] + nodes[22]
            nodes[24] = nodes[5] + nodes[8]
            nodes[25] = nodes[10] + nodes[13]
            nodes[26] = nodes[24] + nodes[25]
            nodes[27] = nodes[23] + nodes[26]
            nodes[28] = nodes[11] + nodes[14]
            nodes[29] = nodes[27] + nodes[28]
        else:
            tree = np.r_[[0], np.repeat(3.0, 16), np.repeat(1.5, 6), 3.2, 0.8,
                         np.repeat(2.0, 3), np.repeat(1.5, 2)]
            pairs = ((1, 2), (5, 8), (3, 6), (10, 13), (4, 7), (11, 14))
            for target, pair in zip(range(17, 23), pairs):
                nodes[target] = nodes[pair[0]] + nodes[pair[1]]
            nodes[23] = nodes[17] + nodes[18]
            nodes[24] = nodes[21] + nodes[22]
            nodes[25] = nodes[19] + nodes[20]
            nodes[26] = nodes[23] + nodes[24]
            nodes[27] = nodes[25] + nodes[26]
            nodes[28] = nodes[9] + nodes[12]
            nodes[29] = nodes[27] + nodes[28]
    else:
        raise ValueError("phylogeny score supports window sizes 3 and 4 only")
    total = float(np.sum(tree * nodes))
    return 0.0 if total == 0 else float(np.sum(tree * (nodes / total) ** 2) ** -1)


def MeHperwindow(pat, start, dis, chrom, D, w, optional, MeH=2, dist=1, strand="f"):
    values = np.asarray(pat, dtype=int)
    indices = np.sum(values * (1 << np.arange(w)), axis=1)
    counts = np.bincount(indices, minlength=2**w).astype(float)
    proportions = counts / float(len(values))
    if MeH == 1:
        score = float(np.sum(proportions**2) ** -1)
    elif MeH == 2:
        interaction = np.multiply.outer(proportions, proportions)
        denominator = float(np.sum(D * interaction))
        score = 0.0 if denominator == 0 else float(
            (np.sum(D * interaction**2) / denominator**2) ** -0.5)
    elif MeH == 3:
        score = _phylogeny_score(counts, w, dist)
    elif MeH == 4:
        nonzero = proportions[proportions > 0]
        score = float(-np.sum(nonzero * np.log2(nonzero)) / w)
    elif MeH == 5:
        score = float(1 - np.sum(proportions**2))
    else:
        raise ValueError("unknown MeH score: {}".format(MeH))
    row = {"chrom": chrom, "pos": start, "MeH": round(score, 5),
           "dis": dis, "strand": strand}
    result = pd.DataFrame([row], columns=RESULT_COLUMNS)
    if not optional:
        return result
    opt = {"chrom": chrom, "pos": start}
    opt.update({"p{:02d}".format(i + 1): int(v) for i, v in enumerate(counts)})
    opt.update({"MeH": round(score, 5), "dis": dis, "strand": strand})
    return result, pd.DataFrame([opt])


def _context_matches(fasta, chrom, pos, context, strand):
    length = fasta.get_reference_length(chrom)
    if strand == "f":
        if pos + 3 > length:
            return False
        sequence = fasta.fetch(chrom, pos, pos + 3).upper()
    else:
        if pos < 2:
            return False
        sequence = fasta.fetch(chrom, pos - 2, pos + 1).upper()
    # BAM headers and a supplied FASTA can disagree about contig length.  A
    # truncated motif cannot be classified safely.
    if len(sequence) != 3:
        return False
    if context == "CG":
        return sequence[:2] == "CG" if strand == "f" else sequence[1:] == "CG"
    if context == "CHG":
        if strand == "f":
            return sequence[0] == "C" and sequence[1] != "G" and sequence[2] == "G"
        return sequence[0] == "C" and sequence[1] != "C" and sequence[2] == "G"
    if strand == "f":
        return sequence[0] == "C" and sequence[1] != "G" and sequence[2] != "G"
    return sequence[2] == "G" and sequence[1] != "C" and sequence[0] != "C"


def _base_calls(column, reverse):
    methylated, unmethylated = ("G", "A") if reverse else ("C", "T")
    calls = {}
    for read in column.pileups:
        if read.is_del or read.is_refskip or read.alignment.is_reverse != reverse:
            continue
        query_position = read.query_position
        sequence = read.alignment.query_sequence
        if query_position is None or sequence is None:
            continue
        base = sequence[query_position].upper()
        value = 1.0 if base == methylated else 0.0 if base == unmethylated else np.nan
        calls.setdefault(read.alignment.query_name, value)
    return calls


def _evaluate_window(entries, chrom, strand, distance_matrix, width, optional,
                     score_type, distance_type, do_impute, rng):
    positions = [entry[0] for entry in entries]
    matrix = pd.DataFrame([entry[1] for entry in entries], index=positions).T.to_numpy(dtype=float)
    if do_impute and enough_reads(matrix, width, complete=False):
        matrix = impute(matrix, width, rng)
    if not enough_reads(matrix, width, complete=True):
        return None
    return MeHperwindow(
        pd.DataFrame(getcomplete(matrix, width)), positions[0], positions[-1] - positions[0],
        chrom, distance_matrix, width, optional, score_type, distance_type, strand)


def _write_csv(frame, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    frame.to_csv(temporary, index=False)
    os.replace(str(temporary), str(path))


def scan_bam(bam_path, fasta_path, output_dir, context, width=4, optional=False,
             methylation_level=False, distance_type=1, score_type=2,
             do_impute=False, quiet=False, seed=0):
    """Scan one BAM and write per-window results."""
    bam_path, fasta_path, output_dir = map(Path, (bam_path, fasta_path, output_dir))
    sample = bam_path.stem
    results, optional_results, levels = [], [], []
    buffers = {"f": deque(maxlen=width), "r": deque(maxlen=width)}
    distance_matrix = PattoDis(pd.DataFrame(_all_patterns(width)), distance_type)
    rng = random.Random("{}:{}:{}".format(seed, sample, context))
    pileup_coverage = context_coverage = 0
    current_chrom = None

    with pysam.AlignmentFile(str(bam_path), "rb") as alignment, \
            pysam.FastaFile(str(fasta_path)) as fasta:
        for column in alignment.pileup(stepper="all", min_base_quality=0):
            pileup_coverage += 1
            chrom = column.reference_name
            if chrom != current_chrom:
                buffers["f"].clear()
                buffers["r"].clear()
                current_chrom = chrom
            if not quiet and column.pos % 2_000_000 == 1:
                LOGGER.info("%s sample=%s chrom=%s pos=%s results=%s", context,
                            sample, chrom, column.pos + 1, len(results))
            for strand, reverse in (("f", False), ("r", True)):
                if not _context_matches(fasta, chrom, column.pos, context, strand):
                    continue
                context_coverage += 1
                calls = _base_calls(column, reverse)
                if methylation_level:
                    observed = np.asarray([v for v in calls.values() if not math.isnan(v)])
                    if observed.size > 3:
                        levels.append({"chrom": chrom, "pos": column.pos + 1,
                                       "ML": float(observed.mean()), "strand": strand,
                                       "depth": int(observed.size)})
                if not calls:
                    continue
                buffers[strand].append((column.pos + 1, calls))
                if len(buffers[strand]) < width:
                    continue
                evaluated = _evaluate_window(
                    list(buffers[strand]), chrom, strand, distance_matrix, width,
                    optional, score_type, distance_type, do_impute, rng)
                if evaluated is None:
                    continue
                if optional:
                    result, composition = evaluated
                    optional_results.append(composition.iloc[0].to_dict())
                else:
                    result = evaluated
                results.append(result.iloc[0].to_dict())

    _write_csv(pd.DataFrame(results, columns=RESULT_COLUMNS),
               output_dir / "{}_{}_0.csv".format(context, sample))
    if optional:
        columns = (["chrom", "pos"] + ["p{:02d}".format(i + 1) for i in range(2**width)]
                   + ["MeH", "dis", "strand"])
        _write_csv(pd.DataFrame(optional_results, columns=columns),
                   output_dir / "{}_opt_{}.csv".format(context, sample))
    if methylation_level:
        _write_csv(pd.DataFrame(levels, columns=ML_COLUMNS),
                   output_dir / "{}_ML_{}_0.csv".format(context, sample))
    return sample, pileup_coverage, context_coverage, context


def _legacy_scan(context, bamfile, w, fa, optional, melv, silence=False,
                 dist=1, MeH=2, imp=True):
    data_dir = Path("MeHdata")
    return scan_bam(data_dir / (bamfile + ".bam"), data_dir / (fa + ".fa"), data_dir,
                    context, w, optional, melv, dist, MeH, imp, silence)


def CGgenome_scr(bamfile, w, fa, optional, melv, silence=False, dist=1, MeH=2, imp=True):
    return _legacy_scan("CG", bamfile, w, fa, optional, melv, silence, dist, MeH, imp)


def CHGgenome_scr(bamfile, w, fa, optional, melv, silence=False, dist=1, MeH=2, imp=True):
    return _legacy_scan("CHG", bamfile, w, fa, optional, melv, silence, dist, MeH, imp)


def CHHgenome_scr(bamfile, w, fa, optional, melv, silence=False, dist=1, MeH=2, imp=True):
    return _legacy_scan("CHH", bamfile, w, fa, optional, melv, silence, dist, MeH, imp)


def _bin_results(frame, value, minimum_sites=1):
    if frame.empty:
        return pd.DataFrame(columns=["chrom", "bin", "strand", value])
    data = frame.copy()
    data["bin"] = ((data["pos"].astype(int) - 1) // 400) * 400 + 200
    if minimum_sites > 1:
        counts = data.groupby(["chrom", "bin", "strand"]).size().reset_index(name="site_count")
        data = data.merge(counts, on=["chrom", "bin", "strand"])
        data = data[data["site_count"] >= minimum_sites]
    return data.groupby(["chrom", "bin", "strand"], as_index=False)[value].mean()


def _merge_context(context, samples, output_dir, include_ml):
    merged = None
    for sample in samples:
        path = output_dir / "{}_{}_0.csv".format(context, sample)
        binned = _bin_results(pd.read_csv(path), "MeH").rename(columns={"MeH": sample})
        merged = binned if merged is None else merged.merge(
            binned, on=["chrom", "bin", "strand"], how="inner")
    _write_csv(merged if merged is not None else pd.DataFrame(),
               output_dir / "{}_Results.csv".format(context))
    if include_ml:
        merged_ml = None
        for sample in samples:
            path = output_dir / "{}_ML_{}_0.csv".format(context, sample)
            binned = _bin_results(pd.read_csv(path), "ML", 5).rename(columns={"ML": sample})
            merged_ml = binned if merged_ml is None else merged_ml.merge(
                binned, on=["chrom", "bin", "strand"], how="inner")
        _write_csv(merged_ml if merged_ml is not None else pd.DataFrame(),
                   output_dir / "{}_ML_Results.csv".format(context))


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-w", "--windowsize", type=int, default=4, choices=range(2, 7),
                        metavar="{2..6}", help="number of cytosines per window")
    parser.add_argument("-c", "--cores", type=int, default=4,
                        help="number of BAM/context jobs run in parallel")
    parser.add_argument("-m", "--MeH", dest="meh", type=int, default=2,
                        choices=range(1, 6), help="1 abundance, 2 pairwise, 3 phylogeny, 4 entropy, 5 epipolymorphism")
    parser.add_argument("-d", "--dist", type=int, default=1, choices=(1, 2),
                        help="1 Hamming, 2 weighted-degree kernel")
    for context in CONTEXTS:
        parser.add_argument("--" + context, action="store_true",
                            help="include {} context".format(context))
    parser.add_argument("--opt", action="store_true", help="output pattern compositions")
    parser.add_argument("--mlv", action="store_true", help="output methylation levels")
    parser.add_argument("--imp", action="store_true", help="enable BSImp imputation")
    parser.add_argument("--data-dir", type=Path, default=Path("MeHdata"),
                        help="directory containing BAM and FASTA files")
    parser.add_argument("--output-dir", type=Path,
                        help="output directory (defaults to data directory)")
    parser.add_argument("--reference", type=Path,
                        help="reference FASTA when data-dir has multiple FASTA files")
    parser.add_argument("--seed", type=int, default=0, help="imputation random seed")
    parser.add_argument("--quiet", action="store_true", help="suppress progress messages")
    parser.add_argument("--version", action="version", version="%(prog)s " + VERSION)
    return parser


def _validate_inputs(parser, args):
    contexts = [context for context in CONTEXTS if getattr(args, context)]
    if not contexts:
        parser.error("select at least one context: --CG, --CHG, or --CHH")
    if args.cores < 1:
        parser.error("--cores must be at least 1")
    if args.meh == 3 and args.windowsize not in (3, 4):
        parser.error("phylogeny score (-m 3) supports window sizes 3 and 4 only")
    data_dir = args.data_dir.resolve()
    if not data_dir.is_dir():
        parser.error("data directory does not exist: {}".format(data_dir))
    bam_paths = sorted(data_dir.glob("*.bam"))
    if not bam_paths:
        parser.error("no BAM files found in {}".format(data_dir))
    candidates = ([args.reference.resolve()] if args.reference else
                  sorted(data_dir.glob("*.fa")) + sorted(data_dir.glob("*.fasta")))
    if len(candidates) != 1:
        parser.error("specify --reference unless data-dir contains exactly one .fa/.fasta")
    if not candidates[0].is_file():
        parser.error("reference FASTA does not exist: {}".format(candidates[0]))
    return contexts, bam_paths, candidates[0], (args.output_dir or data_dir).resolve()


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    contexts, bam_paths, fasta_path, output_dir = _validate_inputs(parser, args)
    output_dir.mkdir(parents=True, exist_ok=True)
    open_log(str(output_dir / "MeHscreening.log"))
    logm("MeH {}: {} BAM file(s), contexts {}".format(
        VERSION, len(bam_paths), ", ".join(contexts)))
    try:
        jobs = [delayed(scan_bam)(
            bam, fasta_path, output_dir, context, args.windowsize, args.opt,
            args.mlv, args.dist, args.meh, args.imp, args.quiet, args.seed)
            for context in contexts for bam in bam_paths]
        stats = Parallel(n_jobs=min(args.cores, len(jobs)))(jobs)
        samples = [path.stem for path in bam_paths]
        for context in contexts:
            _merge_context(context, samples, output_dir, args.mlv)
        for sample, total, context_total, context in stats:
            logm("Sample {} has coverage {} for context {} out of data coverage {}.".format(
                sample, context_total, context, total))
        logm("All done.")
        return 0
    except Exception:
        LOGGER.exception("MeH failed")
        return 1
    finally:
        close_log()


if __name__ == "__main__":
    raise SystemExit(main())
