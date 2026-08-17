#!/usr/bin/env python3
"""Fit gapped-local-alignment extreme-value constants from null histograms.

The fitted model is the BLAST form

  P(max_score < S) = exp(-K * A * exp(-lambda * S))

where A uses the standard finite-length correction with alpha and beta. This
script is a release/calibration tool; neither web nor CLI performs simulation.
"""

from __future__ import annotations

import csv
import json
import math
import sys
import numpy as np
from scipy.optimize import brentq, least_squares


def length_adjustment(m: int, n: int, sequences: int, lam: float, k: float, alpha: float, beta: float) -> float:
    upper = min(float(m - 1), max(0.0, (n - sequences) / max(1, sequences)))
    ell = 0.0
    for _ in range(64):
        area = max(1.0, (m - ell) * (n - sequences * ell))
        candidate = (alpha / lam) * (math.log(k) + math.log(area)) + beta
        candidate = min(upper, max(0.0, candidate))
        if abs(candidate - ell) < 1e-10:
            return candidate
        ell = 0.5 * (ell + candidate)
    return ell


def main(path: str) -> None:
    groups: dict[tuple[str, int, int], dict[str, object]] = {}
    scoring: dict[str, tuple[int, int, int, int]] = {}
    with open(path, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            profile = row["profile"]
            m, n = int(row["query_length"]), int(row["reference_length"])
            scoring[profile] = tuple(int(row[key]) for key in ("match", "mismatch", "gap_open", "gap_extend"))
            group = groups.setdefault((profile, m, n), {"replicates": int(row["replicates"]), "histogram": {}})
            group["histogram"][int(row["score"])] = int(row["count"])

    output: dict[str, object] = {}
    for profile in sorted(scoring):
        match, mismatch, _, _ = scoring[profile]
        ungapped_lambda = float(brentq(
            lambda value: 0.25 * math.exp(match * value)
            + 0.75 * math.exp(mismatch * value) - 1.0,
            1e-8,
            10.0,
        ))
        observations: list[tuple[int, int, int, float, float]] = []
        # Use exact P(max < S), matching the threshold convention in E(S).
        for (candidate, m, n), group in groups.items():
            if candidate != profile:
                continue
            # Transposed matrices have the same null distribution. Fit one
            # orientation and retain the independently simulated transpose as
            # a held-out validation set.
            if m > n:
                continue
            total = int(group["replicates"])
            histogram = group["histogram"]
            cumulative = 0
            for score in range(min(histogram), max(histogram) + 2):
                if score - 1 in histogram:
                    cumulative += histogram[score - 1]
                # Jeffreys smoothing prevents endpoints from becoming infinite.
                probability = (cumulative + 0.5) / (total + 1.0)
                if 0.015 <= probability <= 0.985:
                    # Delta-method precision for log(-log(F)); cap prevents a
                    # single well-sampled center point dominating every length.
                    precision = math.sqrt(total) * (-math.log(probability)) * math.sqrt(probability / (1.0 - probability))
                    observations.append((m, n, score, math.log(-math.log(probability)), min(precision, 80.0)))

        def residual(parameters: np.ndarray) -> np.ndarray:
            lam = math.exp(parameters[0])
            k = math.exp(parameters[1])
            alpha = math.exp(parameters[2])
            beta = parameters[3]
            values = []
            for m, n, score, observed, weight in observations:
                ell = length_adjustment(m, n, 1, lam, k, alpha, beta)
                area = max(1.0, (m - ell) * (n - ell))
                predicted = math.log(k) + math.log(area) - lam * score
                values.append((predicted - observed) * math.sqrt(weight))
            return np.asarray(values)

        fitted = least_squares(
            residual,
            np.asarray([math.log(0.55), math.log(0.1), math.log(1.0), 0.0]),
            # Allowing a gapped lambda above its ungapped limit can fit small-
            # matrix discreteness but is not asymptotically valid.
            bounds=(np.asarray([math.log(0.05), math.log(1e-5), math.log(0.01), -20.0]),
                    np.asarray([math.log(ungapped_lambda * (1.0 - 1e-8)), math.log(10.0), math.log(20.0), 20.0])),
            max_nfev=20000,
        )
        lam, k, alpha = map(math.exp, fitted.x[:3])
        beta = float(fitted.x[3])

        validation = []
        for (candidate, m, n), group in sorted(groups.items()):
            if candidate != profile:
                continue
            total = int(group["replicates"])
            histogram = group["histogram"]
            scores = np.asarray([score for score, count in histogram.items() for _ in range(count)])
            ell = length_adjustment(m, n, 1, lam, k, alpha, beta)
            area = max(1.0, (m - ell) * (n - ell))
            quantiles = []
            for q in (0.5, 0.9, 0.95, 0.99):
                observed_score = int(np.quantile(scores, q, method="higher"))
                predicted_cdf = math.exp(-k * area * math.exp(-lam * (observed_score + 1)))
                quantiles.append({"quantile": q, "score": observed_score, "predicted_cdf_below_next_score": predicted_cdf})
            validation.append({"query_length": m, "reference_length": n, "replicates": total,
                               "split": "fit" if m <= n else "held_out", "quantiles": quantiles})

        output[profile] = {
            "scoring": scoring[profile],
            "lambda": lam,
            "K": k,
            "alpha": alpha,
            "beta": beta,
            "ungapped_lambda": ungapped_lambda,
            "weighted_rmse": float(np.sqrt(np.mean(np.square(residual(fitted.x))))),
            "validation": validation,
        }

    print(json.dumps(output, indent=2, sort_keys=True))


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit(f"usage: {sys.argv[0]} NULL_SCORES.tsv")
    main(sys.argv[1])
