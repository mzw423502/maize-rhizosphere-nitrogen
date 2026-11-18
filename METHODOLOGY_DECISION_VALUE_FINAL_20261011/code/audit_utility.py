"""Independent audit for the fixed-split decision-utility analysis.

This module deliberately does not import utility.py or any prior producer or
certificate module.  It reads immutable saved annual records and reconstructs
the calibration domains, interval-containment decisions, dataset summaries,
and MCSEs.  It never fits a model or resamples a count dataset.
"""
import gzip
import hashlib
import json
import math
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
PRIOR = ROOT.parent / "METHODOLOGY_E2E_ITERATION_20261011"
METHODS = ("LinDA_unit", "CLR_unit")
CAPS = (0.6, 0.7, 0.8, 0.9)
DOMAINS = ("zero", "calibrated_point", "calibrated_interval", "full")
HALFWIDTHS = (0.0, 0.1, 0.25, 0.5, 1.0, 2.0)
ATOL = 1e-10
RTOL = 1e-10


def digest(path):
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def load_json(path):
    with path.open("r", encoding="utf8") as handle:
        return json.load(handle)


def annual_record(cell, replicate):
    path = PRIOR / "analysis" / "formal" / f"c{cell:02d}_r{replicate:04d}.json.gz"
    with gzip.open(path, "rt", encoding="utf8") as handle:
        row = json.load(handle)
    if row["cell"] != cell or row["replicate"] != replicate:
        raise ValueError("saved annual-record identity mismatch")
    return path, row


def ar(obj, year, field):
    return np.asarray(obj["annual"][year][field], dtype=float)


def x_of_rho(rho, s1, s2):
    # This squared-difference form avoids a subtraction cancellation at rho=-1.
    return ((s1 - s2) ** 2 + 2.0 * (1.0 + rho) * s1 * s2) / 4.0


def permitted_x_band(z, d, cap):
    """Independent inverse-window calculation for W(x,z,d)<=cap."""
    cap = np.asarray(cap, dtype=float)
    if not np.all((0.5 <= cap) & (cap < 1.0)):
        raise ValueError("cap outside declared range")
    c = 2.0 * cap - 1.0
    d2 = d * d
    # The two branches meet at z; their union is a single interval.
    lower = np.maximum(0.0, np.minimum(z * (1.0 - cap) / cap, z - c * d2))
    upper = np.maximum(z * cap / (1.0 - cap), z + c * d2)
    return lower, upper


def contained_in_window(rlo, rhi, s1, s2, z, d, cap):
    lo = x_of_rho(rlo, s1, s2)
    hi = x_of_rho(rhi, s1, s2)
    low, high = permitted_x_band(z, d, cap)
    # A zero variance endpoint is a balance-limit calculation, not a regular fit.
    return (lo >= low) & (hi <= high)


def required_cap_by_containment(rlo, rhi, s1, s2, z, d):
    """Find the minimal cap independently by inverse-window bisection."""
    answer = np.empty_like(z, dtype=float)
    # A cap of one may be the exact one-sided supremum at a zero endpoint.
    possible_below_one = contained_in_window(rlo, rhi, s1, s2, z, d, 1.0 - 1e-13)
    answer[~possible_below_one] = 1.0
    active = np.flatnonzero(possible_below_one)
    if active.size:
        left = np.full(active.size, 0.5)
        right = np.full(active.size, 1.0 - 1e-13)
        for _ in range(80):
            middle = (left + right) / 2.0
            ok = contained_in_window(
                np.asarray(rlo)[active], np.asarray(rhi)[active], s1[active], s2[active],
                z[active], d[active], middle,
            )
            right[ok] = middle[ok]
            left[~ok] = middle[~ok]
        answer[active] = right
    return answer


def equalish(a, b):
    if a is None or b is None:
        return a is None and b is None
    if isinstance(a, bool) or isinstance(b, bool):
        return bool(a) is bool(b)
    try:
        return bool(math.isclose(float(a), float(b), rel_tol=RTOL, abs_tol=ATOL))
    except (TypeError, ValueError):
        return a == b


def fail(diffs, location, expected, observed):
    diffs.append({"location": location, "expected": expected, "observed": observed})


def compare_scalar(diffs, location, expected, observed):
    if not equalish(expected, observed):
        fail(diffs, location, expected, observed)


def compare_array(diffs, location, expected, observed):
    if len(expected) != len(observed):
        fail(diffs, location + ".length", len(expected), len(observed))
        return
    for index, (want, got) in enumerate(zip(expected, observed)):
        compare_scalar(diffs, f"{location}[{index}]", want, got)


def truth_masks(truth_a, truth_b):
    same = truth_a * truth_b > 1e-20
    opposite = truth_a * truth_b < -1e-20
    null = (np.abs(truth_a) < 1e-14) & (np.abs(truth_b) < 1e-14)
    one_zero = ~(same | opposite | null)
    return same, opposite, null, one_zero


def ratio(values, denominator):
    return float(np.mean(values[denominator])) if np.any(denominator) else None


def dataset_metrics(selected, truth_a, truth_b, estimate_a, at_zero, full, snr, iv_cap, available):
    same, opposite, null, one_zero = truth_masks(truth_a, truth_b)
    wrong = same & (estimate_a * truth_a < 0.0)
    removed = at_zero & ~full & available
    out = {
        "membership": ratio(selected, available),
        "true_same_retention": ratio(selected, same & available),
        "true_same_wrong_sign_selected_rate": ratio(selected & wrong, same & available),
        "true_opposite_membership": ratio(selected, opposite & available),
        "null_membership": ratio(selected, null & available),
        "one_zero_membership": ratio(selected, one_zero & available),
        "opposite_fraction_among_selected": ratio(opposite, selected),
        "null_fraction_among_selected": ratio(null, selected),
        "paired_difference_from_zero": ratio(selected.astype(int) - at_zero.astype(int), available),
        "full_loss_fraction_of_this_strategy": ratio(~full, selected),
        "true_same_paired_difference_from_zero": ratio(
            selected.astype(int) - at_zero.astype(int), same & available
        ),
        "true_opposite_paired_difference_from_zero": ratio(
            selected.astype(int) - at_zero.astype(int), opposite & available
        ),
        "removed_from_zero_true_same_fraction": ratio(same, removed),
        "removed_from_zero_true_opposite_fraction": ratio(opposite, removed),
        "removed_from_zero_null_fraction": ratio(null, removed),
        "selected_true_same_correct_sign_retention": ratio(selected & ~wrong, same & available),
    }
    for label, mask in (
        ("low", snr < 1.0),
        ("middle", (snr >= 1.0) & (snr < 2.0)),
        ("high", snr >= 2.0),
    ):
        out["true_same_retention_snr_" + label] = ratio(selected, same & mask & available)
    for label, mask in (
        ("balanced", iv_cap <= 0.7),
        ("moderate", (iv_cap > 0.7) & (iv_cap <= 0.9)),
        ("dominated", iv_cap > 0.9),
    ):
        out["true_same_retention_IV_" + label] = ratio(selected, same & mask & available)
    magnitude = np.round(np.minimum(np.abs(truth_a), np.abs(truth_b)), 12)
    for value in sorted(set(magnitude[same])):
        out["true_same_retention_magnitude_" + str(float(value))] = ratio(
            selected, same & (magnitude == value) & available
        )
    den = {
        "all_taxa": len(estimate_a),
        "available": int(available.sum()),
        "unavailable": int((~available).sum()),
        "true_same": int((same & available).sum()),
        "true_opposite": int((opposite & available).sum()),
        "null": int((null & available).sum()),
        "one_zero": int((one_zero & available).sum()),
        "selected": int(selected.sum()),
        "rho0_selected": int((at_zero & available).sum()),
        "full_selected": int((full & available).sum()),
        "rho0_removed_by_full": int(removed.sum()),
    }
    out["full_minus_strategy_membership"] = ratio(full.astype(int) - selected.astype(int), available)
    out["full_minus_strategy_true_same_retention"] = ratio(
        full.astype(int) - selected.astype(int), same & available
    )
    out["full_minus_strategy_true_opposite_membership"] = ratio(
        full.astype(int) - selected.astype(int), opposite & available
    )
    for label, mask in (
        ("low", snr < 1.0),
        ("middle", (snr >= 1.0) & (snr < 2.0)),
        ("high", snr >= 2.0),
    ):
        den["true_same_snr_" + label] = int((same & mask & available).sum())
    for label, mask in (
        ("balanced", iv_cap <= 0.7),
        ("moderate", (iv_cap > 0.7) & (iv_cap <= 0.9)),
        ("dominated", iv_cap > 0.9),
    ):
        den["true_same_IV_" + label] = int((same & mask & available).sum())
    for value in sorted(set(magnitude[same])):
        den["true_same_magnitude_" + str(float(value))] = int(
            (same & (magnitude == value) & available).sum()
        )
    return out, den


class DatasetMoments:
    def __init__(self):
        self.n = 0
        self.mean = 0.0
        self.m2 = 0.0

    def add(self, value):
        if value is None:
            return
        self.n += 1
        delta = value - self.mean
        self.mean += delta / self.n
        self.m2 += delta * (value - self.mean)

    def payload(self):
        return {
            "mean": self.mean if self.n else None,
            "MCSE": math.sqrt(max(0.0, self.m2) / (self.n - 1) / self.n) if self.n > 1 else None,
            "n_datasets": self.n,
        }


def edge_math_checks():
    """Direct likelihood/weight checks on bounded edge cases, independent of files."""
    checks = 0
    for x, z, d in ((0.2, 0.8, 0.0), (0.2, 0.8, 1.5), (2.0, 0.3, -0.7)):
        d2 = d * d
        tau = max(0.0, (d2 - x - z) / 2.0)
        va, vb = x + tau, z + tau
        normalized = max(1.0 / va, 1.0 / vb) / (1.0 / va + 1.0 / vb)
        explicit = 0.5 + abs(x - z) / (2.0 * max(x + z, d2))
        if not math.isclose(normalized, explicit, abs_tol=1e-12, rel_tol=1e-12):
            raise AssertionError("direct normalized-weight identity failed")
        # The actual two-point weighted-residual identity is checked explicitly below.
        mean = ((1.0 / va) * d + (1.0 / vb) * 0.0) / (1.0 / va + 1.0 / vb)
        Q = (1.0 / va) * (d - mean) ** 2 + (1.0 / vb) * mean ** 2
        expected_q = min(d2 / (x + z), 1.0)
        if not math.isclose(Q, expected_q, abs_tol=1e-12, rel_tol=1e-12):
            raise AssertionError("HK residual identity failed")
        # Direct restricted likelihood checks the analytic constrained minimizer.
        def nll(t):
            total = x + z + 2.0 * t
            return math.log(total) + d2 / total
        if tau == 0.0:
            if nll(1e-7) < nll(0.0) - 1e-12:
                raise AssertionError("boundary restricted-likelihood minimum failed")
        else:
            eps_t = min(1e-7, tau / 2.0)
            if nll(tau - eps_t) < nll(tau) - 1e-12 or nll(tau + eps_t) < nll(tau) - 1e-12:
                raise AssertionError("interior restricted-likelihood minimum failed")
        checks += 3
    # x=0 has no regular fit.  Its one-sided balance limit is nevertheless W=1
    # when d=0,z>0; verify approach from positive x without dividing by zero.
    z, d = 1.0, 0.0
    limiting = 0.5 + abs(0.0 - z) / (2.0 * max(z, d * d))
    eps = 1e-12
    tau = max(0.0, (d * d - eps - z) / 2.0)
    near = max(1.0 / (eps + tau), 1.0 / (z + tau)) / (1.0 / (eps + tau) + 1.0 / (z + tau))
    if not (limiting == 1.0 and abs(near - limiting) <= 2e-12):
        raise AssertionError("one-sided zero-variance balance limit failed")
    checks += 1
    # At q=.5, only exact equality of the effective variances is admissible.
    lo, hi = permitted_x_band(np.asarray([1.0]), np.asarray([0.0]), 0.5)
    if not (lo[0] == 1.0 and hi[0] == 1.0):
        raise AssertionError("q=.5 window boundary failed")
    checks += 1
    # Positive rectangular boxes: test both an interval containing zero and one not containing zero.
    for d0, d1 in ((-0.4, 1.0), (0.4, 1.0)):
        x0, x1, z0, z1 = 0.2, 1.7, 0.3, 1.2
        dstar = 0.0 if d0 <= 0.0 <= d1 else min(abs(d0), abs(d1))
        corners = []
        for x in (x0, x1):
            for z in (z0, z1):
                tau = max(0.0, (dstar * dstar - x - z) / 2.0)
                corners.append(max(1.0 / (x + tau), 1.0 / (z + tau)) / (1.0 / (x + tau) + 1.0 / (z + tau)))
        corner_max = max(corners)
        for x in np.linspace(x0, x1, 17):
            for z in np.linspace(z0, z1, 17):
                tau = max(0.0, (dstar * dstar - x - z) / 2.0)
                test = max(1.0 / (x + tau), 1.0 / (z + tau)) / (1.0 / (x + tau) + 1.0 / (z + tau))
                if test > corner_max + 1e-12:
                    raise AssertionError("box-corner bound failed")
        checks += 1
    return checks


def lookup_calibration(rows):
    return {(row["cell"], row["method"]): row for row in rows}


def audit(scale):
    if scale not in ("pilot", "formal"):
        raise ValueError("scale must be pilot or formal")
    out = ROOT / "results" / scale
    result = load_json(out / "RESULT.json")
    if result.get("status") != "AWAITING_INDEPENDENT_AUDIT":
        raise ValueError("producer result is not awaiting an independent audit")
    if result.get("new_fits") or result.get("new_simulated_datasets"):
        raise ValueError("utility stage must reuse immutable saved components")
    expected_cells = [0, 1, 2, 25] if scale == "pilot" else list(range(26))
    evaluation = range(100, 102) if scale == "pilot" else range(100, 200)
    manifest = load_json(PRIOR / "data" / "formal" / "MANIFEST.json")
    seal_path = PRIOR / "provenance" / "formal_POST_AUDIT_SEAL.json"
    seal = load_json(seal_path)["files"]
    diffs = []
    accessed = {}

    # The producer may neither omit a saved component nor silently add one.
    # This is deliberately checked before any row-level reconstruction.
    expected_source_names = {
        f"c{cell:02d}_r{replicate:04d}.json.gz"
        for cell in expected_cells
        for replicate in tuple(range(100)) + tuple(evaluation)
    }
    if set(result.get("source_hashes", {})) != expected_source_names:
        fail(diffs, "RESULT.source_hashes.keyset", sorted(expected_source_names), sorted(result.get("source_hashes", {})))
    if Path(result.get("source_root", "")).resolve() != PRIOR.resolve():
        fail(diffs, "RESULT.source_root", str(PRIOR.resolve()), result.get("source_root"))
    compare_scalar(diffs, "RESULT.prior_seal_sha256", digest(seal_path), result.get("prior_seal_sha256"))
    compare_scalar(diffs, "RESULT.calibration_replicates", [0, 99], result.get("calibration_replicates"))
    compare_scalar(diffs, "RESULT.evaluation_replicates", [100, max(evaluation)], result.get("evaluation_replicates"))
    compare_scalar(diffs, "RESULT.scenario_count", len(expected_cells), result.get("scenario_count"))
    expected_summary_count = len(expected_cells) * len(METHODS) * len(CAPS) * (len(DOMAINS) + len(HALFWIDTHS))
    compare_scalar(diffs, "RESULT.summary_rows", expected_summary_count, result.get("summary_rows"))
    compare_scalar(diffs, "RESULT.original_native_exit_code", None, result.get("original_native_exit_code"))
    compare_scalar(diffs, "RESULT.truth_ledger_separate", True, result.get("truth_ledger_separate"))
    compare_scalar(diffs, "RESULT.FDR_claim", False, result.get("FDR_claim"))

    # Verify that every producer-read immutable source is still byte-identical to the post-audit seal.
    for name, expected_hash in result["source_hashes"].items():
        source = PRIOR / "analysis" / "formal" / name
        actual_hash = digest(source)
        accessed[name] = actual_hash
        compare_scalar(diffs, "source_hashes." + name, expected_hash, actual_hash)
        compare_scalar(diffs, "post_audit_seal.analysis/formal/" + name, seal.get("analysis/formal/" + name), actual_hash)

    raw_calibration = load_json(out / "CALIBRATION.json")
    stored_cal = lookup_calibration(raw_calibration)
    expected_calibration_keys = {(cell, method) for cell in expected_cells for method in METHODS}
    if len(raw_calibration) != len(stored_cal) or set(stored_cal) != expected_calibration_keys:
        fail(
            diffs,
            "CALIBRATION.keyset_and_uniqueness",
            {"keys": sorted(expected_calibration_keys), "rows": len(expected_calibration_keys)},
            {"keys": sorted(stored_cal), "rows": len(raw_calibration)},
        )
    raw_summary = load_json(out / "SUMMARY.json")
    stored_summary = {
        (row["cell"], row["method"], float(row["q"]), row["domain"]): row
        for row in raw_summary
    }
    expected_summary_keys = {
        (cell, method, cap, domain)
        for cell in expected_cells
        for method in METHODS
        for cap in CAPS
        for domain in DOMAINS + tuple("width_" + str(width) for width in HALFWIDTHS)
    }
    if len(raw_summary) != len(stored_summary) or set(stored_summary) != expected_summary_keys:
        fail(
            diffs,
            "SUMMARY.keyset_and_uniqueness",
            {"keys": sorted(expected_summary_keys), "rows": len(expected_summary_keys)},
            {"keys": sorted(stored_summary), "rows": len(raw_summary)},
        )
    stored_cases = load_json(out / "ILLUSTRATIVE_CASES.json")
    grouped = defaultdict(lambda: defaultdict(DatasetMoments))
    den_hist = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
    expected_cases = {}
    decisions_seen = 0
    dataset_seen = 0

    # Rebuild calibration correlations from saved annual estimates/truths.
    calibration = {}
    for cell in expected_cells:
        arrays = {method: {year: [] for year in ("A1", "A2")} for method in METHODS}
        for replicate in range(100):
            path, row = annual_record(cell, replicate)
            h = digest(path)
            accessed[path.name] = h
            compare_scalar(diffs, "post_audit_seal.analysis/formal/" + path.name, seal.get("analysis/formal/" + path.name), h)
            for method in METHODS:
                obj = row["methods"][method]
                for year in ("A1", "A2"):
                    arrays[method][year].append(ar(obj, year, "effect") - ar(obj, year, "truth"))
        for method in METHODS:
            e1 = np.asarray(arrays[method]["A1"], dtype=float)
            e2 = np.asarray(arrays[method]["A2"], dtype=float)
            c1 = e1 - e1.mean(axis=0)
            c2 = e2 - e2.mean(axis=0)
            denom = np.sqrt((c1 * c1).sum(axis=0) * (c2 * c2).sum(axis=0))
            rho = np.divide((c1 * c2).sum(axis=0), denom, out=np.full(100, np.nan), where=denom > 0)
            defined = np.isfinite(rho) & (np.abs(rho) < 1.0)
            lower = np.full(100, np.nan)
            upper = np.full(100, np.nan)
            lower[defined] = np.maximum(-1.0, np.tanh(np.arctanh(rho[defined]) - 1.96 / math.sqrt(97.0)))
            upper[defined] = np.minimum(1.0, np.tanh(np.arctanh(rho[defined]) + 1.96 / math.sqrt(97.0)))
            defined &= np.isfinite(lower) & np.isfinite(upper) & (lower <= upper)
            rho[~defined] = np.nan
            lower[~defined] = np.nan
            upper[~defined] = np.nan
            calibration[(cell, method)] = (rho, lower, upper, defined)
            got = stored_cal.get((cell, method))
            if got is None:
                fail(diffs, f"CALIBRATION[{cell},{method}]", "row", None)
            else:
                want_rho = [None if not np.isfinite(v) else float(v) for v in rho]
                want_lo = [None if not np.isfinite(v) else float(v) for v in lower]
                want_hi = [None if not np.isfinite(v) else float(v) for v in upper]
                compare_array(diffs, f"CALIBRATION[{cell},{method}].rho", want_rho, got["rho"])
                compare_array(diffs, f"CALIBRATION[{cell},{method}].lower", want_lo, got["lower"])
                compare_array(diffs, f"CALIBRATION[{cell},{method}].upper", want_hi, got["upper"])
                compare_array(diffs, f"CALIBRATION[{cell},{method}].availability", defined.tolist(), got["availability"])
                compare_scalar(diffs, f"CALIBRATION[{cell},{method}].n_calibration", 100, got["n_calibration"])
        del arrays

    # Stream producer taxon rows and dataset rows in their deterministic order.
    taxon_path = out / "TAXON_DECISIONS.jsonl.gz"
    dataset_path = out / "DATASET_RECORDS.jsonl"
    with gzip.open(taxon_path, "rt", encoding="utf8") as tf, dataset_path.open("r", encoding="utf8") as df:
        for cell in expected_cells:
            for replicate in evaluation:
                _, source = annual_record(cell, replicate)
                for method in METHODS:
                    line = tf.readline()
                    if not line:
                        fail(diffs, "TAXON_DECISIONS.eof", "row", None)
                        break
                    got_taxon = json.loads(line)
                    decisions_seen += 1
                    compare_scalar(diffs, "taxon.cell", cell, got_taxon.get("cell"))
                    compare_scalar(diffs, "taxon.replicate", replicate, got_taxon.get("replicate"))
                    compare_scalar(diffs, "taxon.method", method, got_taxon.get("method"))
                    obj = source["methods"][method]
                    a = (ar(obj, "A1", "effect") + ar(obj, "A2", "effect")) / 2.0
                    b = ar(obj, "B", "effect")
                    s1 = ar(obj, "A1", "se")
                    s2 = ar(obj, "A2", "se")
                    sb = ar(obj, "B", "se")
                    ta = ar(obj, "A1", "truth")
                    tb = ar(obj, "B", "truth")
                    z = sb * sb
                    d = a - b
                    sign = a * b > 0.0
                    rho, lower, upper, available_cal = calibration[(cell, method)]
                    domain_bounds = {
                        "zero": (np.zeros(100), np.zeros(100), np.ones(100, dtype=bool)),
                        "calibrated_point": (rho, rho, available_cal),
                        "calibrated_interval": (lower, upper, available_cal),
                        "full": (np.full(100, -1.0), np.full(100, 1.0), np.ones(100, dtype=bool)),
                    }
                    expected_req = {}
                    expected_sel = {}
                    for domain, (rlo, rhi, availability) in domain_bounds.items():
                        req = np.full(100, np.nan)
                        indices = np.flatnonzero(availability)
                        if indices.size:
                            req[indices] = required_cap_by_containment(
                                rlo[indices], rhi[indices], s1[indices], s2[indices], z[indices], d[indices]
                            )
                        expected_req[domain] = req
                        # Membership is independently classified by direct window containment.
                        # The bisection cap is audited as a separate continuous quantity.
                        expected_sel[domain] = [
                            sign & availability & contained_in_window(
                                rlo, rhi, s1, s2, z, d, cap
                            )
                            for cap in CAPS
                        ]
                        compare_array(
                            diffs, f"taxon.required_caps.{cell}.{replicate}.{method}.{domain}",
                            [None if not np.isfinite(v) else float(v) for v in req],
                            got_taxon["required_caps"][domain],
                        )
                        for index, cap in enumerate(CAPS):
                            compare_array(
                                diffs, f"taxon.selected.{cell}.{replicate}.{method}.{domain}.{cap}",
                                expected_sel[domain][index].tolist(), got_taxon["selected"][domain][index],
                            )
                    reqiv = required_cap_by_containment(
                        np.full(100, -1.0), np.full(100, 1.0), s1, s2, z, np.zeros(100)
                    )
                    compare_array(diffs, f"taxon.required_cap_IV_full.{cell}.{replicate}.{method}", reqiv.tolist(), got_taxon["required_cap_IV_full"])
                    compare_array(diffs, f"taxon.same_estimated_sign.{cell}.{replicate}.{method}", sign.tolist(), got_taxon["same_estimated_sign"])
                    compare_array(
                        diffs, f"taxon.nonregular.{cell}.{replicate}.{method}",
                        (x_of_rho(-1.0, s1, s2) == 0.0).tolist(), got_taxon["nonregular_full_endpoint"],
                    )
                    compare_array(diffs, f"taxon.truth_A.{cell}.{replicate}.{method}", ta.tolist(), got_taxon["truth_A"])
                    compare_array(diffs, f"taxon.truth_B.{cell}.{replicate}.{method}", tb.tolist(), got_taxon["truth_B"])
                    alpha0 = x_of_rho(0.0, s1, s2)
                    snr = np.minimum(np.abs(a) / np.sqrt(alpha0), np.abs(b) / sb)
                    for qi, cap in enumerate(CAPS):
                        baseline = expected_sel["zero"][qi]
                        full = expected_sel["full"][qi]
                        for domain, (_, _, availability) in domain_bounds.items():
                            expected_m, expected_d = dataset_metrics(
                                expected_sel[domain][qi], ta, tb, a, baseline, full, snr, reqiv, availability
                            )
                            req = expected_req[domain]
                            expected_m["mean_required_cap"] = ratio(req, availability)
                            expected_m["mean_required_cap_IV_full"] = float(reqiv.mean())
                            expected_m["heterogeneity_assisted_membership"] = (
                                float((expected_sel[domain][qi] & (reqiv > cap)).mean()) if domain == "full" else None
                            )
                            row_line = df.readline()
                            if not row_line:
                                fail(diffs, "DATASET_RECORDS.eof", "row", None)
                                continue
                            got_row = json.loads(row_line)
                            dataset_seen += 1
                            for key, value in (("cell", cell), ("replicate", replicate), ("method", method), ("q", cap), ("domain", domain)):
                                compare_scalar(diffs, f"dataset.identity.{key}", value, got_row.get(key))
                            for name, value in expected_m.items():
                                compare_scalar(diffs, f"dataset.metrics.{cell}.{replicate}.{method}.{cap}.{domain}.{name}", value, got_row["metrics"].get(name))
                                grouped[(cell, method, cap, domain)][name].add(value)
                            for name, value in expected_d.items():
                                compare_scalar(diffs, f"dataset.denominators.{cell}.{replicate}.{method}.{cap}.{domain}.{name}", value, got_row["denominators"].get(name))
                                den_hist[(cell, method, cap, domain)][name][str(value)] += 1
                        # Width sensitivity output follows each declared cap after the four core domains.
                        for width in HALFWIDTHS:
                            rlo = np.maximum(-1.0, rho - width)
                            rhi = np.minimum(1.0, rho + width)
                            avail = available_cal
                            req = np.full(100, np.nan)
                            indices = np.flatnonzero(avail)
                            if indices.size:
                                req[indices] = required_cap_by_containment(
                                    rlo[indices], rhi[indices], s1[indices], s2[indices], z[indices], d[indices]
                                )
                            selected = sign & avail & contained_in_window(rlo, rhi, s1, s2, z, d, cap)
                            same, opposite, _, _ = truth_masks(ta, tb)
                            expected_m = {
                                "membership": ratio(selected, avail),
                                "true_same_retention": ratio(selected, same & avail),
                                "true_opposite_membership": ratio(selected, opposite & avail),
                                "mean_required_cap": ratio(req, avail),
                                "mean_domain_width": ratio(rhi - rlo, avail),
                            }
                            row_line = df.readline()
                            if not row_line:
                                fail(diffs, "DATASET_RECORDS.width_eof", "row", None)
                                continue
                            got_row = json.loads(row_line)
                            dataset_seen += 1
                            domain = "width_" + str(width)
                            for key, value in (("cell", cell), ("replicate", replicate), ("method", method), ("q", cap), ("domain", domain)):
                                compare_scalar(diffs, f"width.identity.{key}", value, got_row.get(key))
                            expected_d = {"all_taxa": 100, "available": int(avail.sum()), "unavailable": int((~avail).sum())}
                            for name, value in expected_m.items():
                                compare_scalar(diffs, f"width.metrics.{cell}.{replicate}.{method}.{cap}.{domain}.{name}", value, got_row["metrics"].get(name))
                                grouped[(cell, method, cap, domain)][name].add(value)
                            for name, value in expected_d.items():
                                compare_scalar(diffs, f"width.denominators.{cell}.{replicate}.{method}.{cap}.{domain}.{name}", value, got_row["denominators"].get(name))
                    # Rebuild deterministic illustration selection, including the required calibration availability.
                    same_truth, _, _, _ = truth_masks(ta, tb)
                    for label, test in (
                        ("heterogeneity_assisted", reqiv > 0.8),
                        ("inverse_variance_sufficient", reqiv <= 0.8),
                    ):
                        valid = same_truth & sign & available_cal & (expected_req["full"] > 0.6) & (expected_req["full"] <= 0.8) & test
                        key = method + "_" + label
                        if key not in expected_cases and np.any(valid):
                            j = int(np.flatnonzero(valid)[0])
                            expected_cases[key] = {
                                "cell": cell, "replicate": replicate, "method": method, "taxon": j,
                                "parameters": manifest["cells"][cell], "a": float(a[j]), "b": float(b[j]),
                                "s1": float(s1[j]), "s2": float(s2[j]), "sb": float(sb[j]),
                                "truth_A": float(ta[j]), "truth_B": float(tb[j]),
                                "qstar_REML": float(expected_req["full"][j]), "qstar_IV": float(reqiv[j]),
                                "rho_hat": float(rho[j]), "calibration_interval": [float(lower[j]), float(upper[j])],
                                "selection": "first record under fixed ordering satisfying predeclared illustrative condition",
                            }
        if tf.readline():
            fail(diffs, "TAXON_DECISIONS.trailing_rows", False, True)
        if df.readline():
            fail(diffs, "DATASET_RECORDS.trailing_rows", False, True)

    # Independent summary and deterministic-case comparisons.
    for key, stats in grouped.items():
        cell, method, cap, domain = key
        got = stored_summary.get(key)
        if got is None:
            fail(diffs, f"SUMMARY[{key}]", "row", None)
            continue
        compare_scalar(diffs, f"SUMMARY[{key}].cell", cell, got.get("cell"))
        compare_scalar(diffs, f"SUMMARY[{key}].method", method, got.get("method"))
        compare_scalar(diffs, f"SUMMARY[{key}].q", cap, got.get("q"))
        compare_scalar(diffs, f"SUMMARY[{key}].domain", domain, got.get("domain"))
        if got.get("parameters") != manifest["cells"][cell]:
            fail(diffs, f"SUMMARY[{key}].parameters", manifest["cells"][cell], got.get("parameters"))
        compare_scalar(
            diffs, f"SUMMARY[{key}].uncertainty",
            "dataset MCSE conditional on fixed calibration sample", got.get("uncertainty"),
        )
        for name, moment in stats.items():
            expected = moment.payload()
            for field, value in expected.items():
                compare_scalar(diffs, f"SUMMARY[{key}].metrics.{name}.{field}", value, got["metrics"].get(name, {}).get(field))
        want_hist = {name: dict(values) for name, values in den_hist[key].items()}
        if want_hist != got.get("denominator_histograms"):
            fail(diffs, f"SUMMARY[{key}].denominator_histograms", want_hist, got.get("denominator_histograms"))
    if set(expected_cases) != set(stored_cases):
        fail(diffs, "ILLUSTRATIVE_CASES.keys", sorted(expected_cases), sorted(stored_cases))
    for key, expected in expected_cases.items():
        got = stored_cases.get(key, {})
        for field, value in expected.items():
            if isinstance(value, list):
                compare_array(diffs, f"ILLUSTRATIVE_CASES.{key}.{field}", value, got.get(field, []))
            elif isinstance(value, dict):
                if value != got.get(field):
                    fail(diffs, f"ILLUSTRATIVE_CASES.{key}.{field}", value, got.get(field))
            else:
                compare_scalar(diffs, f"ILLUSTRATIVE_CASES.{key}.{field}", value, got.get(field))

    try:
        math_checks = edge_math_checks()
        math_error = None
    except Exception as error:  # Preserve a real math-check failure as audit evidence.
        math_checks = 0
        math_error = repr(error)
        fail(diffs, "independent_math_edge_checks", "pass", math_error)
    payload = {
        "scale": scale,
        "status": "PASS" if not diffs else "FAIL",
        "differences": len(diffs),
        "examples": diffs[:100],
        "tolerance": {"atol": ATOL, "rtol": RTOL, "independent_likelihood": 1e-6},
        "source_hashes_checked": accessed,
        "source_hashes_unchanged": not any(d["location"].startswith(("source_hashes", "post_audit_seal")) for d in diffs),
        "producer_imported": False,
        "prior_certificate_or_analysis_imported": False,
        "new_fit_or_resampling": False,
        "calibration_reconstructed": True,
        "classification_reconstructed_by_inverse_window": True,
        "dataset_MCSE_reconstructed": True,
        "math_edge_checks": math_checks,
        "math_edge_error": math_error,
        "taxon_decision_rows_checked": decisions_seen,
        "dataset_summary_rows_checked": dataset_seen,
        "original_native_exit_code": None,
    }
    target = out / "AUDIT_RESULT.json"
    with target.open("x", encoding="utf8") as handle:
        json.dump(payload, handle, indent=2, allow_nan=False)
    return 0 if not diffs else 1


if __name__ == "__main__":
    raise SystemExit(audit(sys.argv[1]))
