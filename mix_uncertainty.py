#!/usr/bin/env python3
"""
mix_uncertainty.py — First-order uncertainty of mixing ratios computed with MIX.

For each mixture p:

    C_p      = sum_e lambda_pe^2 * C_ms,e + C_ms,p          (1)
    lambda_p = argmin (c_p - E lambda_p)^T C_p^-1 (c_p - E lambda_p)   (2)
    C_lambda = (E^T C_p^-1 E)^-1                             (3)

End-members (EMs) are treated as fixed (their identity is known) but uncertain
(their concentrations have errors). Input variances are rescaled by a
species-specific factor f_s (Carrera and Neuman 1986, Eq. 6), computed so that
the mean squared normalised residual of each species equals 1.

The script reads the MIX output file (.res) directly and, for every mixture,
reports the mixing ratios, their standard deviations, the +/- k*sigma interval
(clipped to [0, 1]), the correlation matrix and the condition number of the
weighted mixing matrix inverted in (3).

Two modes are available:

  * default ("MIX"): E = EM concentrations estimated by MIX (OPTIMAL SOLUTION),
    lambda = optimal mixing ratios; f_s from the residuals of EMs and mixtures.
  * --fixed-em ("EMMA"): E = measured (input) EM concentrations, lambda = the
    constrained least-squares ratios printed by MIX as INITIAL; f_s from the
    residuals of the mixtures only. This corresponds to EMMA / Gaussian
    mixing models with fixed EMs.

Usage
-----
    python mix_uncertainty.py run.res
    python mix_uncertainty.py run.res --k 2 --out sigma_4EM.csv
    python mix_uncertainty.py run.res --fixed-em --out sigma_EMMA.csv

Requirements: Python >= 3.8, numpy.

Limitations
-----------
* Diagonal covariance matrices (as in MIX input files).
* Mixing ratios fixed by constraints (e.g. rainwater set to zero in some
  samples) are treated as estimated; remove those EMs from the analysis of
  those samples if exact treatment is required.
* First-order (linearised) approximation evaluated at the MIX solution.
  The prior variances of the EMs are included (C_p), but not the posterior
  uncertainty of the EM compositions estimated by MIX, nor conceptual
  (model-structure) uncertainty.
* The sum-to-one and non-negativity constraints are not imposed in the
  covariance; intervals are conservative and clipped to [0, 1].
"""

import argparse
import csv
import re
import sys

import numpy as np

NUM = re.compile(r"^[+-]?(\d+\.?\d*|\.\d+)([eEdD][+-]?\d+)?$")


def is_num(tok):
    return bool(NUM.match(tok))


def to_float(tok):
    return float(tok.replace("D", "E").replace("d", "e"))


# --------------------------------------------------------------------------
# Parsing of the MIX .res file
# --------------------------------------------------------------------------
def read_rows(lines, start, n, ncol):
    """Read n rows 'name v1 ... vncol' starting at line index `start`
    (blank or separator lines are skipped). Returns names, array, next index."""
    names, vals, i = [], [], start
    while len(names) < n and i < len(lines):
        tok = lines[i].split()
        i += 1
        if not tok or set(tok[0]) <= set("-+=:"):
            continue
        nums = [t for t in tok[1:] if is_num(t)]
        if len(nums) < ncol:
            continue
        names.append(tok[0])
        vals.append([to_float(t) for t in nums[:ncol]])
    return names, np.array(vals, float), i


def find(lines, pattern, start=0):
    rx = re.compile(pattern, re.IGNORECASE)
    for i in range(start, len(lines)):
        if rx.search(lines[i]):
            return i
    raise ValueError(f"Pattern not found in .res file: {pattern!r}")


def parse_res(path):
    with open(path, encoding="latin-1") as fh:
        lines = fh.read().splitlines()

    ns = int(lines[find(lines, r"Number\s+of\s+species")].split()[-1])
    nw = int(lines[find(lines, r"Number\s+of\s+wells")].split()[-1])
    ne = int(lines[find(lines, r"Number\s+of\s+sources")].split()[-1])

    # measured concentrations
    i = find(lines, r"^\s*CONCENTRATION\s+DATA")
    species = lines[find(lines, r"^\s*Species", i)].split()[1:1 + ns]
    i = find(lines, r"Sources\s*\(end-members\)", i)
    em_names, E_in, i = read_rows(lines, i + 1, ne, ns)
    i = find(lines, r"Wells\s*\(mixed", i)
    w_names, C, i = read_rows(lines, i + 1, nw, ns)

    # input variances
    i = find(lines, r"VARIANCES\s+OF\s+CONCENTRATION", i)
    i = find(lines, r"Sources\s*\(end-members\)", i)
    _, V_em, i = read_rows(lines, i + 1, ne, ns)
    i = find(lines, r"Wells\s*\(mixed", i)
    _, V_w, i = read_rows(lines, i + 1, nw, ns)

    # initial (EMMA, fixed EMs) mixing ratios
    i = find(lines, r"CONSTRAINED\s+LEAST\s+SQUARES", i)
    lam_ini = []
    while len(lam_ini) < nw:
        i += 1
        tok = lines[i].split()
        if tok and tok[0].lower().startswith("well"):
            nums = [to_float(t) for t in tok if is_num(t)]
            lam_ini.append(nums[1:1 + ne])          # first number = well index
    lam_ini = np.array(lam_ini, float)

    # objective function of each attempt (MIX may converge to local optima)
    attempts, current = [], None
    for k in range(i, len(lines)):
        if re.search(r"Attempt\s+number", lines[k], re.IGNORECASE):
            current = int(lines[k].split()[-1])
        elif re.search(r"OPTIMAL\s+SOLUTION", lines[k], re.IGNORECASE):
            break
        elif current is not None and re.search(r"objective\s+function:", lines[k], re.IGNORECASE):
            attempts.append((current, to_float(lines[k].split()[-1])))
            current = None

    # optimal solution
    i = find(lines, r"OPTIMAL\s+SOLUTION", i)
    fo = to_float(lines[find(lines, r"objective\s+function:", i)].split()[-1])
    j = find(lines, r"Expected\s+value\s+of\s+sources", i)
    _, E_opt, j = read_rows(lines, j + 1, ne, ns)
    j = find(lines, r"Mixing\s+ratios", j)
    _, lam_opt, _ = read_rows(lines, j + 1, nw, ne)

    return dict(species=species, em=em_names, wells=w_names, E_in=E_in, C=C,
                V_em=V_em, V_w=V_w, lam_ini=lam_ini, E_opt=E_opt, lam_opt=lam_opt,
                objective=fo, attempts=attempts)


# --------------------------------------------------------------------------
# Uncertainty (Eqs. 1-3)
# --------------------------------------------------------------------------
def scaling_factors(d, fixed_em):
    """f_s: mean squared normalised residual per species (Carrera & Neuman 1986)."""
    if fixed_em:
        E, lam = d["E_in"], d["lam_ini"]
        r2 = (d["C"] - lam @ E) ** 2 / d["V_w"]
    else:
        E, lam = d["E_opt"], d["lam_opt"]
        r_w = (d["C"] - lam @ E) ** 2 / d["V_w"]
        r_e = (E - d["E_in"]) ** 2 / d["V_em"]
        r2 = np.vstack([r_e, r_w])
    return r2.mean(axis=0)


def uncertainty(d, fixed_em=False):
    E = d["E_in"] if fixed_em else d["E_opt"]
    lam = d["lam_ini"] if fixed_em else d["lam_opt"]
    fs = scaling_factors(d, fixed_em)
    out = []
    for p in range(len(d["wells"])):
        cp = fs * d["V_w"][p] + (lam[p] ** 2) @ (fs * d["V_em"])     # Eq. (1), diagonal
        A = E / np.sqrt(cp)                                          # weighted E (ne x ns)
        M = A @ A.T                                                  # E C_p^-1 E^T
        cov = np.linalg.inv(M)                                       # Eq. (3)
        sig = np.sqrt(np.diag(cov))
        corr = cov / np.outer(sig, sig)
        out.append(dict(well=d["wells"][p], lam=lam[p], sigma=sig, corr=corr,
                        cond=np.linalg.cond(A.T)))
    return out, fs


# --------------------------------------------------------------------------
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("res", help="MIX output file (.res)")
    ap.add_argument("--k", type=float, default=2.0,
                    help="interval half-width in standard deviations (default 2)")
    ap.add_argument("--fixed-em", action="store_true",
                    help="EMMA mode: measured EMs and initial (CLS) mixing ratios")
    ap.add_argument("--out", help="write results to this CSV file")
    ap.add_argument("--corr", action="store_true", help="print correlation matrices")
    a = ap.parse_args(argv)

    d = parse_res(a.res)
    res, fs = uncertainty(d, a.fixed_em)
    em = d["em"]
    mode = "EMMA (fixed EMs)" if a.fixed_em else "MIX (estimated EMs)"

    print(f"File: {a.res}   Mode: {mode}   k = {a.k:g}")
    if d["attempts"]:
        print(f"Objective function (maximised): optimal = {d['objective']:g}; attempts: "
              + ", ".join(f"#{n} {v:g}" for n, v in d["attempts"]))
    print("f_s: " + "  ".join(f"{s}={v:.4g}" for s, v in zip(d["species"], fs)))
    print()
    head = f"{'Mixture':12s} " + " ".join(f"{e:>20s}" for e in em) + f" {'cond':>7s}"
    print(head)
    print("-" * len(head))
    rows = []
    for r in res:
        cells = []
        for e, l, s in zip(em, r["lam"], r["sigma"]):
            lo, hi = max(0.0, l - a.k * s), min(1.0, l + a.k * s)
            cells.append(f"{l:5.3f} ± {a.k * s:5.3f}".rjust(20))
            rows.append([r["well"], e, l, s, lo, hi, r["cond"]])
        print(f"{r['well']:12s} " + " ".join(cells) + f" {r['cond']:7.1f}")
        if a.corr:
            print("  corr:\n" + "\n".join("   " + " ".join(f"{x:6.2f}" for x in row)
                                        for row in r["corr"]))
    print("\ncond = condition number of the weighted mixing matrix C_p^-1/2 E^T")

    if a.out:
        with open(a.out, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["mixture", "end_member", "lambda", "sigma",
                        f"lower_{a.k:g}sigma", f"upper_{a.k:g}sigma", "condition_number"])
            w.writerows(rows)
        print(f"Results written to {a.out}")


if __name__ == "__main__":
    sys.exit(main())
