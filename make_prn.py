#!/usr/bin/env python3
"""
make_prn.py — Build the input file (.prn) for MIX / EMMA-MIX (Carrera et al., 2004)
from a CSV table of concentrations.

Authors: Jesús Carrera, Romina Marcovecchio, Sebastián Dietrich.
Python translation and extension, by R. Marcovecchio, of the R script written by
S. Dietrich (IHLLA-CONICET) to build MIX input files. The layout follows the
MIX Guide Notes (Vázquez-Suñé and Carrera) and reproduces, field by field, the
.prn written by the R script, which MIX reads without errors.

Input CSV
---------
* One row per water: end-members first, then mixtures. MIX keeps this order,
  so the rows of the CSV, the .prn and the .res match. The script asks which
  rows are end-members (or use --sources / --em).
* One column per species, in the order MIX will use them; the LAST column is
  the water name (no spaces; keep it short, <= 10 characters recommended).
* Header row with species names. Names starting with a digit (e.g. "2H") are
  renamed as R does ("X2H").
* Units are not converted: use mmol/L for ions, µS/cm for EC and ‰ for
  isotopes (or any consistent set of units).
* Field separator ";" by default (use --sep "," otherwise). Decimal comma is
  accepted with --decimal ",".

Usage
-----
    python make_prn.py data.csv                 # interactive: asks which rows are end-members
    python make_prn.py data.csv --sources 4     # first 4 rows are end-members
    python make_prn.py data.csv --em 1-3 --exclude 4   # 3 end-members, row 4 left out
    python make_prn.py data.csv --sources 4 --out run4EM.prn --attempts 9
    python make_prn.py data.csv --sources 4 --constraint 5,2,0 --constraint 6,2,0
    python make_prn.py data.csv --sources 4 --em-var-factor 10

Then run MIX with the root name of the .prn file (without extension).

Variances
---------
As in the R script, the variance of every concentration is set equal to the
squared concentration (coefficient of variation of 100 %) for end-members and
mixtures alike. MIX rescales variances by a species-specific factor, so only
their relative values matter. --em-var-factor multiplies the end-member
variances (e.g. 10 makes end-members ~3 times more uncertain, in standard
deviation, than the mixtures). --var-csv gives explicit variances for selected
waters (e.g. a better characterised end-member).

Requirements: Python >= 3.8, pandas, numpy.
"""

import argparse
import re
import sys

import numpy as np
import pandas as pd

HEADER_LABELS = ["Num_spe", "Num_well", "Num_sour", "nº_inten", "iflag", "nseed"]
WIDTH_HEAD = 10      # width of header fields (R: format(width = 10))
WIDTH = 14           # width of data fields   (R: format(width = 14))
NEWLINE = "\r\n"     # R on Windows writes CRLF; MIX was built for Windows


# --------------------------------------------------------------------------
# Number formatting identical to R's format(x) with default options
# (7 significant digits, fixed notation unless scientific is shorter)
# --------------------------------------------------------------------------
def r_format(x, digits=7):
    x = float(x)
    if x == 0:
        return "0"
    if not np.isfinite(x):
        raise ValueError(f"Non-finite value: {x}")
    sci = f"{x:.{digits - 1}e}"
    mant, exp = sci.split("e")
    mant = mant.rstrip("0").rstrip(".")
    e = int(exp)
    sci_str = f"{mant}e{'-' if e < 0 else '+'}{abs(e):02d}"
    nsig = len(mant.replace("-", "").replace(".", ""))
    ndec = max(nsig - 1 - e, 0)
    fixed_str = f"{x:.{ndec}f}"
    return fixed_str if len(fixed_str) <= len(sci_str) else sci_str


def r_make_name(name):
    """Mimic R's make.names() for column headers (e.g. '2H' -> 'X2H')."""
    name = re.sub(r"[^A-Za-z0-9._]", ".", str(name).strip())
    if not name or not (name[0].isalpha() or (name[0] == "." and not name[1:2].isdigit())):
        name = "X" + name
    return name


def row_line(values, name=None):
    cells = [r_format(v).rjust(WIDTH) for v in values]
    line = " ".join(cells)
    if name is not None:
        line += " " + str(name).ljust(WIDTH)
    return line


# --------------------------------------------------------------------------
def read_table(path, sep, decimal):
    df = pd.read_csv(path, sep=sep, decimal=decimal, skipinitialspace=True)
    df.columns = [c.strip() for c in df.columns]
    names = df.iloc[:, -1].astype(str).str.strip()
    data = df.iloc[:, :-1]
    bad = [c for c in data.columns if not pd.api.types.is_numeric_dtype(data[c])]
    if bad:
        raise SystemExit(f"ERROR: non-numeric values in column(s) {bad}. "
                         "Check the separator (--sep) and decimal mark (--decimal).")
    if data.isna().any().any():
        r, c = np.argwhere(data.isna().values)[0]
        raise SystemExit(f"ERROR: empty value in row '{names[r]}', column '{data.columns[c]}'.")
    if names.str.contains(r"\s").any():
        raise SystemExit("ERROR: water names must not contain spaces: "
                         + ", ".join(names[names.str.contains(r'\s')]))
    if names.duplicated().any():
        raise SystemExit("ERROR: duplicated water names: " + ", ".join(names[names.duplicated()]))
    long_names = names[names.str.len() > 10]
    if len(long_names):
        print("WARNING: names longer than 10 characters may be truncated by MIX: "
              + ", ".join(long_names), file=sys.stderr)
    species = [r_make_name(c) for c in data.columns]
    return species, data.values.astype(float), list(names)


def initial_ratios(n_well, n_sour, mode, rng):
    if mode == "equal":
        return np.full((n_well, n_sour), 1.0 / n_sour)
    p = rng.random((n_well, n_sour))
    return p / p.sum(axis=1, keepdims=True)


def parse_constraint(text, n_well, n_sour):
    try:
        w, s, p = text.split(",")
        w, s, p = int(w), int(s), float(p)
    except ValueError:
        raise SystemExit(f"ERROR: constraint '{text}' must be well,source,proportion (e.g. 5,2,0)")
    if not (1 <= w <= n_well and 1 <= s <= n_sour and 0 <= p <= 1):
        raise SystemExit(f"ERROR: constraint '{text}' out of range "
                         f"(wells 1-{n_well}, sources 1-{n_sour}, proportion 0-1)")
    return w, s, p


def build_prn(species, X, names, n_sour, attempts=9, iflag=0, nseed=18201,
              em_var_factor=1.0, constraints=(), init="random", rng=None, var_override=None):
    n_tot, n_spe = X.shape
    n_well = n_tot - n_sour
    if n_sour < 2 or n_well < 1:
        raise SystemExit("ERROR: need at least 2 sources and 1 mixture.")
    if n_spe < n_sour:
        print("WARNING: fewer species than end-members; mixing ratios will be poorly "
              "determined.", file=sys.stderr)

    V = X ** 2                                   # variance = concentration^2
    V[:n_sour] *= em_var_factor
    for name, row in (var_override or {}).items():   # user-given variances (by water name)
        V[names.index(name)] = row
    lam0 = initial_ratios(n_well, n_sour, init, rng or np.random.default_rng())

    out = []
    out.append(" ".join(l.rjust(WIDTH_HEAD) for l in HEADER_LABELS) + " ")
    out.append(" ".join(str(v).rjust(WIDTH_HEAD)
                        for v in (n_spe, n_well, n_sour, attempts, iflag, nseed)) + " ")
    out.append("First sources and then wells ")
    out.append("(CONC(i,j),j=1,Num_spe) nam " + ",".join(species) + " ")
    out += [row_line(X[i], names[i]) for i in range(n_tot)]
    out.append("Covariances ")
    out += [row_line(V[i], names[i]) for i in range(n_tot)]
    # as the R script: "0 " with no constraints; right-aligned count otherwise
    out.append(str(len(constraints)).rjust(WIDTH) if constraints else "0 ")
    out += [row_line((w, s, p)) + " " for w, s, p in constraints]
    out += [row_line(lam0[i]) for i in range(n_well)]
    return NEWLINE.join(out) + NEWLINE


# --------------------------------------------------------------------------
# Interactive questions (used when neither --sources nor --em is given)
# --------------------------------------------------------------------------
def ask(prompt):
    try:
        return input(prompt).strip()
    except EOFError:
        raise SystemExit("\nERROR: no answer received; use --sources N or --em LIST to run non-interactively.")


def parse_selection(txt, n):
    """'1-4' or '1,3,5' or '1-3,7' -> sorted list of 0-based indices."""
    idx = set()
    for part in txt.replace(" ", "").split(","):
        if not part:
            continue
        if "-" in part:
            a, b = part.split("-", 1)
            idx.update(range(int(a), int(b) + 1))
        else:
            idx.add(int(part))
    if not idx or min(idx) < 1 or max(idx) > n:
        raise ValueError
    return sorted(i - 1 for i in idx)


def ask_sources(names):
    """Ask which waters are end-members and which to leave out (0-based indices)."""
    n = len(names)
    print(f"\nThe CSV contains {n} waters:")
    for i, nm in enumerate(names, 1):
        print(f"  {i:3d}  {nm}")
    while True:
        txt = ask("\nWhich rows are end-members? Enter their numbers (e.g. 1-4 or 1,2,3,4): ")
        try:
            sel = parse_selection(txt, n)
        except ValueError:
            print(f"  Please enter row numbers between 1 and {n}, separated by commas or as a range.")
            continue
        txt = ask("Rows to leave out of this run (e.g. 4), or Enter to use all: ")
        try:
            excl = parse_selection(txt, n) if txt else []
        except ValueError:
            print(f"  Please enter row numbers between 1 and {n}.")
            continue
        if set(excl) & set(sel):
            print("  A row cannot be an end-member and be left out at the same time.")
            continue
        rest = [i for i in range(n) if i not in sel and i not in excl]
        if len(sel) < 2 or not rest:
            print("  Choose at least 2 end-members and keep at least 1 mixture.")
            continue
        print(f"  End-members ({len(sel)}): {', '.join(names[i] for i in sel)}")
        print(f"  Mixtures    ({len(rest)}): {', '.join(names[i] for i in rest)}")
        if excl:
            print(f"  Left out    ({len(excl)}): {', '.join(names[i] for i in excl)}")
        if ask("Is this correct? [y/n]: ").lower() in ("y", "yes", "s", "si", "sí"):
            return sel, excl


def ask_constraints(names, k):
    wells, sources = names[k:], names[:k]
    print("\nOptional: fix the proportion of an end-member in a mixture "
          "(e.g. rainwater = 0 in a baseflow sample).")
    print("  Wells (mixtures):     " + "  ".join(f"{i}={w}" for i, w in enumerate(wells, 1)))
    print("  Sources (end-members): " + "  ".join(f"{j}={s}" for j, s in enumerate(sources, 1)))
    cons = []
    while True:
        txt = ask("Constraint as well,source,proportion (e.g. 5,2,0), or Enter to finish: ")
        if not txt:
            return cons
        try:
            cons.append(parse_constraint(txt, len(wells), k))
            w, s, p = cons[-1]
            print(f"  -> {sources[s - 1]} fixed at {p:g} in {wells[w - 1]}")
        except SystemExit as e:
            print("  " + str(e))


def main(argv=None):
    ap = argparse.ArgumentParser(description="Build a MIX .prn file from a CSV table.")
    ap.add_argument("csv", help="input CSV: one row per water, one column per species; last column = name")
    ap.add_argument("--sources", type=int,
                    help="number of end-members, taken as the first rows; if omitted, the script asks")
    ap.add_argument("--em", help='rows that are end-members, e.g. "1-4" or "1,3,6" (any position)')
    ap.add_argument("--exclude", help='rows to leave out of this run, e.g. "4"')
    ap.add_argument("--out", help="output .prn (default: same name as the CSV)")
    ap.add_argument("--attempts", type=int, default=9, help="number of MIX attempts (default 9)")
    ap.add_argument("--nseed", type=int, default=18201, help="MIX random seed (default 18201)")
    ap.add_argument("--sep", default=";", help='CSV field separator (default ";")')
    ap.add_argument("--decimal", default=".", help='decimal mark (default "."; use "," if needed)')
    ap.add_argument("--em-var-factor", type=float, default=1.0,
                    help="multiply end-member variances by this factor (default 1)")
    ap.add_argument("--constraint", action="append", default=[],
                    help="fix a proportion: well,source,proportion, numbered as in MIX (1-based; repeatable)")
    ap.add_argument("--var-csv",
                    help="CSV with user-defined variances for some waters (same columns as the "
                         "data CSV; rows matched by name). Other waters keep variance = conc^2")
    ap.add_argument("--init", choices=["random", "equal"], default="random",
                    help="initial mixing ratios (default random, as in the R script)")
    ap.add_argument("--rng-seed", type=int, help="seed for the random initial ratios (reproducible)")
    a = ap.parse_args(argv)

    species, X, names = read_table(a.csv, a.sep, a.decimal)
    n = len(names)
    interactive = a.sources is None and a.em is None
    excl = []
    if interactive:
        sel, excl = ask_sources(names)
    elif a.em is not None:
        try:
            sel = parse_selection(a.em, n)
        except ValueError:
            raise SystemExit(f"ERROR: --em must list row numbers between 1 and {n}.")
    else:
        sel = list(range(a.sources))
    if a.exclude:
        try:
            excl = parse_selection(a.exclude, n)
        except ValueError:
            raise SystemExit(f"ERROR: --exclude must list row numbers between 1 and {n}.")
    if set(excl) & set(sel):
        raise SystemExit("ERROR: a row cannot be an end-member and be excluded.")
    rest = [i for i in range(n) if i not in sel and i not in excl]
    if len(sel) < 2 or not rest:
        raise SystemExit("ERROR: need at least 2 end-members and 1 mixture.")
    # MIX needs end-members first: reorder rows (end-members, then mixtures, keeping order)
    order = sel + rest
    X, names = X[order], [names[i] for i in order]
    a.sources = len(sel)
    cons = ask_constraints(names, a.sources) if interactive else []
    n_well = len(names) - a.sources
    cons += [parse_constraint(c, n_well, a.sources) for c in a.constraint]
    var_override = None
    if a.var_csv:
        vsp, VX, vnames = read_table(a.var_csv, a.sep, a.decimal)
        if vsp != species:
            raise SystemExit("ERROR: the variance CSV must have the same species columns as the data CSV.")
        unknown = [v for v in vnames if v not in names]
        if unknown:
            raise SystemExit("ERROR: waters in the variance CSV not used in this run: " + ", ".join(unknown))
        if (VX <= 0).any():
            raise SystemExit("ERROR: variances must be positive.")
        var_override = dict(zip(vnames, VX))
        print("User-defined variances for: " + ", ".join(vnames))
    text = build_prn(species, X, names, a.sources, a.attempts, 0, a.nseed,
                     a.em_var_factor, cons, a.init, np.random.default_rng(a.rng_seed), var_override)

    out = a.out or re.sub(r"\.[^.]*$", "", a.csv) + ".prn"
    with open(out, "w", encoding="utf-8", newline="") as fh:
        fh.write(text)

    print(f"Written {out}")
    print(f"  species ({len(species)}): {', '.join(species)}")
    print(f"  sources ({a.sources}): {', '.join(names[:a.sources])}")
    print(f"  mixtures ({n_well}): {', '.join(names[a.sources:])}")
    print(f"  constraints: {len(cons)}   attempts: {a.attempts}   EM variance factor: {a.em_var_factor:g}")


if __name__ == "__main__":
    sys.exit(main())
