# EMMA-MIX workflow: from a CSV file to mixing ratios and their uncertainty

```
data.csv ──► make_prn.py ──► root.prn ──► EMMA-MIX (console) ──► root.res ──► mix_uncertainty.py ──► mixing ratios ± σ
```

| Step | Tool | Input | Output |
|---|---|---|---|
| 1 | `make_prn.py` | `data.csv` | `root.prn` |
| 2 | EMMA-MIX (MIX program, Carrera et al., 2004) | `root.prn` | `root.res`, `rootEMMA.out`, `rooteigen.prn` |
| 3 | `mix_uncertainty.py` | `root.res` | uncertainty of the mixing ratios (screen and `.csv`) |

## Requirements

- Python ≥ 3.8 with `numpy` and `pandas`:
  ```
  pip install numpy pandas
  ```
- The EMMA-MIX (MIX) executable.

Commands are typed in a terminal (e.g. PowerShell, or **Terminal → New Terminal** in VS Code) opened in the folder that contains the scripts and the data. Type one command at a time and press **Enter**.

---

## Step 1. Build the `.prn` file

### Input file: `data.csv`

One row per water (end-members and mixtures), one column per species, and the **name of each water in the last column**.

```
Ca;Mg;Na;K;Cl;NO3;SO4;TIC;F;CE;2H;O18;id
1.7186;1.6901;2.5490;0.3754;0.4545;1.4998;0.1230;7.5622;0.0542;878;-30.3;-5.47;EM1
0.1372;0.0247;0.2653;0.1176;0.0282;0.0290;0.0343;0.4916;0.0053;67;-36.7;-6.24;EM2
...
0.5769;0.5086;1.0492;0.4230;0.1412;0.0371;0.0135;3.3362;0.0230;337;-35.1;-5.88;Sample1
...
```

- **Separator** `;` and **decimal point** `.`. For other formats use `--sep ","` and/or `--decimal ","`.
- **Units** are not converted; use a consistent set (e.g. ions in mmol/L, electrical conductivity in µS/cm, δ²H and δ¹⁸O in ‰ with their sign).
- **Names**: no spaces, unique, preferably ≤ 10 characters.
- **No empty cells.**
- End-members can be in any rows; you choose them when running the script.

### Running `make_prn.py`

There is a single script, `make_prn.py`, which can be used in two ways.

**a) Interactive: only the CSV file is given, and the script asks**

```
python make_prn.py data.csv
```

The script lists the waters and asks:

1. Which rows are end-members (e.g. `1-4`, or `1,2,4`).
2. Which rows to leave out of the run (Enter = none).
3. Whether the selection is correct (`y`/`n`).
4. Optional constraints: the proportion of an end-member fixed in a mixture, written as `mixture,end-member,proportion` (e.g. `5,2,0`). Enter = no constraints.

**b) Non-interactive: the answers are given as options, and the script asks nothing**

```
python make_prn.py data.csv --sources 4
python make_prn.py data.csv --em 1-3 --exclude 4 --out run3EM.prn
python make_prn.py data.csv --sources 4 --constraint 5,2,0 --out run4EM_c.prn
```

| Option | Meaning |
|---|---|
| `--sources N` | the first N rows are end-members |
| `--em LIST` | rows that are end-members, in any position (`1-4`, `1,3,6`) |
| `--exclude LIST` | rows left out of the run |
| `--constraint m,e,p` | fix proportion `p` of end-member `e` in mixture `m` (repeatable) |
| `--var-csv FILE` | user-defined variances for some waters (same columns as `data.csv`, rows matched by name) |
| `--em-var-factor F` | multiply the variances of all end-members by `F` |
| `--attempts N` | number of MIX attempts (default 9) |
| `--rng-seed N` | fixed seed for the initial mixing ratios (identical `.prn` every time) |
| `--out FILE` | name of the `.prn` (default: same name as the CSV, e.g. `data.prn`) |

`python make_prn.py --help` lists all the options.

### What the `.prn` contains

The layout required by MIX (MIX Guide Notes): dimensions; concentrations (end-members first, then mixtures); variances; constraints; initial mixing ratios.

- **Variances:** by default, the variance of each concentration equals its squared value, for end-members and mixtures. MIX rescales them by a factor per species, so only their relative values matter.
- **Initial mixing ratios:** random values, each row summing to 1. They do not change the final solution, because MIX tries several starting points.

---

## Step 2. Run EMMA-MIX

1. Put the **MIX executable and the `.prn` file in the same folder**.
2. Run MIX from the console (or double-click on it).
3. When it asks for the **root name**, type the name of the `.prn` file **without the extension** (e.g. `data`) and press Enter.
4. The output files are written in the same folder. Previous outputs with the same root name are overwritten.

### Outputs

**`root.res`: main results**

| Section | Content / what to check |
|---|---|
| PROBLEM DIMENSIONS, CONCENTRATION DATA, VARIANCES | Data as read by MIX: check numbers, names and order. |
| CONSTRAINS | Constraints, if any. |
| CONSTRAINED LEAST SQUARES MIXING RATIOS (INITIAL) | Mixing ratios with the measured end-members fixed (classical EMMA). |
| Attempts → objective function | One value per attempt. The optimum should be reached by several attempts; otherwise increase `--attempts`. |
| OPTIMAL SOLUTION → Expected value of sources concentrations | End-member compositions estimated by MIX; compare with the measured ones. |
| OPTIMAL SOLUTION → Mixing ratios | Final mixing ratios, one row per mixture. |
| Contribution by species → NormlzdDiff | Normalised differences between estimated and measured concentrations; large values flag poorly fitted species. |

**`rootEMMA.out`: EMMA (principal component analysis of the mixtures)**

- Eigenvalues and percentage of explained variance: the number of end-members is usually the number of retained eigenvectors + 1.
- Projections of mixtures and end-members onto the eigenvectors: plotted (e.g. EV1 vs EV2), the end-members should enclose the mixtures.

**`rooteigen.prn`: eigen-analysis**

Eigenvalues and eigenvectors computed by MIX for the source data and for each mixture, with a performance index of the numerical routine. A warning on the raw source data (e.g. *WARNING ERROR 1 from DEPISF*) may appear; it does not stop the computation.

---

## Step 3. Uncertainty of the mixing ratios

```
python mix_uncertainty.py data.res --out sigma.csv
```

| Option | Meaning |
|---|---|
| `--k K` | interval half-width in standard deviations (default 2: ±2σ) |
| `--corr` | print the correlation matrix of the mixing ratios of each mixture |
| `--fixed-em` | use the measured end-members and the initial (least squares) ratios (classical EMMA) |
| `--out FILE` | save the results in a CSV file |

### Reading the output

```
Mixture          EM1             EM2             EM3             EM4    cond
Sample1   0.000 ± 0.015   0.589 ± 0.104   0.000 ± 0.042   0.411 ± 0.079   10.5
```

- **0.589 ± 0.104**: mixing ratio and its ±kσ interval.
- **cond**: condition number. Values around 10 are fine; hundreds or more indicate nearly collinear end-members and unreliable ratios.
- **f_s** (header line): factor per species that rescales the input variances.
- **corr** (with `--corr`): strongly negative values (e.g. −0.8) between two end-members mean that their proportions can be traded off without changing the fit.

The CSV contains, for every mixture and end-member: `lambda`, `sigma` (1σ), the lower and upper limits of the ±kσ interval (clipped to [0, 1]) and the condition number.

### Method and limitations

For each mixture, C_λ = (Eᵀ C_p⁻¹ E)⁻¹, where C_p is the variance of the mixture plus the variances of the end-members weighted by the squared mixing ratios, and input variances are rescaled per species (Carrera and Neuman, 1986). It is a first-order approximation evaluated at the MIX solution:

- sum-to-one and non-negativity are not imposed in the covariance (intervals are conservative and clipped to [0, 1]);
- the uncertainty of the end-member compositions estimated by MIX and the conceptual uncertainty are not included;
- covariance matrices are diagonal.

---

## Acknowledgements

`make_prn.py` is a Python translation of the R script written by **Sebastián Dietrich** (Instituto de Hidrología de Llanuras "Dr. Eduardo J. Usunoff", IHLLA–CONICET) to build MIX input files. We are grateful for his original implementation, on which this workflow is based.

---

## References

- Carrera J., Vázquez-Suñé E., Castillo O., Sánchez-Vila X. (2004). A methodology to compute mixing ratios with uncertain end-members. *Water Resources Research*, 40, W12101. https://doi.org/10.1029/2003WR002263
- Carrera J., Neuman S.P. (1986). Estimation of aquifer parameters under transient and steady state conditions: 1. Maximum likelihood method incorporating prior information. *Water Resources Research*, 22(2), 199–210.
