# EMMA-MIX workflow: from a CSV table to EMMA-MIX mixing ratios and their uncertainty

<!-- DOI badge: replace after the Zenodo release, e.g.
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX) -->

```
data.csv ──► make_prn.py ──► root.prn ──► EMMA-MIX (console) ──► root.res ──► mix_uncertainty.py ──► mixing ratios ± σ
```

| Step | Tool | Input | Output |
|---|---|---|---|
| 1 | `make_prn.py` | `data.csv` | `root.prn` |
| 2 | EMMA-MIX (MIX program, Carrera et al., 2004) | `root.prn` | `root.res`, `rootEMMA.out`, `rooteigen.prn` |
| 3 | `mix_uncertainty.py` | `root.res` | uncertainty of the mixing ratios (screen and `.csv`) |

## Contents

| File / folder | Description |
|---|---|
| `make_prn.py` | Builds the MIX input file (`.prn`) from a CSV table |
| `mix_uncertainty.py` | Uncertainty of the mixing ratios from the MIX output (`.res`) |
| `data.csv` | Example input file |
| `MIX/` | EMMA-MIX (MIX) executable (`mix_06_lowdim.exe`) and *MIX Guide Notes* (`Guide_Notes_MIX_2-2.pdf`; Vázquez-Suñé and Carrera) |
| `CITATION.cff` | How to cite this software |
| `LICENSE` | MIT license (Python scripts and documentation; see *License*) |

## Requirements

- Python ≥ 3.8 with `numpy` and `pandas`:
  ```
  pip install numpy pandas
  ```
- The EMMA-MIX (MIX) executable for Windows, `mix_06_lowdim.exe`, provided in the `MIX/` folder. The format of its input and output files is described in the *MIX Guide Notes* in the same folder.

Commands are typed in a terminal opened in the folder that contains the scripts and the data: on Windows, Command Prompt (`cmd`), PowerShell, or **Terminal → New Terminal** in VS Code. Type one command at a time and press **Enter**.

---

## Step 1. Build the `.prn` file

### Input file: `data.csv`

One row per water, one column per species, and the **name of each water in the last column**. **End-members go in the first rows, followed by the mixtures.** MIX keeps this order, so the rows of `data.csv`, the `.prn` and the `.res` match. In the example file (values shortened here), the first four rows are the end-members and the rest are mixtures:

```
Ca;Mg;Na;K;Cl;NO3;SO4;TIC;F;CE;2H;O18;id
1.7186;1.6901;2.5490;0.3754;0.4545;1.4998;0.1230;7.5622;0.0542;878;-30.3;-5.47;P29
0.1372;0.0247;0.2653;0.1176;0.0282;0.0290;0.0343;0.4916;0.0053;67;-36.7;-6.24;Rainfall
0.4588;0.5303;11.4472;0.3022;0.4083;0.1645;0.1599;12.0851;0.1995;1136;-27.8;-5.08;P38
1.9626;1.2619;3.0493;0.0908;0.2201;0.0461;0.0661;7.2804;0.0256;698;-23.8;-4.44;UZ
0.5769;0.5086;1.0492;0.4230;0.1412;0.0371;0.0135;3.3362;0.0230;337;-35.1;-5.88;S18122023
...
```

- **Columns separated by** `;` and **decimal point** `.` (each row on its own line). Files with other separators can be read with `--sep ","` and/or `--decimal ","`.
- **Units** are not converted; use a consistent set (e.g. ions in mmol/L, electrical conductivity in µS/cm, δ²H and δ¹⁸O in ‰ with their sign).
- **Names**: no spaces, unique, preferably ≤ 10 characters.
- **No empty cells.**

### Running `make_prn.py`

There is a single script, `make_prn.py`, which can be used in two ways.

**a) Interactive: only the CSV file is given, and the script asks**

```
python make_prn.py data.csv
```

The script lists the waters and asks:

1. Which rows are end-members (e.g. `1-4`).
2. Which rows to leave out of the run (Enter = none).
3. Whether the selection is correct (`y`/`n`).
4. Optional constraints: fix the proportion of a source (end-member) in a well (mixture), written as `well,source,proportion` (e.g. `5,2,0`). Enter = no constraints. Constraints are explained in the *MIX Guide Notes* (`MIX/Guide_Notes_MIX_2-2.pdf`, p. 3; examples on pp. 7, 9 and 10).

Constraints use the MIX terms and numbering (*MIX Guide Notes*, p. 3): mixtures are *wells* and end-members are *sources*, numbered separately. The well number counts only the mixtures (1 = first mixture) and the source number only the end-members (1 = first end-member); neither is the row of the CSV. The script shows both numbered lists before asking. In the example file, `5,2,0` sets the proportion of `Rainfall` (source 2) to 0 in `S13012022` (well 5, row 9 of the CSV).

**b) Non-interactive: the answers are given as options, and the script asks nothing**

```
python make_prn.py data.csv --sources 4
python make_prn.py data.csv --em 1-3 --exclude 4 --out run3EM.prn
python make_prn.py data.csv --sources 4 --constraint 5,2,0 --out run4EM_c.prn
```

| Option | Meaning |
|---|---|
| `--sources N` | the first N rows are end-members |
| `--em LIST` | rows that are end-members (e.g. `1-3`); with `--exclude`, to run with a subset of the end-members |
| `--exclude LIST` | rows left out of the run |
| `--constraint w,s,p` | fix proportion `p` of source `s` in well `w`, numbered as in MIX (see above; repeatable) |
| `--var-csv FILE` | user-defined variances for some waters (same columns as `data.csv`, rows matched by name) |
| `--em-var-factor F` | multiply the variances of all end-members by `F` |
| `--attempts N` | number of MIX attempts (default 9) |
| `--rng-seed N` | fixed seed for the initial mixing ratios (identical `.prn` every time) |
| `--out FILE` | name of the `.prn` (default: same name as the CSV, e.g. `data.prn`) |

`python make_prn.py --help` lists all the options.

### What the `.prn` contains

The layout required by MIX (MIX Guide Notes): dimensions; concentrations (end-members first, then mixtures); variances; constraints; initial mixing ratios.

- **Variances:** by default, the variance of each concentration equals its squared value (σ = c), for end-members and mixtures. See *Choosing the variances*.
- **Constraints:** first the number of constraints, then one line per constraint with the well, the source and the fixed proportion. With no constraints, the block is just `0`. For example, two constraints (rainfall = 0 in wells 4 and 5), written with `--constraint 4,2,0 --constraint 5,2,0`:
  ```
               2
               4              2              0
               5              2              0
  ```
- **Initial mixing ratios:** random values, each row summing to 1. They do not change the final solution, because MIX tries several starting points: the input ratios, the least squares ratios, equal proportions and, from the fourth attempt on, the input ratios with random noise (*MIX Guide Notes*, p. 4).

### Choosing the variances

The variances act as weights: they set how much each species counts in the fit and how much MIX can move each end-member. Multiplying all of them by the same factor does not change the results, but their relative values do.

By default, `make_prn.py` sets the variance of each concentration equal to its squared value (σ = c), for end-members and mixtures alike. Variances can also be built according to the conceptual model of each study. For example, mixtures can be given the analytical error, and end-members the analytical error plus the natural variability of the source they represent (σ² = σ²_analytical + σ²_natural). For stable isotopes, σ can be given in ‰ independently of the δ value.

User-defined variances are loaded with `--var-csv FILE`: a CSV with the same columns as `data.csv` (species, and the water name in the last column), containing variances instead of concentrations. Only the waters listed in the file are changed; the others keep the default.

---

## Step 2. Run EMMA-MIX

1. Put the **MIX executable (`mix_06_lowdim.exe`, from the `MIX/` folder) and the `.prn` file in the same folder**.
2. Double-click on `mix_06_lowdim.exe`, or run it from the terminal: type `mix_06_lowdim.exe` in Command Prompt, or `.\mix_06_lowdim.exe` in PowerShell and in the VS Code terminal.
3. When it asks for the **root name**, type the name of the `.prn` file **without the extension** (e.g. `data`) and press Enter.
4. The output files are written in the same folder. Previous outputs with the same root name are overwritten (*MIX Guide Notes*, pp. 11–12).

### Outputs

MIX writes three files in the same folder:

| File | Content |
|---|---|
| `root.res` | Main results: input data as read by MIX, initial (least squares) mixing ratios, each attempt, and the optimal solution (estimated end-member compositions, mixing ratios and contribution of each species to the objective function). Described in the *MIX Guide Notes* (pp. 12–19). |
| `rooteigen.prn` | Eigen-analysis of the source data and of each mixture (*MIX Guide Notes*, p. 20, where it is named `rooteigen.res`). A warning such as *WARNING ERROR 1 from DEPISF* may appear; it does not stop the computation. |
| `rootEMMA.out` | EMMA (principal component analysis of the mixtures): eigenvalues, explained variance and projections of mixtures and end-members onto the eigenvectors, listed without names in the order of the `.prn` (mixtures first, then end-members). |

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
Mixture                       P29             Rainfall                  P38                   UZ    cond
--------------------------------------------------------------------------------------------------------
S18122023           0.000 ± 0.015        0.589 ± 0.104        0.000 ± 0.042        0.411 ± 0.079    10.5
```

- **0.589 ± 0.104**: mixing ratio and its ±kσ interval.
- **cond**: condition number. Values around 10 are fine; hundreds or more indicate nearly collinear end-members and unreliable ratios.
- **f_s** (line above the table): factor per species that rescales the input variances.
- **corr** (with `--corr`): strongly negative values (e.g. −0.8) between two end-members mean that their proportions can be traded off without changing the fit.

The CSV contains, for every mixture and end-member: `lambda`, `sigma` (1σ), the lower and upper limits of the ±kσ interval (clipped to [0, 1]) and the condition number.

### Method and limitations

The uncertainty formulation is due to J. Carrera. For each mixture, C_λ = (Eᵀ C_p⁻¹ E)⁻¹, where C_p is the variance of the mixture plus the variances of the end-members weighted by the squared mixing ratios. Input variances are rescaled by a factor per species (f_s) so that the mean squared normalised residual of each species equals 1. This is a first-order approximation, evaluated at the MIX solution, with the following limitations:

- **Sum to one and non-negativity are not imposed.** The mixing ratios of a mixture are treated as independent, so the intervals are wider than needed (conservative). They are clipped to [0, 1].
- **Errors are assumed independent between species** (diagonal covariance matrices, as in MIX).

---

## Applications of EMMA-MIX

EMMA-MIX has been used, among others, to:

- identify the sources of urban groundwater recharge in Barcelona (Vázquez-Suñé et al., 2010);
- compute mixing ratios and identify hydrochemical reactions in an aquifer recharged by the Besòs River (Tubau et al., 2014);
- characterise mountain system recharge in the Sierra Nevada, California (Armengol et al., 2024).

Pelizardi et al. (2017) extended the method to reactive mixing by using conservative chemical components instead of individual species. Goyetche et al. (2022) applied this reactive approach (rEMMA) to quantify mixing ratios and the extent of chemical reactions in a coastal aquifer.

---

## Authors

| Author | Contribution |
|---|---|
| **Jesús Carrera** | MIX method and formulation of the uncertainty of the mixing ratios |
| **Romina Marcovecchio** | Python code (`make_prn.py`, `mix_uncertainty.py`), workflow and documentation |
| **Sebastián Dietrich** | Original R script to build the MIX input files, translated and extended in `make_prn.py` |

The MIX program and its *Guide Notes* were developed by the Hydrogeology Group (GHS, UPC–CSIC; Carrera et al., 2004).

## How to cite

See `CITATION.cff`, or use **"Cite this repository"** on the right side of the GitHub page.

## License

The Python scripts and the documentation are released under the MIT license (`LICENSE`). The MIX executable and the *MIX Guide Notes* in `MIX/` are not covered by this license: they belong to the Hydrogeology Group (GHS, UPC–CSIC) and are included with the permission of J. Carrera. Please cite Carrera et al. (2004) when using MIX.

---

## References

- Armengol S., Ajami H., Acero Triana J.S., Sickman J.O., Ortega L. (2024). Isogeochemical characterization of mountain system recharge processes in the Sierra Nevada, California. *Water Resources Research*, 60, e2023WR035719. https://doi.org/10.1029/2023WR035719
- Carrera J., Vázquez-Suñé E., Castillo O., Sánchez-Vila X. (2004). A methodology to compute mixing ratios with uncertain end-members. *Water Resources Research*, 40, W12101. https://doi.org/10.1029/2003WR002263
- Goyetche T., Luquot L., Carrera J., Martínez-Pérez L., Folch A. (2022). Identification and quantification of chemical reactions in a coastal aquifer to assess submarine groundwater discharge composition. *Science of the Total Environment*, 838, 155978. https://doi.org/10.1016/j.scitotenv.2022.155978
- Pelizardi F., Bea S.A., Carrera J., Vives L. (2017). Identifying geochemical processes using End Member Mixing Analysis to decouple chemical components for mixing ratio calculations. *Journal of Hydrology*, 550, 144–156. https://doi.org/10.1016/j.jhydrol.2017.04.010
- Tubau I., Vázquez-Suñé E., Jurado A., Carrera J. (2014). Using EMMA and MIX analysis to assess mixing ratios and to identify hydrochemical reactions in groundwater. *Science of the Total Environment*, 470–471, 1120–1131. https://doi.org/10.1016/j.scitotenv.2013.10.121
- Vázquez-Suñé E., Carrera J. *MIX Guide Notes* (v2.2). Included in `MIX/`.
- Vázquez-Suñé E., Carrera J., Tubau I., Sánchez-Vila X., Soler A. (2010). An approach to identify urban groundwater recharge. *Hydrology and Earth System Sciences*, 14, 2085–2097. https://doi.org/10.5194/hess-14-2085-2010
