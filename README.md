# Robot Polishing Tool-Orientation Optimization

This repository contains the MATLAB implementation associated with a study on
tool-orientation planning for robotic polishing of concave cavity surfaces under
point-wise heterogeneous feasible domains.

The repository focuses on two core computational stages:

1. construction of the feasible orientation domain at each cutter-contact point
   by combining spherical-contact and collision constraints; and
2. multi-objective optimization of the tool-orientation sequence using a
   non-dominated sorting whale optimization algorithm (NSWOA).

## Repository structure

```text
Robot-Polishing-Tool-Orientation-Optimization/
├── README.md
├── LICENSE
├── CITATION.cff
├── .gitignore
├── code/
│   ├── 01_compute_feasible_orientation_domain.m
│   └── 02_run_NSWOA_tool_orientation.m
├── data/
│   ├── tool_points.csv
│   └── feasible_orientation_bounds.csv
├── results/
│   └── README.md
└── docs/
    └── algorithm_description.md
```

## Requirements

- MATLAB
- No additional third-party toolbox is required by the two core scripts.

The code was prepared as a standalone reproducibility package for the
optimization procedure reported in the associated paper.

## Workflow

### Step 1 — Compute point-wise feasible orientation domains

Run:

```matlab
01_compute_feasible_orientation_domain
```

This script:

- defines the polishing tool/workpiece geometric parameters;
- loads the 31 cutter-contact points;
- calculates the spherical-contact constraints;
- calculates the collision constraints;
- intersects the two constraint intervals at each point; and
- outputs the final feasible bounds for `alpha_i` and `beta_i`.

The script writes:

```text
tool_orientation_range_rad.xlsx
tool_orientation_range_deg.xlsx
theta_result.xlsx
```

These files are generated locally and are excluded from version control by
default because they can be reproduced from the source code.

### Step 2 — Run NSWOA multi-objective optimization

Run:

```matlab
02_run_NSWOA_tool_orientation
```

The decision vector contains 62 variables:

```text
[alpha_1, ..., alpha_31, beta_1, ..., beta_31]
```

with point-dependent lower and upper bounds.

The implementation combines the basic whale optimization algorithm with:

- non-dominated sorting;
- crowding-distance based environmental selection; and
- Pareto-front based leader selection.

The two objectives are:

1. maximize the posture-dependent material-removal capability indicator
   (implemented as minimization of the negative mean indicator); and
2. minimize the mean squared spatial angle between adjacent tool-axis vectors.

The script generates the Pareto front, convergence figures, the selected
compromise solution, the optimized `alpha`/`beta` sequences, and adjacent
tool-axis angles.

## Data files

`data/tool_points.csv` contains the 31 cutter-contact points and the associated
`theta_i` values.

`data/feasible_orientation_bounds.csv` contains the point-wise lower and upper
bounds for `alpha_i` and `beta_i` used by the NSWOA script.

These CSV files are provided to make the numerical inputs easy to inspect
without reading the MATLAB source.

## Reproducibility note

The NSWOA implementation contains stochastic operations. Therefore, individual
Pareto solutions may vary between runs unless a random seed is fixed.

For exact run-to-run reproducibility, add a MATLAB random-number seed before
population initialization, for example:

```matlab
rng(1);
```

The original source code is preserved in this repository without adding a fixed
seed, so that the released implementation remains consistent with the supplied
research code.

## Main parameters

The current source code uses:

```text
Tool radius RT = 4
Compression parameter c0 = 1
Number of cutter-contact points = 31
NSWOA population size = 100
Maximum iterations = 200
Number of objectives = 2
```

Additional geometric parameters are defined near the beginning of
`01_compute_feasible_orientation_domain.m`.

## Associated paper

**Paper title:** Material-Removal-Driven Multi-Objective Tool Orientation Planning for Robotic Polishing under Point-Wise Heterogeneous Feasible Orientation Domains

**Authors:** Kang Fu, Wenjun Yu, Xiaogao Li, Li Zhang, Fanqiang Bu, Xingguo Wang

**Journal:** Advanced Engineering Informatics

**DOI:** 

If you use this code, please cite the associated paper.

## License

This project is released under the MIT License. See `LICENSE`.

## Contact

For questions about the implementation or reproducibility, please contact the
corresponding author listed in the associated paper.
