# Algorithm description

## 1. Point-wise feasible orientation domain

For each cutter-contact point `P_i`, the tool orientation is represented by two
angles, `alpha_i` and `beta_i`.

The first MATLAB script evaluates two classes of constraints:

- spherical-contact constraints; and
- collision constraints.

The final feasible interval at each point is obtained by intersecting the
corresponding intervals:

```text
alpha_i,min = max(alpha_i,min^sph, alpha_i,min^col)
alpha_i,max = min(alpha_i,max^sph, alpha_i,max^col)

beta_i,min  = max(beta_i,min^sph,  beta_i,min^col)
beta_i,max  = min(beta_i,max^sph,  beta_i,max^col)
```

This produces a point-wise heterogeneous feasible domain along the polishing
trajectory.

## 2. Multi-objective orientation-sequence optimization

The optimization variable is the complete orientation sequence:

```text
x = [alpha_1, ..., alpha_N, beta_1, ..., beta_N]
```

where `N = 31`.

Each variable is constrained by its corresponding point-wise feasible bounds.

## 3. Objective functions

### Objective 1: material-removal capability

The code evaluates a posture-dependent material-removal capability indicator
`I` and minimizes its negative mean value:

```text
F1 = -mean(I)
```

Thus, minimizing `F1` is equivalent to maximizing the average removal
capability.

### Objective 2: orientation continuity

The spatial angle `gamma_i` between adjacent tool-axis vectors is calculated
from the dot product of adjacent normalized tool-axis vectors.

The second objective is:

```text
F2 = mean(gamma_i^2)
```

which penalizes abrupt changes between neighboring orientations.

## 4. NSWOA procedure

The supplied implementation uses a basic whale optimization algorithm combined
with non-dominated sorting.

At each iteration:

1. leaders are selected from the current first Pareto front;
2. standard WOA encircling, exploration, and spiral-update mechanisms generate
   offspring;
3. decision variables are clipped to their point-dependent bounds;
4. parent and offspring populations are merged;
5. non-dominated sorting is applied; and
6. crowding distance is used when only part of a Pareto front can be retained.

No mutation operator or additional advanced WOA variant is introduced.

## 5. Compromise solution

After the final Pareto front is obtained, the code normalizes the objective
values and selects the Pareto solution with the minimum Euclidean distance to
the ideal point as a representative compromise solution.
