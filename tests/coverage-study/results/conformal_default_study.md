# Conformal default study
reps=150, alpha=0.10, target coverage=0.90, N=25 T=20 T0=15 r=2, true effect 0.

## Part 1: scale sweep (single treated, block)

| dgp | scale | cov | med.width | bad | joint.point | joint.sim |
|---|---|---|---|---|---|---|
| iid | none | 0.940 | 2.15 | 0 | 0.59 | 0.84 |
| iid | sd | 0.927 | 2.41 | 0 | 0.71 | 0.83 |
| iid | rmspe | 0.927 | 2.41 | 0 | 0.71 | 0.83 |
| iid | mad | 0.913 | 2.63 | 0 | 0.69 | 0.83 |
| iid | diff | 0.927 | 2.51 | 0 | 0.67 | 0.83 |
| ar1 | none | 0.913 | 3.07 | 0 | 0.61 | 0.83 |
| ar1 | sd | 0.927 | 3.34 | 0 | 0.67 | 0.87 |
| ar1 | rmspe | 0.927 | 3.34 | 0 | 0.67 | 0.87 |
| ar1 | mad | 0.920 | 3.66 | 0 | 0.73 | 0.89 |
| ar1 | diff | 0.893 | 3.27 | 0 | 0.69 | 0.87 |
| hetero | none | 1.000 | 5.43 | 0 | 0.97 | 0.99 |
| hetero | sd | 0.920 | 2.82 | 0 | 0.68 | 0.86 |
| hetero | rmspe | 0.920 | 2.82 | 0 | 0.68 | 0.86 |
| hetero | mad | 0.920 | 3.16 | 0 | 0.69 | 0.87 |
| hetero | diff | 0.907 | 3.03 | 0 | 0.66 | 0.85 |
| nonstat | none | 0.920 | 4.51 | 0 | 0.63 | 0.85 |
| nonstat | sd | 0.933 | 4.05 | 0 | 0.67 | 0.85 |
| nonstat | rmspe | 0.933 | 4.05 | 0 | 0.67 | 0.85 |
| nonstat | mad | 0.900 | 4.57 | 0 | 0.66 | 0.80 |
| nonstat | diff | 0.940 | 4.41 | 0 | 0.65 | 0.85 |

## Scale ranking (worst-case scalar coverage, then width)

| scale | worst.cov | mean.cov | mean.width | max.bad |
|---|---|---|---|---|
| sd | 0.920 | 0.927 | 3.16 | 0 |
| rmspe | 0.920 | 0.927 | 3.16 | 0 |
| none | 0.913 | 0.943 | 3.79 | 0 |
| mad | 0.900 | 0.913 | 3.50 | 0 |
| diff | 0.893 | 0.917 | 3.30 | 0 |

**Provisional default scale: sd** (most robust at smallest width).

## Part 2: weight sweep (staggered + heteroskedastic treated)

| weight | cov | med.width | bad |
|---|---|---|---|
| cell | 0.873 | 1.99 | 0 |
| unit | 0.873 | 1.99 | 0 |
| precision | 0.993 | 1.99 | 0 |
