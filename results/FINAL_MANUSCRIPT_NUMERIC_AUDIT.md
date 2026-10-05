# Final Manuscript Numeric Audit

All checks below are computed from generated result files.

- Symbolic candidates: 599
- Rank-1 candidates: 86
- Rank-1 cells: 24/48
- Spatial/cruise validation rows: 40

## Rank-1 temporal and external RMSE
| basin          |   validation_RMSE |   external_RMSE |   parameter_count |   expression_complexity |
|:---------------|------------------:|----------------:|------------------:|------------------------:|
| Atlantic       |           20.1715 |        17.672   |                 5 |                       8 |
| Indian         |           14.4825 |        18.4129  |                 5 |                       8 |
| Pacific        |           16.0016 |        13.8208  |                 5 |                       8 |
| Southern Ocean |           10.5753 |         8.98541 |                 5 |                       8 |

## Rank-2 deltas
| basin          |   validation_delta |   external_delta |
|:---------------|-------------------:|-----------------:|
| Atlantic       |          1.98358   |        1.91441   |
| Indian         |          0.254342  |        0.444669  |
| Pacific        |         -0.289342  |        0.266558  |
| Southern Ocean |         -0.0348738 |       -0.0905896 |

## Blocked validation
| basin          | block_type   | model           |   n_test |   n_groups_total |   RMSE_mean |   RMSE_sd |   R2_mean |    Bias_mean |
|:---------------|:-------------|:----------------|---------:|-----------------:|------------:|----------:|----------:|-------------:|
| Atlantic       | cruise       | frozen_symbolic |   137198 |              192 |    17.5187  |  4.68891  |  0.892813 | -1.72421     |
| Atlantic       | cruise       | full_quadratic  |   137198 |              192 |    15.6289  |  2.77652  |  0.914445 |  0.238396    |
| Atlantic       | cruise       | linear          |   137198 |              192 |    16.6669  |  4.33576  |  0.90325  |  0.067073    |
| Atlantic       | cruise       | rank1           |   137198 |              192 |    16.9096  |  3.88643  |  0.900469 |  0.383379    |
| Atlantic       | cruise       | rank2           |   137198 |              192 |    15.5794  |  2.50227  |  0.91468  |  0.181039    |
| Atlantic       | spatial      | frozen_symbolic |   137198 |              159 |    17.927   |  1.97726  |  0.887704 | -1.72422     |
| Atlantic       | spatial      | full_quadratic  |   137198 |              159 |    15.5468  |  2.32314  |  0.915392 |  0.129989    |
| Atlantic       | spatial      | linear          |   137198 |              159 |    16.9258  |  2.62604  |  0.899591 |  0.0648965   |
| Atlantic       | spatial      | rank1           |   137198 |              159 |    16.2991  |  2.10955  |  0.907336 |  0.0954677   |
| Atlantic       | spatial      | rank2           |   137198 |              159 |    15.4638  |  2.11337  |  0.916529 |  0.0114394   |
| Indian         | cruise       | frozen_symbolic |    36262 |               40 |    14.0636  |  3.69781  |  0.983805 |  0.837146    |
| Indian         | cruise       | full_quadratic  |    36262 |               40 |    15.3747  |  4.68654  |  0.980751 | -0.322644    |
| Indian         | cruise       | linear          |    36262 |               40 |    17.5543  |  4.49232  |  0.975269 | -0.277819    |
| Indian         | cruise       | rank1           |    36262 |               40 |    17.3937  |  5.59291  |  0.97514  | -0.286973    |
| Indian         | cruise       | rank2           |    36262 |               40 |    16.9806  |  4.92726  |  0.976551 | -0.197074    |
| Indian         | spatial      | frozen_symbolic |    36262 |               88 |    14.3539  |  1.83012  |  0.984117 |  0.834998    |
| Indian         | spatial      | full_quadratic  |    36262 |               88 |    14.9191  |  2.55322  |  0.982562 | -0.0311584   |
| Indian         | spatial      | linear          |    36262 |               88 |    16.989   |  1.87858  |  0.97788  | -0.131085    |
| Indian         | spatial      | rank1           |    36262 |               88 |    16.6088  |  2.09327  |  0.97874  | -0.11757     |
| Indian         | spatial      | rank2           |    36262 |               88 |    16.0034  |  1.805    |  0.980355 | -0.000528867 |
| Pacific        | cruise       | frozen_symbolic |   169993 |              383 |    15.5198  |  2.43724  |  0.987005 | -0.368102    |
| Pacific        | cruise       | full_quadratic  |   169993 |              383 |    16.0967  |  3.72954  |  0.985683 | -0.248321    |
| Pacific        | cruise       | linear          |   169993 |              383 |    18.5933  |  1.86382  |  0.981676 | -0.0408604   |
| Pacific        | cruise       | rank1           |   169993 |              383 |    16.0312  |  2.32724  |  0.986191 | -0.0451089   |
| Pacific        | cruise       | rank2           |   169993 |              383 |    15.8656  |  2.9976   |  0.986287 | -0.0956999   |
| Pacific        | spatial      | frozen_symbolic |   169993 |              241 |    15.6343  |  1.21682  |  0.987172 | -0.368103    |
| Pacific        | spatial      | full_quadratic  |   169993 |              241 |    16.4284  |  2.94763  |  0.985676 | -0.0980412   |
| Pacific        | spatial      | linear          |   169993 |              241 |    18.6227  |  1.0481   |  0.981795 |  0.063015    |
| Pacific        | spatial      | rank1           |   169993 |              241 |    16.0307  |  1.27989  |  0.986508 |  0.0167783   |
| Pacific        | spatial      | rank2           |   169993 |              241 |    16.2557  |  2.54682  |  0.986048 | -0.123989    |
| Southern Ocean | cruise       | frozen_symbolic |    88528 |              114 |     9.05183 |  0.70281  |  0.97751  | -0.15094     |
| Southern Ocean | cruise       | full_quadratic  |    88528 |              114 |     9.29527 |  0.933874 |  0.976208 | -0.0147427   |
| Southern Ocean | cruise       | linear          |    88528 |              114 |     9.41838 |  0.748223 |  0.975643 |  0.00687404  |
| Southern Ocean | cruise       | rank1           |    88528 |              114 |     9.22486 |  0.796829 |  0.976614 |  0.0129623   |
| Southern Ocean | cruise       | rank2           |    88528 |              114 |     9.19451 |  0.892687 |  0.97673  |  0.000625239 |
| Southern Ocean | spatial      | frozen_symbolic |    88528 |              232 |     9.04058 |  0.863095 |  0.977156 | -0.150999    |
| Southern Ocean | spatial      | full_quadratic  |    88528 |              232 |     9.38101 |  1.61664  |  0.974746 | -0.0174392   |
| Southern Ocean | spatial      | linear          |    88528 |              232 |     9.34213 |  0.854123 |  0.97568  |  0.0191396   |
| Southern Ocean | spatial      | rank1           |    88528 |              232 |     9.15729 |  0.939105 |  0.97652  |  0.014456    |
| Southern Ocean | spatial      | rank2           |    88528 |              232 |     9.1378  |  1.1468   |  0.976441 |  0.000289698 |

## AAIW and abyssal failure
| Model    |    RMSE |    Bias |        R2 |   n |
|:---------|--------:|--------:|----------:|----:|
| Rank1    | 28.4257 | 27.1427 | -0.114958 | 428 |
| FullQuad | 22.3334 | 21.1264 |  0.311749 | 428 |

Southern abyss global R2=-0.173896, local R2=0.785601, local RMSE=3.182377.

## Consistency checks
- candidate_count_599: PASS
- rank1_count_86: PASS
- rank1_cells_24_of_48: PASS
- blocked_rows_pre2018: PASS
- cruise_and_spatial_present: PASS
- abyss_global_negative: PASS
- abyss_local_positive: PASS
