# Layer 5: Computed Summary Tables from CSVs

All numbers below are computed from the CSV outputs, not hand-typed.

## Temperature Sweep (vs 37°C Turner2004 baseline)

### Pooled Retention (total retained / total reference stems)

| tier   |   37C_baseline |    25.0C |    30.0C |    42.0C |
|:-------|---------------:|---------:|---------:|---------:|
| FIRM   |       1        | 0.949985 | 0.949985 | 1        |
| FLOPPY |       0.363618 | 0.272718 | 0.363655 | 0.363618 |
| SOFT   |       0.73334  | 0.600007 | 0.600007 | 0.73334  |

### Per-Case Mean Retention (average across 5 cases)

| tier   |   37C_baseline |   25.0C |   30.0C |   42.0C |
|:-------|---------------:|--------:|--------:|--------:|
| FIRM   |        1       | 0.97142 | 0.97142 | 1       |
| FLOPPY |        0.2857  | 0.26785 | 0.3393  | 0.2857  |
| SOFT   |        0.63334 | 0.53334 | 0.53334 | 0.63334 |

## Parameter Set Sweep (vs Turner2004 baseline)

### Pooled Retention

| tier   |   Turner2004_baseline |   Andronescu2007 |   Langdon2018 |
|:-------|----------------------:|-----------------:|--------------:|
| FIRM   |              1        |        0.79999   |     0.600005  |
| FLOPPY |              0.363618 |        0.454527  |     0.0909091 |
| SOFT   |              0.73334  |        0.0666667 |     0.13334   |

### Per-Case Mean Retention

| tier   |   Turner2004_baseline |   Andronescu2007 |   Langdon2018 |
|:-------|----------------------:|-----------------:|--------------:|
| FIRM   |               1       |          0.76476 |       0.64286 |
| FLOPPY |               0.2857  |          0.4107  |       0.125   |
| SOFT   |               0.63334 |          0.05    |       0.13334 |

## Window Context (vs 0-nt baseline)

### Pooled Retention

| tier   |        0 |       25 |       50 |      100 |
|:-------|---------:|---------:|---------:|---------:|
| FIRM   | 1        | 0.449995 | 0.7      | 0.7      |
| FLOPPY | 0.363618 | 0.545436 | 0.363655 | 0.454527 |
| SOFT   | 0.73334  | 0.2      | 0.2      | 0.133333 |

### Per-Case Mean Retention

| tier   |       0 |      25 |     50 |    100 |
|:-------|--------:|--------:|-------:|-------:|
| FIRM   | 1       | 0.53714 | 0.68   | 0.68   |
| FLOPPY | 0.2857  | 0.5357  | 0.3393 | 0.4107 |
| SOFT   | 0.63334 | 0.15    | 0.15   | 0.1    |

## Window Context: Mean MEA BP Distance

| Flank Size (nt) | Mean BP Distance |
|-----------------|------------------|
| 0 | 0.0 |
| 25 | 30.4 |
| 50 | 19.6 |
| 100 | 22.4 |
