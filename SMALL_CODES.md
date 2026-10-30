# Small effective codes for every prime p <= 97

*2026-10-30. Companion to `COVERAGE.md` (which ranks by gamma). Here the question is the practical one: at each prime,
which small codes exist, how good must the raw magic states be (threshold eps*), how fast does one round suppress
error (eps_out ~ c eps^d), and how many blocks are accepted (P_acc). All thresholds and coefficients are EXACT:
Reed-Solomon codes from the closed-form MDS engine (`qmsd.rs_family`, `runs/coverage/rs_family_all.json`), the
rest from the enumerating engine (`qmsd/data/performance/`). eps is the convention of arXiv:2610.08012 (probability
of a nontrivial Z error per input).*

## Reading guide

- **k = 1 Reed-Solomon codes [[p-1, 1, r+1]]_p** are the best single-output codes at every prime >= 11: thresholds
  0.5-0.7 and leading coefficients that fall like p^{-(d-1)}; e.g. [[22,1,8]]_23 has eps* = 0.64 and c = 1.3e-4.
- Larger k trades distance for rate along the family (d = r - k + 2); the minimum-gamma member is the last column.
- Below p = 11 the RS family is trivial (d = 2) and the small codes are the cyclic / doubling / paper codes.

## Reed-Solomon family, exact (per-output threshold eps*, coefficient c of eps_out ~ c eps^d, P_acc at eps = 1e-2)

| p | r | k = 1 code | eps* | c | P_acc(1e-2) | min-gamma member | gamma | eps* | all k |
|---|---|---|---|---|---|---|---|---|---|
| 5 | 1 | [[4,1,2]] | 0.363 | 1.50e+00 | 0.96 | [[4,1,2]] | 2.000 | 0.363 | [[4,1,2]] eps*=0.36 |
| 7 | 1 | [[6,1,2]] | 0.232 | 2.50e+00 | 0.94 | [[6,1,2]] | 2.585 | 0.232 | [[6,1,2]] eps*=0.23 |
| 11 | 3 | [[10,1,4]] | 0.532 | 2.10e-01 | 0.90 | [[9,2,3]] | 1.369 | 0.409 | [[10,1,4]] eps*=0.53, [[9,2,3]] eps*=0.41, [[8,3,2]] eps*=0.19 |
| 13 | 3 | [[12,1,4]] | 0.479 | 2.86e-01 | 0.89 | [[11,2,3]] | 1.552 | 0.353 | [[12,1,4]] eps*=0.48, [[11,2,3]] eps*=0.35, [[10,3,2]] eps*=0.15 |
| 17 | 5 | [[16,1,6]] | 0.602 | 7.64e-03 | 0.85 | [[13,4,3]] | 1.073 | 0.327 | [[16,1,6]] eps*=0.60, [[15,2,5]] eps*=0.54, [[14,3,4]] eps*=0.45, [[13,4,3]] eps*=0.33, [[12,5,2]] eps*=0.13 |
| 19 | 5 | [[18,1,6]] | 0.569 | 9.82e-03 | 0.83 | [[15,4,3]] | 1.203 | 0.293 | [[18,1,6]] eps*=0.57, [[17,2,5]] eps*=0.50, [[16,3,4]] eps*=0.42, [[15,4,3]] eps*=0.29, [[14,5,2]] eps*=0.11 |
| 23 | 7 | [[22,1,8]] | 0.643 | 1.28e-04 | 0.80 | [[18,5,4]] | 0.924 | 0.396 | [[22,1,8]] eps*=0.64, [[21,2,7]] eps*=0.60, [[20,3,6]] eps*=0.55, [[19,4,5]] eps*=0.48, [[18,5,4]] eps*=0.40, [[17,6,3]] eps*=0.27, [[16,7,2]] eps*=0.10 |
| 29 | 9 | [[28,1,10]] | 0.671 | 1.24e-06 | 0.75 | [[22,7,4]] | 0.826 | 0.353 | [[28,1,10]] eps*=0.67, [[27,2,9]] eps*=0.64, [[26,3,8]] eps*=0.60, [[25,4,7]] eps*=0.56, [[24,5,6]] eps*=0.51, [[23,6,5]] eps*=0.44, [[22,7,4]] eps*=0.35, [[21,8,3]] eps*=0.24, [[20,9,2]] eps*=0.08 |
| 31 | 9 | [[30,1,10]] | 0.652 | 1.53e-06 | 0.74 | [[25,6,5]] | 0.887 | 0.416 | [[30,1,10]] eps*=0.65, [[29,2,9]] eps*=0.62, [[28,3,8]] eps*=0.58, [[27,4,7]] eps*=0.54, [[26,5,6]] eps*=0.48, [[25,6,5]] eps*=0.42, [[24,7,4]] eps*=0.33, [[23,8,3]] eps*=0.22, [[22,9,2]] eps*=0.07 |
| 37 | 11 | [[36,1,12]] | 0.676 | 9.51e-09 | 0.70 | [[29,8,5]] | 0.800 | 0.383 | [[36,1,12]] eps*=0.68, [[35,2,11]] eps*=0.65, [[34,3,10]] eps*=0.62, [[33,4,9]] eps*=0.59, [[32,5,8]] eps*=0.55, [[31,6,7]] eps*=0.50, [[30,7,6]] eps*=0.45, [[29,8,5]] eps*=0.38, [[28,9,4]] eps*=0.30, [[27,10,3]] eps*=0.20, [[26,11,2]] eps*=0.06 |
| 41 | 13 | [[40,1,14]] | 0.708 | 3.46e-11 | 0.67 | [[31,10,5]] | 0.703 | 0.371 | [[40,1,14]] eps*=0.71, [[39,2,13]] eps*=0.69, [[38,3,12]] eps*=0.67, [[37,4,11]] eps*=0.64, [[36,5,10]] eps*=0.61, [[35,6,9]] eps*=0.58, [[34,7,8]] eps*=0.54, [[33,8,7]] eps*=0.49, [[32,9,6]] eps*=0.44, [[31,10,5]] eps*=0.37, [[30,11,4]] eps*=0.29, [[29,12,3]] eps*=0.19, [[28,13,2]] eps*=0.06 |
| 43 | 13 | [[42,1,14]] | 0.694 | 4.18e-11 | 0.66 | [[34,9,6]] | 0.742 | 0.420 | [[42,1,14]] eps*=0.69, [[41,2,13]] eps*=0.67, [[40,3,12]] eps*=0.65, [[39,4,11]] eps*=0.62, [[38,5,10]] eps*=0.59, [[37,6,9]] eps*=0.56, [[36,7,8]] eps*=0.52, [[35,8,7]] eps*=0.47, [[34,9,6]] eps*=0.42, [[33,10,5]] eps*=0.36, [[32,11,4]] eps*=0.28, [[31,12,3]] eps*=0.18, [[30,13,2]] eps*=0.05 |
| 47 | 15 | [[46,1,16]] | 0.721 | 1.13e-13 | 0.63 | [[36,11,6]] | 0.662 | 0.409 | [[46,1,16]] eps*=0.72, [[45,2,15]] eps*=0.70, [[44,3,14]] eps*=0.69, [[43,4,13]] eps*=0.66, [[42,5,12]] eps*=0.64, [[41,6,11]] eps*=0.61, [[40,7,10]] eps*=0.58, [[39,8,9]] eps*=0.55, [[38,9,8]] eps*=0.51, [[37,10,7]] eps*=0.46, [[36,11,6]] eps*=0.41, [[35,12,5]] eps*=0.34, [[34,13,4]] eps*=0.27, [[33,14,3]] eps*=0.17, [[32,15,2]] eps*=0.05 |
| 53 | 17 | [[52,1,18]] | 0.733 | 2.87e-16 | 0.59 | [[40,13,6]] | 0.627 | 0.384 | [[52,1,18]] eps*=0.73, [[51,2,17]] eps*=0.72, [[50,3,16]] eps*=0.70, [[49,4,15]] eps*=0.68, [[48,5,14]] eps*=0.66, [[47,6,13]] eps*=0.64, [[46,7,12]] eps*=0.62, [[45,8,11]] eps*=0.59, [[44,9,10]] eps*=0.56, [[43,10,9]] eps*=0.52, [[42,11,8]] eps*=0.48, [[41,12,7]] eps*=0.44, [[40,13,6]] eps*=0.38, [[39,14,5]] eps*=0.32, [[38,15,4]] eps*=0.25, [[37,16,3]] eps*=0.16, [[36,17,2]] eps*=0.04 |
| 59 | 19 | [[58,1,20]] | 0.742 | 5.77e-19 | 0.56 | [[45,14,7]] | 0.600 | 0.415 | [[58,1,20]] eps*=0.74, [[57,2,19]] eps*=0.73, [[56,3,18]] eps*=0.72, [[55,4,17]] eps*=0.70, [[54,5,16]] eps*=0.68, [[53,6,15]] eps*=0.66, [[52,7,14]] eps*=0.64, [[51,8,13]] eps*=0.62, [[50,9,12]] eps*=0.60, [[49,10,11]] eps*=0.57, [[48,11,10]] eps*=0.54, [[47,12,9]] eps*=0.50, [[46,13,8]] eps*=0.46, [[45,14,7]] eps*=0.42, [[44,15,6]] eps*=0.36, [[43,16,5]] eps*=0.30, [[42,17,4]] eps*=0.23, [[41,18,3]] eps*=0.15, [[40,19,2]] eps*=0.04 |
| 61 | 19 | [[60,1,20]] | 0.733 | 6.88e-19 | 0.55 | [[47,14,7]] | 0.622 | 0.404 | [[60,1,20]] eps*=0.73, [[59,2,19]] eps*=0.72, [[58,3,18]] eps*=0.71, [[57,4,17]] eps*=0.69, [[56,5,16]] eps*=0.67, [[55,6,15]] eps*=0.65, [[54,7,14]] eps*=0.63, [[53,8,13]] eps*=0.61, [[52,9,12]] eps*=0.58, [[51,10,11]] eps*=0.56, [[50,11,10]] eps*=0.52, [[49,12,9]] eps*=0.49, [[48,13,8]] eps*=0.45, [[47,14,7]] eps*=0.40, [[46,15,6]] eps*=0.35, [[45,16,5]] eps*=0.29, [[44,17,4]] eps*=0.22, [[43,18,3]] eps*=0.14, [[42,19,2]] eps*=0.04 |
| 67 | 21 | [[66,1,22]] | 0.742 | 1.12e-21 | 0.52 | [[51,16,7]] | 0.596 | 0.385 | [[66,1,22]] eps*=0.74, [[65,2,21]] eps*=0.73, [[64,3,20]] eps*=0.72, [[63,4,19]] eps*=0.70, [[62,5,18]] eps*=0.69, [[61,6,17]] eps*=0.67, [[60,7,16]] eps*=0.65, [[59,8,15]] eps*=0.63, [[58,9,14]] eps*=0.61, [[57,10,13]] eps*=0.59, [[56,11,12]] eps*=0.56, [[55,12,11]] eps*=0.54, [[54,13,10]] eps*=0.50, [[53,14,9]] eps*=0.47, [[52,15,8]] eps*=0.43, [[51,16,7]] eps*=0.38, [[50,17,6]] eps*=0.33, [[49,18,5]] eps*=0.28, [[48,19,4]] eps*=0.21, [[47,20,3]] eps*=0.13, [[46,21,2]] eps*=0.03 |
| 71 | 23 | [[70,1,24]] | 0.758 | 1.28e-24 | 0.49 | [[53,18,7]] | 0.555 | 0.377 | [[70,1,24]] eps*=0.76, [[69,2,23]] eps*=0.75, [[68,3,22]] eps*=0.74, [[67,4,21]] eps*=0.72, [[66,5,20]] eps*=0.71, [[65,6,19]] eps*=0.70, [[64,7,18]] eps*=0.68, [[63,8,17]] eps*=0.67, [[62,9,16]] eps*=0.65, [[61,10,15]] eps*=0.63, [[60,11,14]] eps*=0.61, [[59,12,13]] eps*=0.58, [[58,13,12]] eps*=0.56, [[57,14,11]] eps*=0.53, [[56,15,10]] eps*=0.50, [[55,16,9]] eps*=0.46, [[54,17,8]] eps*=0.42, [[53,18,7]] eps*=0.38, [[52,19,6]] eps*=0.33, [[51,20,5]] eps*=0.27, [[50,21,4]] eps*=0.20, [[49,22,3]] eps*=0.13, [[48,23,2]] eps*=0.03 |
| 73 | 23 | [[72,1,24]] | 0.750 | 1.52e-24 | 0.48 | [[56,17,8]] | 0.573 | 0.412 | [[72,1,24]] eps*=0.75, [[71,2,23]] eps*=0.74, [[70,3,22]] eps*=0.73, [[69,4,21]] eps*=0.72, [[68,5,20]] eps*=0.70, [[67,6,19]] eps*=0.69, [[66,7,18]] eps*=0.67, [[65,8,17]] eps*=0.66, [[64,9,16]] eps*=0.64, [[63,10,15]] eps*=0.62, [[62,11,14]] eps*=0.60, [[61,12,13]] eps*=0.57, [[60,13,12]] eps*=0.55, [[59,14,11]] eps*=0.52, [[58,15,10]] eps*=0.49, [[57,16,9]] eps*=0.45, [[56,17,8]] eps*=0.41, [[55,18,7]] eps*=0.37, [[54,19,6]] eps*=0.32, [[53,20,5]] eps*=0.26, [[52,21,4]] eps*=0.20, [[51,22,3]] eps*=0.12, [[50,23,2]] eps*=0.03 |
| 79 | 25 | [[78,1,26]] | 0.757 | 1.74e-27 | 0.46 | [[60,19,8]] | 0.553 | 0.395 | [[78,1,26]] eps*=0.76, [[77,2,25]] eps*=0.75, [[76,3,24]] eps*=0.74, [[75,4,23]] eps*=0.73, [[74,5,22]] eps*=0.71, [[73,6,21]] eps*=0.70, [[72,7,20]] eps*=0.69, [[71,8,19]] eps*=0.67, [[70,9,18]] eps*=0.66, [[69,10,17]] eps*=0.64, [[68,11,16]] eps*=0.62, [[67,12,15]] eps*=0.60, [[66,13,14]] eps*=0.58, [[65,14,13]] eps*=0.56, [[64,15,12]] eps*=0.53, [[63,16,11]] eps*=0.50, [[62,17,10]] eps*=0.47, [[61,18,9]] eps*=0.43, [[60,19,8]] eps*=0.40, [[59,20,7]] eps*=0.35, [[58,21,6]] eps*=0.30, [[57,22,5]] eps*=0.25, [[56,23,4]] eps*=0.19, [[55,24,3]] eps*=0.12, [[54,25,2]] eps*=0.03 |
| 83 | 27 | [[82,1,28]] | 0.770 | 1.43e-30 | 0.44 | [[62,21,8]] | 0.521 | 0.389 | [[82,1,28]] eps*=0.77, [[81,2,27]] eps*=0.76, [[80,3,26]] eps*=0.75, [[79,4,25]] eps*=0.74, [[78,5,24]] eps*=0.73, [[77,6,23]] eps*=0.72, [[76,7,22]] eps*=0.71, [[75,8,21]] eps*=0.70, [[74,9,20]] eps*=0.68, [[73,10,19]] eps*=0.67, [[72,11,18]] eps*=0.65, [[71,12,17]] eps*=0.63, [[70,13,16]] eps*=0.61, [[69,14,15]] eps*=0.59, [[68,15,14]] eps*=0.57, [[67,16,13]] eps*=0.55, [[66,17,12]] eps*=0.52, [[65,18,11]] eps*=0.49, [[64,19,10]] eps*=0.46, [[63,20,9]] eps*=0.43, [[62,21,8]] eps*=0.39, [[61,22,7]] eps*=0.35, [[60,23,6]] eps*=0.30, [[59,24,5]] eps*=0.24, [[58,25,4]] eps*=0.18, [[57,26,3]] eps*=0.11, [[56,27,2]] eps*=0.03 |
| 89 | 29 | [[88,1,30]] | 0.776 | 1.21e-33 | 0.41 | [[67,22,9]] | 0.507 | 0.412 | [[88,1,30]] eps*=0.78, [[87,2,29]] eps*=0.77, [[86,3,28]] eps*=0.76, [[85,4,27]] eps*=0.75, [[84,5,26]] eps*=0.74, [[83,6,25]] eps*=0.73, [[82,7,24]] eps*=0.72, [[81,8,23]] eps*=0.71, [[80,9,22]] eps*=0.70, [[79,10,21]] eps*=0.68, [[78,11,20]] eps*=0.67, [[77,12,19]] eps*=0.65, [[76,13,18]] eps*=0.64, [[75,14,17]] eps*=0.62, [[74,15,16]] eps*=0.60, [[73,16,15]] eps*=0.58, [[72,17,14]] eps*=0.56, [[71,18,13]] eps*=0.53, [[70,19,12]] eps*=0.51, [[69,20,11]] eps*=0.48, [[68,21,10]] eps*=0.45, [[67,22,9]] eps*=0.41, [[66,23,8]] eps*=0.37, [[65,24,7]] eps*=0.33, [[64,25,6]] eps*=0.29, [[63,26,5]] eps*=0.23, [[62,27,4]] eps*=0.18, [[61,28,3]] eps*=0.11, [[60,29,2]] eps*=0.03 |
| 97 | 31 | [[96,1,32]] | 0.775 | 1.05e-36 | 0.38 | [[73,24,9]] | 0.506 | 0.391 | [[96,1,32]] eps*=0.77, [[95,2,31]] eps*=0.77, [[94,3,30]] eps*=0.76, [[93,4,29]] eps*=0.75, [[92,5,28]] eps*=0.74, [[91,6,27]] eps*=0.73, [[90,7,26]] eps*=0.72, [[89,8,25]] eps*=0.71, [[88,9,24]] eps*=0.70, [[87,10,23]] eps*=0.69, [[86,11,22]] eps*=0.68, [[85,12,21]] eps*=0.66, [[84,13,20]] eps*=0.65, [[83,14,19]] eps*=0.63, [[82,15,18]] eps*=0.61, [[81,16,17]] eps*=0.60, [[80,17,16]] eps*=0.58, [[79,18,15]] eps*=0.56, [[78,19,14]] eps*=0.53, [[77,20,13]] eps*=0.51, [[76,21,12]] eps*=0.48, [[75,22,11]] eps*=0.46, [[74,23,10]] eps*=0.42, [[73,24,9]] eps*=0.39, [[72,25,8]] eps*=0.35, [[71,26,7]] eps*=0.31, [[70,27,6]] eps*=0.27, [[69,28,5]] eps*=0.22, [[68,29,4]] eps*=0.16, [[67,30,3]] eps*=0.10, [[66,31,2]] eps*=0.02 |

## Search and literature codes with small blocks (p <= 29)

| p | code | source | c (block unless noted) | eps* | notes |
|---|---|---|---|---|---|
| 3 | [[15,1,3]]_3 | doubling (2610.29457) | - | - | |
| 3 | [[19,1,3]]_3 | cyclic B (2610.08012) | 5.25 | 0.226 | |
| 3 | [[8,1,2]]_3 | RM | 2.0 | 0.211 | |
| 3 | [[37,1,5]]_3 | doubling | - | - | |
| 3 | [[72,9,3]]_3 | paper Table 3 | - | - | |
| 3 | [[206,37,4]]_3 | paper Table 3 | 15.2 (per output) | 0.0616 | |
| 5 | [[4,1,2]]_5 | RS | 1.5 | 0.363 | |
| 5 | [[11,1,3]]_5 | cyclic B | 1.44 | 0.346 | |
| 5 | [[16,1,4]]_5 | cyclic | 3.38 | 0.363 | |
| 5 | [[33,1,5]]_5 | doubling | - | - | |
| 5 | [[20,5,2]]_5 | paper Table 3 | - | - | |
| 5 | [[112,13,3]]_5 | paper Table 3 | - | - | |
| 7 | [[6,1,2]]_7 | RS | 2.5 | 0.232 | |
| 7 | [[7,1,3]]_7 | cyclic B | 0.97 | 0.446 | |
| 7 | [[12,1,3]]_7 | cyclic A | 1.11 | 0.366 | |
| 7 | [[20,1,4]]_7 | doubling | - | - | |
| 7 | [[44,5,3]]_7 | m = 2 search (catalog) | 1.11 | 0.173 | |
| 11 | [[103,18,4]]_11 | m = 2 search (catalog `imported-p11-n103-k18-48190679c64d`) | - | - | d exact (MITM certificate); exact curve not enumerable (dim A too large): P_acc = (1-eps)^n within 1% (COVERAGE.md), leading order needs the A_d census, rigorous bounds via `qmsd assess` |
| 11 | [[102,19,4]]_11 | m = 2 search (catalog `imported-p11-n102-k19-1dee9c45003b`) | - | - | d exact (MITM certificate); exact curve not enumerable (dim A too large): P_acc = (1-eps)^n within 1% (COVERAGE.md), leading order needs the A_d census, rigorous bounds via `qmsd assess` |
| 11 | [[101,20,4]]_11 | m = 2 search (catalog `imported-p11-n101-k20-2ef14bfb1b4d`) | - | - | d exact (MITM certificate); exact curve not enumerable (dim A too large): P_acc = (1-eps)^n within 1% (COVERAGE.md), leading order needs the A_d census, rigorous bounds via `qmsd assess` |
| 11 | [[100,21,4]]_11 | m = 2 search (catalog `imported-p11-n100-k21-20afd92fd535`) | - | - | d exact (MITM certificate); exact curve not enumerable (dim A too large): P_acc = (1-eps)^n within 1% (COVERAGE.md), leading order needs the A_d census, rigorous bounds via `qmsd assess` |
| 11 | [[99,22,3]]_11 | m = 2 search (catalog `imported-p11-n99-k22-1a07d37d9bc4`) | - | - | d exact (MITM certificate); exact curve not enumerable (dim A too large): P_acc = (1-eps)^n within 1% (COVERAGE.md), leading order needs the A_d census, rigorous bounds via `qmsd assess` |
| 11 | [[98,23,3]]_11 | m = 2 search (catalog `imported-p11-n98-k23-50d647de3876`) | - | - | d exact (MITM certificate); exact curve not enumerable (dim A too large): P_acc = (1-eps)^n within 1% (COVERAGE.md), leading order needs the A_d census, rigorous bounds via `qmsd assess` |
| 11 | [[97,24,3]]_11 | m = 2 search (catalog `imported-p11-n97-k24-a67756e8a78c`) | - | - | d exact (MITM certificate); exact curve not enumerable (dim A too large): P_acc = (1-eps)^n within 1% (COVERAGE.md), leading order needs the A_d census, rigorous bounds via `qmsd assess` |
| 11 | [[95,26,2]]_11 | m = 2 search (catalog `imported-p11-n95-k26-0c17be2979eb`) | - | - | d exact (MITM certificate); exact curve not enumerable (dim A too large): P_acc = (1-eps)^n within 1% (COVERAGE.md), leading order needs the A_d census, rigorous bounds via `qmsd assess` |
| 13 | [[147,22,5]]_13 | m = 2 search (catalog `imported-p13-n147-k22-0645a2561edb`) | - | - | d exact (MITM certificate); exact curve not enumerable (dim A too large): P_acc = (1-eps)^n within 1% (COVERAGE.md), leading order needs the A_d census, rigorous bounds via `qmsd assess` |
| 13 | [[145,24,4]]_13 | m = 2 search (catalog `imported-p13-n145-k24-233fc9c338f5`) | - | - | d exact (MITM certificate); exact curve not enumerable (dim A too large): P_acc = (1-eps)^n within 1% (COVERAGE.md), leading order needs the A_d census, rigorous bounds via `qmsd assess` |
| 13 | [[144,25,4]]_13 | m = 2 search (catalog `imported-p13-n144-k25-27375004be9c`) | - | - | d exact (MITM certificate); exact curve not enumerable (dim A too large): P_acc = (1-eps)^n within 1% (COVERAGE.md), leading order needs the A_d census, rigorous bounds via `qmsd assess` |
| 13 | [[143,26,4]]_13 | m = 2 search (catalog `imported-p13-n143-k26-1b214cf0a034`) | - | - | d exact (MITM certificate); exact curve not enumerable (dim A too large): P_acc = (1-eps)^n within 1% (COVERAGE.md), leading order needs the A_d census, rigorous bounds via `qmsd assess` |
| 13 | [[142,27,4]]_13 | m = 2 search (catalog `imported-p13-n142-k27-292c2753dcb5`) | - | - | d exact (MITM certificate); exact curve not enumerable (dim A too large): P_acc = (1-eps)^n within 1% (COVERAGE.md), leading order needs the A_d census, rigorous bounds via `qmsd assess` |
| 13 | [[141,28,4]]_13 | m = 2 search (catalog `imported-p13-n141-k28-23159eb36bae`) | - | - | d exact (MITM certificate); exact curve not enumerable (dim A too large): P_acc = (1-eps)^n within 1% (COVERAGE.md), leading order needs the A_d census, rigorous bounds via `qmsd assess` |
| 13 | [[139,30,3]]_13 | m = 2 search (catalog `imported-p13-n139-k30-113c64c72bc7`) | - | - | d exact (MITM certificate); exact curve not enumerable (dim A too large): P_acc = (1-eps)^n within 1% (COVERAGE.md), leading order needs the A_d census, rigorous bounds via `qmsd assess` |
| 13 | [[137,32,3]]_13 | m = 2 search (catalog `imported-p13-n137-k32-d87b912ef43a`) | - | - | d exact (MITM certificate); exact curve not enumerable (dim A too large): P_acc = (1-eps)^n within 1% (COVERAGE.md), leading order needs the A_d census, rigorous bounds via `qmsd assess` |
| 17 | [[237,52,6]]_17 | m = 2 search (catalog `new-p17-d6`) | - | - | d exact (MITM certificate); exact curve not enumerable (dim A too large): P_acc = (1-eps)^n within 1% (COVERAGE.md), leading order needs the A_d census, rigorous bounds via `qmsd assess` |
| 19 | [[299,62,6]]_19 | m = 2 search (catalog `imported-p19-n299-k62-177e7ce30c1d`) | - | - | d exact (MITM certificate); exact curve not enumerable (dim A too large): P_acc = (1-eps)^n within 1% (COVERAGE.md), leading order needs the A_d census, rigorous bounds via `qmsd assess` |
| 19 | [[296,65,6]]_19 | m = 2 search (catalog `imported-p19-n296-k65-211e1fbb4b1a`) | - | - | d exact (MITM certificate); exact curve not enumerable (dim A too large): P_acc = (1-eps)^n within 1% (COVERAGE.md), leading order needs the A_d census, rigorous bounds via `qmsd assess` |
| 19 | [[293,68,5]]_19 | m = 2 search (catalog `new-p19-d5`) | - | - | d exact (MITM certificate); exact curve not enumerable (dim A too large): P_acc = (1-eps)^n within 1% (COVERAGE.md), leading order needs the A_d census, rigorous bounds via `qmsd assess` |
| 23 | [[429,100,6]]_23 | m = 2 search (catalog `imported-p23-n429-k100-994d448ce6fc`) | - | - | d exact (MITM certificate); exact curve not enumerable (dim A too large): P_acc = (1-eps)^n within 1% (COVERAGE.md), leading order needs the A_d census, rigorous bounds via `qmsd assess` |
| 23 | [[424,105,6]]_23 | m = 2 search (catalog `imported-p23-n424-k105-9d5204d5ccf4`) | - | - | d exact (MITM certificate); exact curve not enumerable (dim A too large): P_acc = (1-eps)^n within 1% (COVERAGE.md), leading order needs the A_d census, rigorous bounds via `qmsd assess` |
| 29 | [[666,175,6]]_29 | m = 2 search (catalog `imported-p29-n666-k175-4f8ad426d19b`) | - | - | d exact (MITM certificate); exact curve not enumerable (dim A too large): P_acc = (1-eps)^n within 1% (COVERAGE.md), leading order needs the A_d census, rigorous bounds via `qmsd assess` |

## How to extend

- `scratchpad/rs_family_table.py` regenerates the RS table; `scratchpad/small_codes_build.py` rebuilds this file and the catalog/performance records.
- m = 2 searches at small k (high d, n ~ 100-300) for p = 11, 13 run on loranne (`runs/coverage/small/`); p >= 17 at d >= 6 needs the 3+3 kernel memory (lenore).
- Cost per output (raw states per accepted output at a target error) comes from `qmsd protocol` / `qmsd assess`; not yet tabulated here.
