# Independent review of the historical ozone simulations

All three integrations completed 1850–2014 from the common `icA` state. Independent checks covered 495 annual checkpoints, 5,940 model months and all 14,850 monthly NetCDF files. Each trajectory advanced 1,928,480 steps. The archived annual table, summary and figure agree with temperatures recomputed directly from the model output. [Annual temperatures](fig1/results/annual_temperature.csv) · [Summary](fig1/results/summary.json) · [Temperature figure](fig1/results/fig1_surface_temperature.png)

## Scientific validity

I recomputed each annual global near-surface air temperature directly from `tas`, using Gaussian latitude areas and Gregorian month lengths. I also recalculated Antarctic September–November temperatures in the top sigma layer by intersecting Gaussian latitude bands with the region south of 60°S. The table uses each experiment’s inclusive 1850–1900 mean as its baseline.

| Experiment | 1850–1900 baseline (K) | 2014 temperature (K) | 2014 anomaly (K) | 2005–2014 mean anomaly (K) | Antarctic spring top-layer change, 1980–2014 (K/decade) |
|---|---:|---:|---:|---:|---:|
| Historical CO₂ and ozone | 286.5453 | 287.2230 | +0.6777 | +0.6970 | -0.9011 |
| Historical ozone, fixed CO₂ | 286.4779 | 286.4123 | -0.0655 | -0.0996 | -0.3635 |
| Historical CO₂, fixed ozone | 286.5490 | 287.0202 | +0.4712 | +0.4285 | -0.5157 |

For 2005–2014 the mean raw joint-case difference relative to ozone-only was +0.8640 K and relative to CO₂-only was +0.2647 K. These differences compare single trajectories. Internal variability, nonlinear adjustment and drift can contribute. The shared initial ocean and atmosphere had been equilibrated with idealized ozone. Quantifying the interaction of both forcings requires an additional combination.

Reconstruction from the four native input4MIPs archives checked the spherical area integrals of 130,680 monthly pressure levels. Vertical transfer used the pressure overlap of each native pressure cell with the model sigma layer. The maximum matched-column difference across 4,055,040 target columns was 0.00002161 DU. The repeated seasonal file matched each of the twelve historical 1850 fields throughout. The L10 top layer held 74.6–81.0% of the global prescribed column. At that resolution the top-layer temperature averages a wide range of pressures. [Forcing validation and limitations](README.md#forcing-validation) · [Monthly anchor columns](fig2/results/monthly_columns.csv)

I also checked all 12,165,120 online monthly grid-cell columns against the bracketing prescribed ozone fields and the actual monthly surface pressure. Positive pressure bounds each month’s mean ozone mass by the minimum and maximum of the three applicable anchor profiles integrated over sigma and multiplied by mean pressure. The largest excess beyond this range was 0.000000 DU with a 0.02 DU allowance for floating-point accumulation. The online figure of merit tests the pressure weighting during the integration.

## Repository and submission

I verified a one-commit clean source candidate over upstream `f29ee8e` and the fork integration tree containing that source. The fork preserves its previous history. The two-rank T21/L10/O16 binary came from the integration tree used by the run script, giving all three experiments the same model. The ozone files, build and CLI switches, desktop option, monthly diagnostics, tests and documentation follow the repository’s feature entry points. [Production source](https://github.com/zsy1207/PlaSiC/commit/6f71e172f539831f66e29737e60057854f888dc9) · [Single-commit production source candidate](https://github.com/zsy1207/PlaSiC/commit/6fbb04ba0193f4c696f9824ea11a095ae7e17f04)

## Upstream maintainer requirements

The original PR specifies four historical UReading-CCMI-1-0 ozone files, three forcing combinations, the 1850–2014 period, and a comparison figure. I verified the four names and checksums, the current upstream sea ice and ocean baseline, and the three run commands. The ozone-only trajectory used 284.316853 ppmv CO₂, while CO₂-only repeated the external historical 1850 seasonal ozone cycle. All three used the same upstream warm restart. Annual values cover every year through 2014. Centred nine-year curves cover the 157 complete windows per case. [Maintainer comments](https://github.com/sunmoumou1/PlaSiC/pull/9) · [Source-file checksums](SHA256SUMS) · [Vector temperature figure](fig1/results/fig1_surface_temperature.pdf)

## Code, execution and results quality

The first draft of this rerun could deadlock if a two-rank library caller retried a valid short interval after an error in a longer interval. I reproduced the failure with an invalid March value on the southern rank. On the corrected production build both ranks rejected the extended interval at the same clock and successfully retried January. The 31-test fork discovery run passed. Two-step runs with the production binary verified the first 1850 radiation column of 287.615 DU and the expected distinction between historical and fixed CO₂.

I verified each incoming restart against its preceding output hash, each segment’s model and forcing identities, its calendar length and the manifest hashes of its restart and thirty monthly files. Between cases the incoming solar fields agreed and the two transient CO₂ fields agreed for all 165 years. I checked every annual absolute temperature, anomaly, paired difference and Antarctic spring temperature against the archived table, and recomputed the reported summary statistics. The printed PDF reproduces the 495 annual and 471 smoothed data points within the validation tolerance of 0.35 points, with matching curve and legend colours.

| Production item | SHA256 |
|---|---|
| MPI executable | `052be6d78dfa0adbc2d286db2daa999ba349e68359dc64468ab462aad53e3f85` |
| Initial state | `ebd4fc394b233a5fead4fd221ba842f3e4c51711a43b69a01c5a67f07eab8825` |
| Historical ozone | `ec881675b14bef25a3e586cf16f22ee83300e9ed1c4a729266fb6fd63c99a397` |
| Repeating 1850 ozone | `49b9092f36dec92853710c2483c84bee48aab222e7adffc1b2ad1917d989b62a` |
| Historical CO₂ and ozone final restart | `9528ea370a814fe7364dde49ae933095998bd1f139f0819ec9c0b1ece63a1b33` |
| Historical ozone, fixed CO₂ final restart | `dc14c71ab9d475c998e519c46483584db8e9302df778424d35764bae0cacc40e` |
| Historical CO₂, fixed ozone final restart | `bb5cdabf1d38ac01db4514333fbe46efa12c957b98239b52c5e554e383233cad` |
| fig1_surface_temperature.png | `d11b2a0d391da1eb99fe0f73914a736b4924729b8bd79421c78824e7023d2d2e` |
| fig1_surface_temperature.pdf | `bf57f995c4690738042d7a89a1c8fa59dcec9efcd96f4ff70022727b0b71c85a` |
| summary.json | `1b7cb02b8e9b6db8a814ff20ab8e9f11df0be4db38ec7d6035e7dc58d21a4939` |
| annual_temperature.csv | `c0379d54827ba58a49190c8e8184ed114b0b0367b4b6a6401e1edac6b7cb9746` |

Independent integrated review completed 2026-11-24 UTC.
