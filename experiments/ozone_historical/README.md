# Historical CO₂ and ozone forcing, 1850–2014

Three experiments compare annual global near-surface air temperature under changing CO₂ and ozone. They follow the [maintainer's five requests](https://github.com/sunmoumou1/PlaSiC/pull/9#issuecomment-7079546681) and use the grid, ocean and warm initial condition from the [reference historical experiment](../historical_3d_500y/).

| Experiment | CO₂ | Ozone |
|---|---|---|
| CO₂ + ozone | Historical UoM monthly concentrations | Historical UReading monthly fields |
| Ozone only | Fixed at 284.316853 ppmv | Historical UReading monthly fields |
| CO₂ only | Historical UoM monthly concentrations | Repeating external 1850 seasonal cycle |

Fixed CO₂ is a single global concentration. Fixed ozone retains the latitude, longitude and monthly pattern of 1850. All three start from `icA.restart`, spun up with idealized ozone. The common switch to external ozone can cause adjustment. Each configuration has one member. Own-series baselines remove a mean offset, but do not identify forced trends separately from initial adjustment, model drift or internal variability. The three trajectories describe the model under these configurations.

## Forcing data

The four original `vmro3` files are from CMIP6 input4MIPs, UReading-CCMI-1-0, native grid `gn`, version `v20160711`. They cover 1850–1899, 1900–1949, 1950–1999 and 2000–2014. Together they contain 1,980 monthly fields on 66 pressure levels and a 96 × 144 horizontal grid, totalling 5,927,876,651 bytes. SHA256 values in [SHA256SUMS](SHA256SUMS) match the ORNL and Argonne catalogue replicas. [Dataset identifier](https://hdl.handle.net/21.14100/ab27892c-eaa5-3430-a3c5-f9b0549a556f) · [ESGF search requested by the maintainer](https://github.com/sunmoumou1/PlaSiC/pull/9#issuecomment-7089574072)

## Reproduce the three experiments

Run from the repository root. The pipeline verifies the four ESGF checksums, builds T21/L10 with a 16-level ocean and two MPI ranks, prepares the common pressure climatology, checks the converted ozone columns, integrates the three cases and writes their annual temperature CSV, summary and figures. It uses Python 3.10 or later with NumPy, SciPy, netCDF4 and Matplotlib, plus curl and the [model build dependencies](../../INSTALL.md).

```bash
python experiments/ozone_historical/run_pipeline.py --jobs 2
```

The default download and output directories are `experiments/ozone_historical/data` and `experiments/ozone_historical/runs`. Use `--data` and `--output` to choose another disk. `--model` reuses an existing two-rank executable. `--reference-ps` reuses an existing twelve-month pressure file. The runner and reference-pressure calculation both select a 45-minute time step. With current builds they pass `--step-minutes 45`, overriding the general T21 default of 30 minutes. The archived T21 executable predates the option and uses 45 minutes by default. The default pressure reference is the first Gregorian year from `icA.restart` with idealized ozone, fixed 1850 CO₂ and 32 steps per day. All three historical experiments reuse this monthly pressure climatology.

The individual stages can also be run separately. The [converter](prepare_ozone.py) takes the four UReading files, `--reference-ps` and `--output`. Add `--fixed-1850` for the repeating seasonal cycle. The [experiment runner](run_experiments.py) accepts `--historical-ozone` and `--fixed-ozone`. The [temperature analysis](fig1/compute_temperature.py) and [plot](fig1/plot_temperature.py) read completed runs from `PLASIC_OZONE_RUNS`.

For example, with the four archived inputs and the pressure reference in the default locations:

```bash
exp=experiments/ozone_historical
python "$exp/prepare_ozone.py" "$exp"/data/vmro3_input4MIPs_ozone_CMIP_UReading-CCMI-1-0_gn_*.nc \
  --reference-ps "$exp/data/reference/monthly/ps.nc" --output "$exp/data/ozone_historical.nc"
```

The launcher can be replaced in full, for example with `--mpirun 'srun -n 2'` under an appropriate Slurm allocation. Each case takes 1,928,480 steps. To exercise a short run with the experiment runner, set `--end-year 1850`. To continue it, rerun with the same input files and output directory. Each completed year is atomically published at `<case>/<year>/` with its restart, monthly NetCDF files and provenance manifest. The current format checks the SHA256 of the restart and every monthly output on resume. Older manifests verify the restart and temperature file. An incomplete or changed checkpoint stops the run.

The converted files retain the integer year and month read by the model and add a CF time coordinate with monthly bounds, sigma layer bounds, the actual reference-pressure grid and source-code fingerprints. The pressure coordinate is `p = ptop + sigma × (ps − ptop)`, with `ptop = 0 Pa`. Older prepared files retain their original inputs and can still be used by the model.

## Forcing validation

The source-grid and model-grid monthly ozone columns cover every input anchor from January 1850 through December 2014. The source-grid mean integrates the 66 pressure cells to the reference surface pressure interpolated onto the original grid, weighted by the areas of the source cells. A second calculation area-averages the source profiles over each target cell and integrates their original pressure cells to the target reference pressure. On each native pressure level the area-weighted global mixing ratio is checked before and after horizontal remapping. The target uses the stored model layer means and Gaussian area weights. All three columns use standard gravity and the model ozone density at STP. [Monthly columns (CSV)](fig2/results/monthly_columns.csv) · [Column calculation](fig2/compute_columns.py)

The global model column spans 278.240–306.733 DU. Its difference from the original-grid mean ranges from −0.00505 to −0.00097 DU. At matched horizontal cell areas and reference pressure the maximum difference over all 4,055,040 grid columns is 0.00002160 DU, including float32 storage.

![Prescribed zonal-mean ozone column by latitude and calendar month](fig2/results/zonal_seasonal_cycle.png)

The figure averages the prescribed column over longitudes and over the 165 occurrences of each calendar month. The ozone column peaks in northern spring at high latitudes. [Zonal climatology (CSV)](fig2/results/zonal_seasonal_cycle.csv) · [Vector figure (PDF)](fig2/results/zonal_seasonal_cycle.pdf)

## Remapping and interpretation

- Source values are ozone mole fractions. The model converts mole to mass fraction using 47.9982/28.9644 for ozone/dry air. At each source pressure level the converter weights latitude cells by their overlap in sine of latitude and longitude cells by their periodic overlap. Gaussian quadrature weights define the target latitude-band areas. Vertically each pressure cell contributes in proportion to its overlap with the target layer. Use the deepest source value between the bottom source edge and local ground where necessary.
- Model interfaces use the same monthly reference surface pressure in every year. This conserves the vertical column on that reference pressure grid. The radiation calculation uses the model's instantaneous pressure. The converter reads the Gaussian horizontal grid from the pressure file. Set `--levels` to match atmospheric builds other than L10.
- Both external files use identical time interpolation. January 1850 holds its first value before mid-month, and December 2014 holds its last value after mid-month. In other months the model interpolates between adjacent monthly centres. The interpolated field generally has a different calendar-month average from its input anchor. The fixed file repeats the twelve 1850 monthly fields throughout. This retains the spatial and seasonal structure of the common initial ozone forcing. Reading separately lists a CMIP6 control ozone product. [CMIP6 provider catalogue](https://blogs.reading.ac.uk/ccmi/forcing-databases-in-support-of-cmip6/)
- At fixed October reference pressure south of 60°S, monthly averaging of the interpolated fields raises the mean annual October spatial minimum by 9.76% over 2000–2014. The corresponding increase in the polar-cap area mean is 1.80%. The October 16 midnight field uses the original anchor.
- In a 60-month sensitivity comparison, switching between bilinear sampling and area-conservative horizontal remapping changed the global column by at most 0.018 DU and individual cells by up to 3.70 DU. In the preliminary bilinear calculation, substituting a ten-year monthly pressure climatology from the CO₂-only case changed global columns by up to 0.080 DU and individual cells by up to 6.30 DU. These are offline sensitivities, with both methods evaluated at a common pressure in each comparison.
- The L10 top layer contains 74.6–81.0% of the global ozone column over the 1,980 months. Across all months, using layer-midpoint pressure in the longwave calculation overestimates the top-layer pressure-weighted ozone path by 42.9–48.8% and the full-column path by 11.1–15.5%. These percentages describe the optical-path approximation. They do not quantify the resulting flux or temperature errors.
- In 1851–1860, the global mean net downward top-of-atmosphere flux minus full-depth ocean heat storage was +4.86, +4.31, +4.23 W m⁻² for the combined, ozone-only and CO₂-only runs, respectively. This partial budget excludes storage in other reservoirs. A one-year reference with idealized ozone also showed net downward top-of-atmosphere flux and decreasing ocean heat content.

## Instantaneous radiation diagnostic

At the first radiation call on 1 January 1850 at 00:00, three ozone substitutions use the same `icA` meteorology, CO₂, solar geometry and cloud fields. The January 2014 field is evaluated on this same 1850 clock. Global means use Gaussian area weights, with positive flux downward at the top of the atmosphere. All other kernel inputs were checked for byte equality across substitutions on both MPI ranks.

| Ozone substitution | Net shortwave change (W m⁻²) | Net longwave change (W m⁻²) | Net change (W m⁻²) |
|---|---:|---:|---:|
| January 1850 external minus model idealized | −0.1644 | −0.9348 | −1.0992 |
| January 2014 external minus January 1850 external | −0.0627 | +0.0932 | +0.0305 |

A separate offline longwave calculation integrates the first pressure moment of ozone across the published source cells within each L10 layer. It replaces the midpoint approximation only for the longwave ozone path, holding other radiation inputs fixed. Correcting the top layer alone reduces net downward longwave flux by 0.5297 W m⁻² in the 1850 snapshot and 0.4983 W m⁻² for the 2014 field. Correcting all layers gives a 2014-minus-1850 longwave difference of +0.1295 W m⁻² and a net difference of +0.0668 W m⁻². This diagnostic quantifies one vertical approximation in one atmospheric state. The 165-year coupled experiments use the standard L10 radiation scheme.

These instantaneous differences cannot identify a tropospheric or stratospheric contribution to the temperature trends. Checa-Garcia et al. report +0.30 ± 0.17 W m⁻² for a different diagnostic using 1850–1859 and 2010–2014 climatologies, 60 vertical levels and stratospheric temperature adjustment. [Checa-Garcia et al. (2018)](https://doi.org/10.1002/2017GL076770)

## Temperature analysis and validation

Use the model’s near-surface air temperature (`tas`), Gaussian quadrature in latitude, equal longitude weights and actual calendar-month lengths. The CSV contains each annual mean in kelvin and its anomaly relative to that case's inclusive 1850–1900 mean. The figure shows annual values and a separate centred nine-year mean for each case. Means span 1854–2010 with no endpoint padding. Model diagnostics, the CSV and PNG/PDF are written under `$PLASIC_OZONE_RUNS`. The summary JSON records absolute 2014 temperature, the 1850–1900 mean and both 2014 and 2005–2014 anomalies for each case. The CSV also gives raw annual paired differences between the combined case and each single-forcing case. These are conditional comparisons between single trajectories. A common additive offset cancels, but nonlinear adjustment and internal variability may differ. With three cases the CO₂–ozone interaction is unidentifiable.

For the archived atmospheric temperature field `ta`, the analysis averages the top layer (σ = 0 to 0.0766) south of 60°S over September–November. Gaussian latitude-band areas are clipped at 60°S and months are weighted by duration. The CSV gives annual spring values, and the summary compares 1960–1979, 1980–1999 and 2000–2014 and reports the 1980–2014 linear slope. This is a layer-mean temperature diagnostic for one model realization. The reference historical plotting code and these calculations agree on annual weighting and the full nine-year windows. [Reference plotting script](../historical_3d_500y/plots/plot_historical_surface_temperature.py)

## Historical temperature results

Each trajectory covers January 1850 through December 2014. Anomalies use the corresponding case’s 1850–1900 mean. The table uses the model’s near-surface air temperature diagnostic (`tas`). The production runs used 45-minute time steps (32 steps per day), as did the published reference experiment. [Production revision](https://github.com/zsy1207/PlaSiC/commit/6f71e172f539831f66e29737e60057854f888dc9) · [Upstream time-step update made during the runs](https://github.com/sunmoumou1/PlaSiC/commit/9568fa9fdd6ffb4e5a9c4a811b0e4e8216a2c5de)

After the upstream time-step update, the rebuilt model replayed the combined case for 1850 from the common initial state and the ozone-only case for 2014 from its archived end-of-2013 restart. At 45 minutes, all 30 monthly NetCDF files and the final restart were byte-identical to the production archive in each year.

| Experiment | 1850–1900 mean (K) | 2005–2014 mean anomaly (K) | 2014 anomaly (K) |
|---|---:|---:|---:|
| CO₂ + ozone | 286.545 | +0.697 | +0.678 |
| Ozone only | 286.478 | −0.100 | −0.066 |
| CO₂ only | 286.549 | +0.429 | +0.471 |

The 2005–2014 mean of the raw annual combined-minus-CO₂-only temperature difference was +0.265 K. The corresponding combined-minus-ozone-only difference was +0.864 K. Each difference compares two single trajectories.

![Annual global near-surface air temperature anomalies under three historical forcing combinations](fig1/results/fig1_surface_temperature.png)

Thin curves show all 165 annual values. Thick curves show complete centred nine-year means. [Annual values (CSV)](fig1/results/annual_temperature.csv) · [Statistics (JSON)](fig1/results/summary.json) · [Vector figure (PDF)](fig1/results/fig1_surface_temperature.pdf) · [Independent review](FINAL_REVIEW.md)

```bash
make -C src test-data test-ozone test-ozone-mpi
```

To test the annual restart boundary, run `python experiments/ozone_historical/verify_restart.py` after completing the combined case through 1851. It replays both years continuously with the recorded inputs and checks every stored monthly array and the complete binary end restart against the two annual segments. Set `PLASIC_OZONE_RUNS` or pass `--output` to select the archive. With the actual joint-forcing inputs, the 23,360-step continuous 1850–1851 integration matched both segmented years across 30 monthly files and all 132 records of the binary end restart (SHA256 `05788b87649ebc2ebc107633561ffea67b8d54db78770f5bc3b41fed93feccaf`).

Both ozone schemes write `monthly/toz.nc` with active radiation. This field integrates the ozone mass in each layer using the instantaneous surface pressure, then averages the column over all steps in a complete calendar month. The stored unit of $10^{-5}$ m at STP equals one Dobson unit. Each file records the ozone source. For synthetic ozone the independent January 1850 check reproduced all 2,048 grid cells within 0.00065 DU. Its area-weighted monthly mean was 305.0429 DU. A separate external-field test prescribes a global mean of 300 DU.

For prescribed ozone each model invocation reports the first radiation step's Gaussian-area-weighted column in DU. Startup scans all monthly inputs that can affect the requested run and reports suspicious ozone magnitudes on any MPI latitude slab. The ozone-specific tests cover remapping and native-grid units, calendar edges, startup validation, two-rank error propagation and radiation responses. The general upstream `test-runtime` target in `f29ee8e` references a missing test script. GCC 14 also issues a pre-existing `-Wmaybe-uninitialized` warning in the ocean halo code under strict builds. The default model and ozone tests build without promoting that upstream warning to an error.
