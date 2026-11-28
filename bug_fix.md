# PlaSiC Bug Fix Log

## 1. The sea-ice concentration was reset to the binary ice mask at every step

**Fixed on** 2026-11-21 19:43:31 (commit `09d27c9`)

### Original issue

At the start of `ice_step()`, a `memcpy()` copied `ice.ice_cover`, the previous step's 0/1 ice mask, into `ice.diagnosed_ice_cover`, the continuous compactness $`c_i^*`$. At a grid point below the 0.5 threshold, this reset $`c_i^*`$ to zero even when ice remained. Small increases could not accumulate across steps, so the ice edge could retreat but struggled to advance. In the old T42/L16 restart, only 50 of 800 sea grid points with positive ice thickness were inside the ice mask. The old T85/L25 restart had no points inside the mask.

### Fix

The model now carries the continuous value in `ice.compactness` across time steps. A cold start initializes it once from the climatological ice cover. The restart reader and writer use `xicecc` to preserve it across runs. At each step, `ice_compactness()` updates $`c_i^*`$ from the change in ice thickness. In prognostic ice mode, `ice_step()` then derives the binary mask from the updated value, assigning 1 when $`c_i^* \ge 0.5`$ and 0 otherwise. This lets small gains in compactness accumulate until a grid point crosses the threshold.


## 2. Reflected longwave radiation was missing from the surface flux

**Fixed on** 2026-11-28 11:20:10

### Original issue

The code left out reflected radiation when calculating the upward longwave flux at the surface. As a result, the net longwave flux into the surface was too high.

### Fix

The code first stores the surface emissivity in `zeps[horizontal]`, then overwrites it with the reflected longwave flux $`(1-\varepsilon)F_\downarrow`$. Because upward flux is negative in `rad_longwave()`, the fix includes the full reflected flux at the surface by executing `upward_flux[surface] -= zeps[horizontal];`. The corrected net surface longwave flux is

```math
F_{\mathrm{net},s}
= F_\downarrow-\varepsilon B_s-(1-\varepsilon)F_\downarrow
= \varepsilon(F_\downarrow-B_s).
```

Here, $`F_\downarrow`$ is the downward longwave flux at the surface, $`\varepsilon`$ is the surface emissivity, and $`B_s=\sigma_{\mathrm{SB}}T_s^4`$ is the blackbody flux at the surface temperature.

## 3. Dissipative heating from vertical mixing counted only half the kinetic energy transport

**Fixed on** 2026-11-28 11:36:37

### Original issue

In the `ndheat > 0` branch of `vdiff()`, `zu` and `zv` are the wind components before vertical diffusion, and `zun` and `zvn` are the wind components computed by the implicit diffusion solver. The variable `zke` stores the kinetic energy per unit mass, calculated as `0.5f * zabsu[index] * zabsu[index]`. The variable `zken` is the result of applying scalar diffusion to `zke` with the momentum diffusion coefficient `zkdiffm`, so `zken - zke` is the change in kinetic energy due to transport between layers.

Let $`\Delta S`$ be `zun² + zvn² - zu² - zv²`, $`\Delta K_{\mathrm{tr}}`$ be `zken - zke`, and $`\Delta t`$ be `deltsec2`. The original code took the difference between `wind_energy_change` ($`\Delta S/\Delta t`$) and `transported_energy_change` ($`\Delta K_{\mathrm{tr}}/\Delta t`$), then multiplied the entire difference by $`1/2`$. The factor of $`1/2`$ converts $`\Delta S`$ to the change in kinetic energy of the wind field. The factor of $`1/2`$ in $`\Delta K_{\mathrm{tr}}`$ comes from the definition of kinetic energy. The original expression therefore counted only half of the kinetic energy transported between layers and miscalculated the energy contributing to dissipative heating in each layer.

### Fix

The dissipative heating rate is the kinetic energy tendency due to transport minus the actual kinetic energy tendency of the wind field. Dividing this rate by the specific heat capacity of moist air, $`c_p(q)`$, gives the temperature tendency `zdtdt[index]`, which is then added to `dtdt[index]`. In the code, $`c_p(q)`$ is `cpd * (1.0f + cpv_cpd_minus1 * dq[index])`.

```math
\left(\frac{dT}{dt}\right)_{\mathrm{diss}}
= \frac{\Delta K_{\mathrm{tr}}-\tfrac12\Delta S}{c_p(q)\,\Delta t}.
```

## 4. The pressure vertical velocity calculation used the wrong discretization

**Fixed on** 2026-11-28 11:41:26

### Original issue

`pressure_vertical_velocity` should store $`\omega=dp/dt`$ ($`\mathrm{Pa\,s^{-1}}`$) at layer centers. The original code used $`p_s\Omega[P_k-(C(D+P))_k]`$ for the interior and bottom layers.

- $`P_k-(C(D+P))_k`$ is the layer average of $`\omega/(p\Omega)`$ used for compressional heating. Using this average in place of the layer-center value introduces a discretization error.
- Converting the layer-center value of $`\omega/(p\Omega)`$ to $`\omega`$ requires multiplying by the local pressure $`p_k=\sigma_kp_s`$ and by $`\Omega`$. Using $`p_s`$ in place of the local pressure introduces a separate factor of $`1/\sigma_k`$, which amplifies the result more strongly near the model top.

The top layer used a separate expression, $`+p_s\Omega D_1`$, which had the wrong sign and omitted the half-layer thickness $`\Delta\sigma_1/2`$.

### Fix

In $`\sigma_kP_k-\int_0^{\sigma_k}(D+P)\,d\sigma`$, split the integral into the full layers above the current layer and the upper half of the current layer. Combining the $`P_k`$ terms gives the layer-center formula

```math
\frac{\omega_k}{p_s\Omega}
=\sigma_{k-1/2}P_k
-\sum_{j<k}\Delta\sigma_j(D_j+P_j)
-\frac{\Delta\sigma_k}{2}D_k.
```

Here, $`\sigma_{k-1/2}`$ is the upper interface of the current layer. The upper interface of the top layer is at zero, and there are no full layers above it, so $`\omega_1=-p_s\Omega(\Delta\sigma_1/2)D_1`$. The matrix $`C`$ is still used for compressional heating, and the parameter was renamed to `compression_matrix`.

