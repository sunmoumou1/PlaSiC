#include "plasicmod_state.h"

#include <stddef.h>

// Fixed physical constants shared by the whole model. They are const and are
// therefore not touched by plasic_state_init(); keeping them here gives every
// global model quantity a single home.
const float plasic_two_pi = 6.28318530718f;
const float plasic_water_gas_constant = 461.51f;
const float plasic_water_heat_capacity = 1869.46f;

int32_t plasic_nstep = 0;
// Active calendar kind (0 = Gregorian, 1 = 12x30-day 360-day calendar).
// It decides how nstep maps to a date, is restored from the restart file on a
// warm start, and is always saved back into the restart file.
int32_t plasic_calendar = 0;
float plasic_mpstep = 0; // minutes per timestep (may be fractional, e.g. 7.5)
int32_t plasic_ntspd = 0;

// Current model date & time array:
// [0] year, [1] month, [2] day, [3] hour, [4] minute,
// [5] weekday (0-6), [6] leap-year flag (0 = common, 1 = leap)
int32_t plasic_ndatim[(7)] = {0};

// Accumulation interval in timesteps (0 = once per day): the daily-mean
// accumulators are reset every plasic_nafter steps; set to plasic_ntspd at init.
int32_t plasic_nafter = 0;

// Number of initial (startup) timesteps using explicit Euler; 0 means the model
// has entered the normal leapfrog integration phase.
int32_t plasic_nkits = 0;

// Switch for radiation: 1 = enabled, 0 = disabled.
int32_t plasic_nrad = 0;
// Switch for advection: 1 = enabled, 0 = disabled.
int32_t plasic_nadv = 0;
// Use equidistant sigma levels: 1 = yes, 0 = no.
int32_t plasic_neqsig = 0;
// Running counter of samples accumulated in the current output window;
// outputs are divided by this value to form time means.
int32_t plasic_naccuout = 0;
// Critical wavenumber for horizontal diffusion.
int32_t plasic_nhdiff = 0;
// Switch for heating due to momentum dissipation: 1 = enabled, 0 = disabled.
int32_t plasic_ndheat = 0;
// Switch for the top-of-model sponge layer: 1 = enabled, 0 = disabled.
int32_t plasic_nsponge = 0;

// Latent heat of sublimation (J/kg).
float plasic_ls = 0;
// Latent heat of vaporization (J/kg).
float plasic_lv = 0;
// Planetary vorticity (Coriolis parameter at the reference latitude, 1/s).
float plasic_plavor = 0;
// Timestep in seconds (solar_day / ntspd).
float plasic_deltsec = 0;
// Twice the timestep in seconds (2 * deltsec), used by leapfrog physics updates.
float plasic_deltsec2 = 0;
// Non-dimensional timestep = deltsec * Omega (rotation rate); computed as 2*pi/ntspd.
float plasic_delt = 0;
// Twice the non-dimensional timestep (2 * delt).
float plasic_delt2 = 0;
// Ground temperature of the mean atmospheric profile (K).
float plasic_tgr = 0;
// Global mean surface pressure (Pa).
float plasic_psurf = 0;
// Melting point of water (K).
float plasic_tmelt = 0;
// Damping time in days for the top sponge layer.
float plasic_dampsp = 0;

float plasic_sd[(PLASIC_NRSP) * (PLASIC_NLEV)] = {0};
float plasic_st[(PLASIC_NRSP) * (PLASIC_NLEV)] = {0};
float plasic_sz[(PLASIC_NRSP) * (PLASIC_NLEV)] = {0};
float plasic_sq[(PLASIC_NRSP) * (PLASIC_NLEV)] = {0};
float plasic_sp[(PLASIC_NRSP)] = {0};
float plasic_so[(PLASIC_NRSP)] = {0};
float plasic_sr[(PLASIC_NRSP) * (PLASIC_NLEV)] = {0};
float plasic_sdp[(PLASIC_NSPP) * (PLASIC_NLEV)] = {0};
float plasic_stp[(PLASIC_NSPP) * (PLASIC_NLEV)] = {0};
float plasic_szp[(PLASIC_NSPP) * (PLASIC_NLEV)] = {0};
float plasic_sqp[(PLASIC_NSPP) * (PLASIC_NLEV)] = {0};
float plasic_spp[(PLASIC_NSPP)] = {0};
float plasic_sop[(PLASIC_NSPP)] = {0};
float plasic_srp[(PLASIC_NSPP) * (PLASIC_NLEV)] = {0};
float plasic_sdt[(PLASIC_NSPP) * (PLASIC_NLEV)] = {0};
float plasic_stt[(PLASIC_NSPP) * (PLASIC_NLEV)] = {0};
float plasic_szt[(PLASIC_NSPP) * (PLASIC_NLEV)] = {0};
float plasic_sqt[(PLASIC_NSPP) * (PLASIC_NLEV)] = {0};
float plasic_spt[(PLASIC_NSPP)] = {0};
float plasic_sdm[(PLASIC_NSPP) * (PLASIC_NLEV)] = {0};
float plasic_stm[(PLASIC_NSPP) * (PLASIC_NLEV)] = {0};
float plasic_szm[(PLASIC_NSPP) * (PLASIC_NLEV)] = {0};
float plasic_sqm[(PLASIC_NSPP) * (PLASIC_NLEV)] = {0};
float plasic_spm[(PLASIC_NSPP)] = {0};

// Scale-selective spectral diffusion factor L_{n,k} over the full spectral
// array (layout: [mode + level * PLASIC_NRSP]); it is 0 below the critical
// wavenumber plasic_nhdiff and ramps to 1 at the truncation.
float plasic_Lnk[(PLASIC_NRSP) * (PLASIC_NLEV)] = {0};
// Per-rank copy of plasic_Lnk containing only the PLASIC_NSPP local spectral
// coefficients owned by this MPI rank (layout: [local_mode + level * PLASIC_NSPP]);
// consumed by spectrald_dissipation during time stepping.
float plasic_Lnkpp[(PLASIC_NSPP) * (PLASIC_NLEV)] = {0};

// Per-spectral-mode normalization factors (+/-1/sqrt(2), with the sign
// alternating with zonal wavenumber m) for each Re/Im coefficient pair.
float plasic_spnorm[(PLASIC_NRSP)] = {0};
// Total spherical wavenumber n of each spectral mode (both Re/Im slots share
// the same n); used for the 1/[n(n+1)] Laplacian term in the semi-implicit solver.
// The array is padded to PLASIC_NESP so that the last rank can index its whole
// NSPP-slot block (including padding) without reading past the end.
int32_t plasic_nindex[(PLASIC_NESP)] = {0};
// Per-level order (exponent) of the scale-selective horizontal diffusion:
// ((n - plasic_nhdiff) / (N - plasic_nhdiff))^ndel[level].
int32_t plasic_ndel[(PLASIC_NLEV)] = {0};

float plasic_gd[(PLASIC_NHOR) * (PLASIC_NLEV)] = {0};
float plasic_gt[(PLASIC_NHOR) * (PLASIC_NLEV)] = {0};
float plasic_gz[(PLASIC_NHOR) * (PLASIC_NLEV)] = {0};
float plasic_gq[(PLASIC_NHOR) * (PLASIC_NLEV)] = {0};
float plasic_gu[(PLASIC_NHOR) * (PLASIC_NLEV)] = {0};
float plasic_gv[(PLASIC_NHOR) * (PLASIC_NLEV)] = {0};
float plasic_gtdt[(PLASIC_NHOR) * (PLASIC_NLEV)] = {0};
float plasic_gqdt[(PLASIC_NHOR) * (PLASIC_NLEV)] = {0};
float plasic_gudt[(PLASIC_NHOR) * (PLASIC_NLEV)] = {0};
float plasic_gvdt[(PLASIC_NHOR) * (PLASIC_NLEV)] = {0};

// Grid-point log surface pressure anomaly pi = ln(ps / plasic_psurf),
// synthesized from the spectral surface pressure plasic_sp.
float plasic_gp[(PLASIC_NHOR)] = {0};
// Robert-form (cos-latitude scaled) northward gradient of ln(ps).
float plasic_gpj[(PLASIC_NHOR)] = {0};

float plasic_rcsq[(PLASIC_NHOR)] = {0};
float plasic_dt[(PLASIC_NHOR) * (PLASIC_NLEP)] = {0};
float plasic_dq[(PLASIC_NHOR) * (PLASIC_NLEP)] = {0};
float plasic_du[(PLASIC_NHOR) * (PLASIC_NLEP)] = {0};
float plasic_dv[(PLASIC_NHOR) * (PLASIC_NLEP)] = {0};
float plasic_dp[(PLASIC_NHOR)] = {0};
float plasic_dqsat[(PLASIC_NHOR) * (PLASIC_NLEP)] = {0};

// Humidity changes caused by dynamic advection
float plasic_dqt[(PLASIC_NHOR) * (PLASIC_NLEP)] = {0};
float plasic_dcc[(PLASIC_NHOR) * (PLASIC_NLEP)] = {0};
float plasic_dql[(PLASIC_NHOR) * (PLASIC_NLEP)] = {0};

// pressure_velocity Pa s-1
float plasic_dw[(PLASIC_NHOR) * (PLASIC_NLEV)] = {0};

float plasic_dtdt[(PLASIC_NHOR) * (PLASIC_NLEP)] = {0};
float plasic_dqdt[(PLASIC_NHOR) * (PLASIC_NLEP)] = {0};
float plasic_dudt[(PLASIC_NHOR) * (PLASIC_NLEP)] = {0};
float plasic_dvdt[(PLASIC_NHOR) * (PLASIC_NLEP)] = {0};

// Dimensional gridpoint snapshots of the state at time t, saved before the
// physics tendencies modify the working fields (gridpointa_snapshots):
// surface pressure ps = psurf * (ps/psurf) [Pa].
float plasic_dp0[(PLASIC_NHOR)] = {0};
// Zonal wind u = cv * gu / cos(lat) [m/s] at time t.
float plasic_du0[(PLASIC_NHOR) * (PLASIC_NLEP)] = {0};
// Meridional wind v = cv * gv / cos(lat) [m/s] at time t.
float plasic_dv0[(PLASIC_NHOR) * (PLASIC_NLEP)] = {0};

float plasic_dalb[(PLASIC_NHOR)] = {0};
float plasic_dswfl[(PLASIC_NHOR) * (PLASIC_NLEP)] = {0};
float plasic_dlwfl[(PLASIC_NHOR) * (PLASIC_NLEP)] = {0};
float plasic_dflux[(PLASIC_NHOR) * (PLASIC_NLEP)] = {0};

// Shortwave (solar) radiation flux in W/m2, positive downward; upward-going
// flux is stored negative, so the net flux is dfu + dfd.
float plasic_dfu[(PLASIC_NHOR) * (PLASIC_NLEP)] = {0};
float plasic_dfd[(PLASIC_NHOR) * (PLASIC_NLEP)] = {0};

// Longwave (thermal) radiation flux in W/m2, positive downward; upward-going
// flux is stored negative, so the net flux is dftu + dftd.
float plasic_dftu[(PLASIC_NHOR) * (PLASIC_NLEP)] = {0};
float plasic_dftd[(PLASIC_NHOR) * (PLASIC_NLEP)] = {0};

float plasic_dwetfac[(PLASIC_NHOR)] = {0};
float plasic_dls[(PLASIC_NHOR)] = {0};
// Aerodynamic roughness length z0 (m), used by the surface flux formulation.
float plasic_dz0[(PLASIC_NHOR)] = {0};

float plasic_diced[(PLASIC_NHOR)] = {0};
float plasic_dicec[(PLASIC_NHOR)] = {0};
float plasic_dtaux[(PLASIC_NHOR)] = {0};
float plasic_dtauy[(PLASIC_NHOR)] = {0};
float plasic_dust3[(PLASIC_NHOR)] = {0};
float plasic_dshfl[(PLASIC_NHOR)] = {0};
float plasic_dlhfl[(PLASIC_NHOR)] = {0};
float plasic_devap[(PLASIC_NHOR)] = {0};
float plasic_dtsa[(PLASIC_NHOR)] = {0};
float plasic_dmld[(PLASIC_NHOR)] = {0};
float plasic_dshdt[(PLASIC_NHOR)] = {0};
float plasic_dlhdt[(PLASIC_NHOR)] = {0};
float plasic_dprc[(PLASIC_NHOR)] = {0};
float plasic_dprl[(PLASIC_NHOR)] = {0};
float plasic_dprs[(PLASIC_NHOR)] = {0};
float plasic_dqvi[(PLASIC_NHOR)] = {0};

// Maximum soil water holding capacity of the surface bucket (field capacity),
// in kg/m2.
float plasic_dwmax[(PLASIC_NHOR)] = {0};
// Actual soil water content of the surface bucket, in kg/m2.
float plasic_dwatc[(PLASIC_NHOR)] = {0};
// Surface snow depth expressed as snow water equivalent, in kg/m2
float plasic_dsnow[(PLASIC_NHOR)] = {0};
// Snow melt rate, in kg/m2/s, positive when snow is melting.
float plasic_dsmelt[(PLASIC_NHOR)] = {0};
// Tendency of snow depth (snow water equivalent), in kg/m2/s.
float plasic_dsndch[(PLASIC_NHOR)] = {0};

float plasic_dtsoil[(PLASIC_NHOR)] = {0};
// Diagnostic soil temperature of soil layers 2 to 5, in K.
float plasic_dtd2[(PLASIC_NHOR)] = {0};
float plasic_dtd3[(PLASIC_NHOR)] = {0};
float plasic_dtd4[(PLASIC_NHOR)] = {0};
float plasic_dtd5[(PLASIC_NHOR)] = {0};
// Glacier mask: 1 for glacier-covered land points, 0 otherwise.
float plasic_dglac[(PLASIC_NHOR)] = {0};

float plasic_aevap[(PLASIC_NHOR)] = {0};
float plasic_aprl[(PLASIC_NHOR)] = {0};
float plasic_aprc[(PLASIC_NHOR)] = {0};
float plasic_aprs[(PLASIC_NHOR)] = {0};
float plasic_ashfl[(PLASIC_NHOR)] = {0};
float plasic_alhfl[(PLASIC_NHOR)] = {0};
float plasic_asmelt[(PLASIC_NHOR)] = {0};
float plasic_asndch[(PLASIC_NHOR)] = {0};
float plasic_acc[(PLASIC_NHOR)] = {0};

// Accumulated net surface shortwave (solar) radiation flux, in W/m2.
float plasic_assol[(PLASIC_NHOR)] = {0};
// Accumulated net surface longwave (thermal) radiation flux, in W/m2.
float plasic_asthr[(PLASIC_NHOR)] = {0};
// Accumulated net top-of-atmosphere shortwave radiation flux, in W/m2.
float plasic_atsol[(PLASIC_NHOR)] = {0};
// Accumulated net top-of-atmosphere longwave radiation flux, in W/m2.
float plasic_atthr[(PLASIC_NHOR)] = {0};
// Accumulated surface upward shortwave radiation flux, in W/m2.
float plasic_assolu[(PLASIC_NHOR)] = {0};
// Accumulated surface upward longwave radiation flux, in W/m2.
float plasic_asthru[(PLASIC_NHOR)] = {0};
// Accumulated top-of-atmosphere upward shortwave radiation flux, in W/m2.
float plasic_atsolu[(PLASIC_NHOR)] = {0};

float plasic_ataux[(PLASIC_NHOR)] = {0};
float plasic_atauy[(PLASIC_NHOR)] = {0};
float plasic_aqvi[(PLASIC_NHOR)] = {0};
float plasic_atsa[(PLASIC_NHOR)] = {0};
float plasic_atsama[(PLASIC_NHOR)] = {0};
float plasic_atsami[(PLASIC_NHOR)] = {0};
float plasic_ats0[(PLASIC_NHOR)] = {0};

// Sine of the Gaussian grid latitude, sin(phi), per latitude band.
double plasic_sid[(PLASIC_NLAT)] = {0};
// Gaussian quadrature weights associated with each latitude band.
double plasic_gwd[(PLASIC_NLAT)] = {0};
// Square of the cosine of latitude, cos^2(phi).
float plasic_csq[(PLASIC_NLAT)] = {0};
// Cosine of latitude, cos(phi) = sqrt(csq).
float plasic_cola[(PLASIC_NLAT)] = {0};
// Reciprocal of the cosine of latitude, 1/cos(phi).
float plasic_rcs[(PLASIC_NLAT)] = {0};

float plasic_tdissd[(PLASIC_NLEV)] = {0};
float plasic_tdissz[(PLASIC_NLEV)] = {0};
float plasic_tdisst[(PLASIC_NLEV)] = {0};
float plasic_tdissq[(PLASIC_NLEV)] = {0};

float plasic_restim[(PLASIC_NLEV)] = {0};

// Reference (restoration) temperature profile T0(k), in K.
float plasic_t0[(PLASIC_NLEV)] = {0};
// Rayleigh friction (sponge-layer) time scale of each level, in seconds;
// shorter time means stronger damping and is nonzero only near the model top.
float plasic_tfrc[(PLASIC_NLEV)] = {0};
// User-supplied sigma half-level interfaces, used when neqsig selects a custom
// vertical grid.
float plasic_sigmah_configured[(PLASIC_NLEV)] = {0};
// Newtonian temperature restoration (relaxation) rate toward the reference
// profile T0.
float plasic_damp[(PLASIC_NLEV)] = {0};
// Sigma-coordinate layer thickness, dsigma(k) = sigmah(k) - sigmah(k-1).
float plasic_dsigma[(PLASIC_NLEV)] = {0};
// Reciprocal of twice the layer thickness, 1/(2*dsigma) = 0.5/dsigma. Beware:
// this is NOT 1/dsigma; the factor 2 comes from the centered vertical
// difference used by the reference-temperature vertical advection.
float plasic_rdsig[(PLASIC_NLEV)] = {0};
// Sigma-coordinate layer centers, sigma = p/ps, decreasing from near 1 at the
// surface to near 0 at the model top.
float plasic_sigma[(PLASIC_NLEV)] = {0};
// Sigma-coordinate layer interfaces (half levels), sigmah.
float plasic_sigmah[(PLASIC_NLEV)] = {0};
// Reference temperature difference between adjacent layers, T0(k+1) - T0(k);
// the surface-most entry has no layer below it and is set to zero.
float plasic_t01s2[(PLASIC_NLEV)] = {0};
// Reference temperature scaled by the adiabatic exponent, kappa * T0(k).
float plasic_tkp[(PLASIC_NLEV)] = {0};

float plasic_c[(PLASIC_NLEV) * (PLASIC_NLEV)] = {0};
float plasic_g[(PLASIC_NLEV) * (PLASIC_NLEV)] = {0};
float plasic_tau[(PLASIC_NLEV) * (PLASIC_NLEV)] = {0};
float plasic_bm1[(PLASIC_NLEV) * (PLASIC_NLEV) * (PLASIC_NTRU)] = {0};

int32_t plasic_mypid = 0;

// Random number generator seed (8 integers) used to initialise the model's RNG
// for reproducible run; if the first entry is zero a clock-based seed is used.
// Saved and restored with the restart files so that restarts continue the same
// random sequence.s
int32_t plasic_seed[(8)] = {0};
float plasic_kap = 0;
float plasic_ga = 0;
float plasic_gascon = 0;
float plasic_plarad = 0;
// Robert-Asselin time filter coefficient, typically ~0.1
float plasic_pnu = 0;
float plasic_sidereal_day = 0;
float plasic_solar_day = 0;
float plasic_ww = 0;
float plasic_ra1 = 0;
float plasic_ra2 = 0;
float plasic_ra4 = 0;
float plasic_cpd = 0;
float plasic_cpv_cpd_minus1 = 0;
float plasic_cv = 0;
float plasic_ct = 0;
// Companion Robert-Asselin time filter weight, pnu21 = 1 - 2*pnu; the filtered
// state is pnu21*X_new + pnu*X_saved.
float plasic_pnu21 = 0;
float plasic_rdbrv = 0;

static void plasic_fill_int32(int32_t *values, size_t count, int32_t value)
{
    for (size_t index = 0; index < count; ++index) {
        values[index] = value;
    }
}

static void plasic_fill_float(float *values, size_t count, float value)
{
    for (size_t index = 0; index < count; ++index) {
        values[index] = value;
    }
}

void plasic_state_init(void)
{
    static int initialized = 0;

    if (initialized != 0) {
        return;
    }
    initialized = 1;

    plasic_nstep = (int32_t)(0);
    plasic_calendar = (int32_t)(0);
    plasic_mpstep = (float)(0.0);
    plasic_ntspd = (int32_t)(0);
    plasic_fill_int32(plasic_ndatim, sizeof(plasic_ndatim) / sizeof(plasic_ndatim[0]), (int32_t)(-1));
    plasic_nafter = (int32_t)(0);
    plasic_nkits = (int32_t)(3);
    plasic_nrad = (int32_t)(1);
    plasic_nadv = (int32_t)(1);
    plasic_neqsig = (int32_t)(0);
    plasic_naccuout = (int32_t)(0);
    plasic_nhdiff = (int32_t)(15);
    plasic_ndheat = (int32_t)(1);
    plasic_nsponge = (int32_t)(0);
    plasic_ls = (float)(2.8345E6);
    plasic_lv = (float)(2.5008E6);
    plasic_plavor = (float)(PLASIC_Y10_PLANETARY_VORTICITY_FACTOR);
    plasic_deltsec = (float)(0.0);
    plasic_deltsec2 = (float)(0.0);
    plasic_tgr = (float)(288.0);
    plasic_psurf = (float)(101100.0);
    plasic_tmelt = (float)(273.16);
    plasic_dampsp = (float)(0.0);
    plasic_fill_float(plasic_sd, sizeof(plasic_sd) / sizeof(plasic_sd[0]), (float)(0.0));
    plasic_fill_float(plasic_st, sizeof(plasic_st) / sizeof(plasic_st[0]), (float)(0.0));
    plasic_fill_float(plasic_sz, sizeof(plasic_sz) / sizeof(plasic_sz[0]), (float)(0.0));
    plasic_fill_float(plasic_sq, sizeof(plasic_sq) / sizeof(plasic_sq[0]), (float)(0.0));
    plasic_fill_float(plasic_sp, sizeof(plasic_sp) / sizeof(plasic_sp[0]), (float)(0.0));
    plasic_fill_float(plasic_so, sizeof(plasic_so) / sizeof(plasic_so[0]), (float)(0.0));
    plasic_fill_float(plasic_sr, sizeof(plasic_sr) / sizeof(plasic_sr[0]), (float)(0.0));
    plasic_fill_float(plasic_sdp, sizeof(plasic_sdp) / sizeof(plasic_sdp[0]), (float)(0.0));
    plasic_fill_float(plasic_stp, sizeof(plasic_stp) / sizeof(plasic_stp[0]), (float)(0.0));
    plasic_fill_float(plasic_szp, sizeof(plasic_szp) / sizeof(plasic_szp[0]), (float)(0.0));
    plasic_fill_float(plasic_sqp, sizeof(plasic_sqp) / sizeof(plasic_sqp[0]), (float)(0.0));
    plasic_fill_float(plasic_spp, sizeof(plasic_spp) / sizeof(plasic_spp[0]), (float)(0.0));
    plasic_fill_float(plasic_sop, sizeof(plasic_sop) / sizeof(plasic_sop[0]), (float)(0.0));
    plasic_fill_float(plasic_srp, sizeof(plasic_srp) / sizeof(plasic_srp[0]), (float)(0.0));
    plasic_fill_float(plasic_sdt, sizeof(plasic_sdt) / sizeof(plasic_sdt[0]), (float)(0.0));
    plasic_fill_float(plasic_stt, sizeof(plasic_stt) / sizeof(plasic_stt[0]), (float)(0.0));
    plasic_fill_float(plasic_szt, sizeof(plasic_szt) / sizeof(plasic_szt[0]), (float)(0.0));
    plasic_fill_float(plasic_sqt, sizeof(plasic_sqt) / sizeof(plasic_sqt[0]), (float)(0.0));
    plasic_fill_float(plasic_spt, sizeof(plasic_spt) / sizeof(plasic_spt[0]), (float)(0.0));
    plasic_fill_float(plasic_sdm, sizeof(plasic_sdm) / sizeof(plasic_sdm[0]), (float)(0.0));
    plasic_fill_float(plasic_stm, sizeof(plasic_stm) / sizeof(plasic_stm[0]), (float)(0.0));
    plasic_fill_float(plasic_szm, sizeof(plasic_szm) / sizeof(plasic_szm[0]), (float)(0.0));
    plasic_fill_float(plasic_sqm, sizeof(plasic_sqm) / sizeof(plasic_sqm[0]), (float)(0.0));
    plasic_fill_float(plasic_spm, sizeof(plasic_spm) / sizeof(plasic_spm[0]), (float)(0.0));
    plasic_fill_float(plasic_Lnk, sizeof(plasic_Lnk) / sizeof(plasic_Lnk[0]), (float)(0.0));
    plasic_fill_float(plasic_Lnkpp, sizeof(plasic_Lnkpp) / sizeof(plasic_Lnkpp[0]), (float)(0.0));
    plasic_fill_float(plasic_spnorm, sizeof(plasic_spnorm) / sizeof(plasic_spnorm[0]), (float)(0.0));
    plasic_fill_int32(plasic_nindex, sizeof(plasic_nindex) / sizeof(plasic_nindex[0]), (int32_t)(PLASIC_NTRU));
    plasic_fill_int32(plasic_ndel, sizeof(plasic_ndel) / sizeof(plasic_ndel[0]), (int32_t)(2));
    plasic_fill_float(plasic_gd, sizeof(plasic_gd) / sizeof(plasic_gd[0]), (float)(0.));
    plasic_fill_float(plasic_gt, sizeof(plasic_gt) / sizeof(plasic_gt[0]), (float)(0.));
    plasic_fill_float(plasic_gz, sizeof(plasic_gz) / sizeof(plasic_gz[0]), (float)(0.));
    plasic_fill_float(plasic_gq, sizeof(plasic_gq) / sizeof(plasic_gq[0]), (float)(0.));
    plasic_fill_float(plasic_gu, sizeof(plasic_gu) / sizeof(plasic_gu[0]), (float)(0.));
    plasic_fill_float(plasic_gv, sizeof(plasic_gv) / sizeof(plasic_gv[0]), (float)(0.));
    plasic_fill_float(plasic_gtdt, sizeof(plasic_gtdt) / sizeof(plasic_gtdt[0]), (float)(0.));
    plasic_fill_float(plasic_gqdt, sizeof(plasic_gqdt) / sizeof(plasic_gqdt[0]), (float)(0.));
    plasic_fill_float(plasic_gudt, sizeof(plasic_gudt) / sizeof(plasic_gudt[0]), (float)(0.));
    plasic_fill_float(plasic_gvdt, sizeof(plasic_gvdt) / sizeof(plasic_gvdt[0]), (float)(0.));
    plasic_fill_float(plasic_gp, sizeof(plasic_gp) / sizeof(plasic_gp[0]), (float)(0.));
    plasic_fill_float(plasic_gpj, sizeof(plasic_gpj) / sizeof(plasic_gpj[0]), (float)(0.));
    plasic_fill_float(plasic_rcsq, sizeof(plasic_rcsq) / sizeof(plasic_rcsq[0]), (float)(0.));
    plasic_fill_float(plasic_dt, sizeof(plasic_dt) / sizeof(plasic_dt[0]), (float)(0.));
    plasic_fill_float(plasic_dq, sizeof(plasic_dq) / sizeof(plasic_dq[0]), (float)(0.));
    plasic_fill_float(plasic_du, sizeof(plasic_du) / sizeof(plasic_du[0]), (float)(0.));
    plasic_fill_float(plasic_dv, sizeof(plasic_dv) / sizeof(plasic_dv[0]), (float)(0.));
    plasic_fill_float(plasic_dp, sizeof(plasic_dp) / sizeof(plasic_dp[0]), (float)(0.));
    plasic_fill_float(plasic_dqsat, sizeof(plasic_dqsat) / sizeof(plasic_dqsat[0]), (float)(0.));
    plasic_fill_float(plasic_dqt, sizeof(plasic_dqt) / sizeof(plasic_dqt[0]), (float)(0.));
    plasic_fill_float(plasic_dcc, sizeof(plasic_dcc) / sizeof(plasic_dcc[0]), (float)(0.));
    plasic_fill_float(plasic_dql, sizeof(plasic_dql) / sizeof(plasic_dql[0]), (float)(0.));
    plasic_fill_float(plasic_dw, sizeof(plasic_dw) / sizeof(plasic_dw[0]), (float)(0.));
    plasic_fill_float(plasic_dtdt, sizeof(plasic_dtdt) / sizeof(plasic_dtdt[0]), (float)(0.));
    plasic_fill_float(plasic_dqdt, sizeof(plasic_dqdt) / sizeof(plasic_dqdt[0]), (float)(0.));
    plasic_fill_float(plasic_dudt, sizeof(plasic_dudt) / sizeof(plasic_dudt[0]), (float)(0.));
    plasic_fill_float(plasic_dvdt, sizeof(plasic_dvdt) / sizeof(plasic_dvdt[0]), (float)(0.));
    plasic_fill_float(plasic_dp0, sizeof(plasic_dp0) / sizeof(plasic_dp0[0]), (float)(0.));
    plasic_fill_float(plasic_du0, sizeof(plasic_du0) / sizeof(plasic_du0[0]), (float)(0.));
    plasic_fill_float(plasic_dv0, sizeof(plasic_dv0) / sizeof(plasic_dv0[0]), (float)(0.));
    plasic_fill_float(plasic_dwetfac, sizeof(plasic_dwetfac) / sizeof(plasic_dwetfac[0]), (float)(0.));
    plasic_fill_float(plasic_dls, sizeof(plasic_dls) / sizeof(plasic_dls[0]), (float)(1.));
    plasic_fill_float(plasic_dz0, sizeof(plasic_dz0) / sizeof(plasic_dz0[0]), (float)(0.));
    plasic_fill_float(plasic_diced, sizeof(plasic_diced) / sizeof(plasic_diced[0]), (float)(0.));
    plasic_fill_float(plasic_dicec, sizeof(plasic_dicec) / sizeof(plasic_dicec[0]), (float)(0.));
    plasic_fill_float(plasic_dtaux, sizeof(plasic_dtaux) / sizeof(plasic_dtaux[0]), (float)(0.));
    plasic_fill_float(plasic_dtauy, sizeof(plasic_dtauy) / sizeof(plasic_dtauy[0]), (float)(0.));
    plasic_fill_float(plasic_dust3, sizeof(plasic_dust3) / sizeof(plasic_dust3[0]), (float)(0.));
    plasic_fill_float(plasic_dshfl, sizeof(plasic_dshfl) / sizeof(plasic_dshfl[0]), (float)(0.));
    plasic_fill_float(plasic_dlhfl, sizeof(plasic_dlhfl) / sizeof(plasic_dlhfl[0]), (float)(0.));
    plasic_fill_float(plasic_devap, sizeof(plasic_devap) / sizeof(plasic_devap[0]), (float)(0.));
    plasic_fill_float(plasic_dtsa, sizeof(plasic_dtsa) / sizeof(plasic_dtsa[0]), (float)(0.));
    plasic_fill_float(plasic_dmld, sizeof(plasic_dmld) / sizeof(plasic_dmld[0]), (float)(0.));
    plasic_fill_float(plasic_dwmax, sizeof(plasic_dwmax) / sizeof(plasic_dwmax[0]), (float)(0.0));
    plasic_fill_float(plasic_dwatc, sizeof(plasic_dwatc) / sizeof(plasic_dwatc[0]), (float)(0.));
    plasic_fill_float(plasic_dsnow, sizeof(plasic_dsnow) / sizeof(plasic_dsnow[0]), (float)(0.));
    plasic_fill_float(plasic_dsmelt, sizeof(plasic_dsmelt) / sizeof(plasic_dsmelt[0]), (float)(0.));
    plasic_fill_float(plasic_dsndch, sizeof(plasic_dsndch) / sizeof(plasic_dsndch[0]), (float)(0.));
    plasic_fill_float(plasic_dtsoil, sizeof(plasic_dtsoil) / sizeof(plasic_dtsoil[0]), (float)(0.));
    plasic_fill_float(plasic_dtd2, sizeof(plasic_dtd2) / sizeof(plasic_dtd2[0]), (float)(0.));
    plasic_fill_float(plasic_dtd3, sizeof(plasic_dtd3) / sizeof(plasic_dtd3[0]), (float)(0.));
    plasic_fill_float(plasic_dtd4, sizeof(plasic_dtd4) / sizeof(plasic_dtd4[0]), (float)(0.));
    plasic_fill_float(plasic_dtd5, sizeof(plasic_dtd5) / sizeof(plasic_dtd5[0]), (float)(0.));
    plasic_fill_float(plasic_dglac, sizeof(plasic_dglac) / sizeof(plasic_dglac[0]), (float)(0.));
    plasic_fill_float(plasic_aevap, sizeof(plasic_aevap) / sizeof(plasic_aevap[0]), (float)(0.));
    plasic_fill_float(plasic_aprl, sizeof(plasic_aprl) / sizeof(plasic_aprl[0]), (float)(0.));
    plasic_fill_float(plasic_aprc, sizeof(plasic_aprc) / sizeof(plasic_aprc[0]), (float)(0.));
    plasic_fill_float(plasic_aprs, sizeof(plasic_aprs) / sizeof(plasic_aprs[0]), (float)(0.));
    plasic_fill_float(plasic_ashfl, sizeof(plasic_ashfl) / sizeof(plasic_ashfl[0]), (float)(0.));
    plasic_fill_float(plasic_alhfl, sizeof(plasic_alhfl) / sizeof(plasic_alhfl[0]), (float)(0.));
    plasic_fill_float(plasic_asmelt, sizeof(plasic_asmelt) / sizeof(plasic_asmelt[0]), (float)(0.));
    plasic_fill_float(plasic_asndch, sizeof(plasic_asndch) / sizeof(plasic_asndch[0]), (float)(0.));
    plasic_fill_float(plasic_acc, sizeof(plasic_acc) / sizeof(plasic_acc[0]), (float)(0.));
    plasic_fill_float(plasic_assol, sizeof(plasic_assol) / sizeof(plasic_assol[0]), (float)(0.));
    plasic_fill_float(plasic_asthr, sizeof(plasic_asthr) / sizeof(plasic_asthr[0]), (float)(0.));
    plasic_fill_float(plasic_atsol, sizeof(plasic_atsol) / sizeof(plasic_atsol[0]), (float)(0.));
    plasic_fill_float(plasic_atthr, sizeof(plasic_atthr) / sizeof(plasic_atthr[0]), (float)(0.));
    plasic_fill_float(plasic_assolu, sizeof(plasic_assolu) / sizeof(plasic_assolu[0]), (float)(0.));
    plasic_fill_float(plasic_asthru, sizeof(plasic_asthru) / sizeof(plasic_asthru[0]), (float)(0.));
    plasic_fill_float(plasic_atsolu, sizeof(plasic_atsolu) / sizeof(plasic_atsolu[0]), (float)(0.));
    plasic_fill_float(plasic_ataux, sizeof(plasic_ataux) / sizeof(plasic_ataux[0]), (float)(0.));
    plasic_fill_float(plasic_atauy, sizeof(plasic_atauy) / sizeof(plasic_atauy[0]), (float)(0.));
    plasic_fill_float(plasic_aqvi, sizeof(plasic_aqvi) / sizeof(plasic_aqvi[0]), (float)(0.));
    plasic_fill_float(plasic_atsa, sizeof(plasic_atsa) / sizeof(plasic_atsa[0]), (float)(0.));
    plasic_fill_float(plasic_atsama, sizeof(plasic_atsama) / sizeof(plasic_atsama[0]), (float)(0.));
    plasic_fill_float(plasic_atsami, sizeof(plasic_atsami) / sizeof(plasic_atsami[0]), (float)(0.));
    plasic_fill_float(plasic_ats0, sizeof(plasic_ats0) / sizeof(plasic_ats0[0]), (float)(0.));
    plasic_fill_float(plasic_tdissd, sizeof(plasic_tdissd) / sizeof(plasic_tdissd[0]), (float)(0.20));
    plasic_fill_float(plasic_tdissz, sizeof(plasic_tdissz) / sizeof(plasic_tdissz[0]), (float)(1.10));
    plasic_fill_float(plasic_tdisst, sizeof(plasic_tdisst) / sizeof(plasic_tdisst[0]), (float)(5.60));
    plasic_fill_float(plasic_tdissq, sizeof(plasic_tdissq) / sizeof(plasic_tdissq[0]), (float)(0.1));
    plasic_fill_float(plasic_restim, sizeof(plasic_restim) / sizeof(plasic_restim[0]), (float)(0.0));
    plasic_fill_float(plasic_t0, sizeof(plasic_t0) / sizeof(plasic_t0[0]), (float)(250.0));
    plasic_fill_float(plasic_tfrc, sizeof(plasic_tfrc) / sizeof(plasic_tfrc[0]), (float)(0.0));
    plasic_fill_float(plasic_sigmah_configured, sizeof(plasic_sigmah_configured) / sizeof(plasic_sigmah_configured[0]), (float)(0.0));
    plasic_mypid = (int32_t)(0);
    plasic_fill_int32(plasic_seed, sizeof(plasic_seed) / sizeof(plasic_seed[0]), (int32_t)(0));
    plasic_kap = (float)(0.0);
    plasic_ga = (float)(0.0);
    plasic_gascon = (float)(0.0);
    plasic_plarad = (float)(0.0);
    plasic_pnu = (float)(0.0);
    plasic_sidereal_day = (float)(0.0);
    plasic_solar_day = (float)(0.0);
    plasic_ww = (float)(0.0);
    plasic_ra1 = (float)(0.0);
    plasic_ra2 = (float)(0.0);
    plasic_ra4 = (float)(0.0);
    plasic_cpd = (float)(0.0);
    plasic_cpv_cpd_minus1 = (float)(0.0);
    plasic_cv = (float)(0.0);
    plasic_ct = (float)(0.0);
    plasic_pnu21 = (float)(0.0);
    plasic_rdbrv = (float)(0.0);
}
