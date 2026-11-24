#ifndef PLASICMOD_STATE_H
#define PLASICMOD_STATE_H

#include <stdint.h>

#if !defined(PLASIC_NLAT) || !defined(PLASIC_NLEV) || !defined(PLASIC_NPRO)
#error "Define PLASIC_NLAT, PLASIC_NLEV, and PLASIC_NPRO"
#endif

/*
 * 目前支持的水平分辨率及其默认时间步长（秒）；其他 NLAT 直接编译失败。
 * 这些默认步长足够短，保证长时间气候积分稳定、不发生数值发散。用户用
 * --step-minutes（JSON 配置中的 step_minutes）给出的步长会直接替代默认值，
 * 程序不检验其稳定性。
 *
 * Supported horizontal resolutions and their default time step in seconds;
 * any other NLAT fails to compile. The defaults are short enough to keep long
 * climate integrations stable. A step given with --step-minutes (step_minutes
 * in the JSON config) replaces the default and is not checked for stability.
 */
#if PLASIC_NLAT == 32 /* T21 */
#define PLASIC_DEFAULT_STEP_SECONDS 1800 /* 30 min, 48 steps per day */
#elif PLASIC_NLAT == 48 /* T31 */
#define PLASIC_DEFAULT_STEP_SECONDS 1080 /* 18 min, 80 steps per day */
#elif PLASIC_NLAT == 64 /* T42 */
#define PLASIC_DEFAULT_STEP_SECONDS 720 /* 12 min, 120 steps per day */
#elif PLASIC_NLAT == 128 /* T85 */
#define PLASIC_DEFAULT_STEP_SECONDS 450 /* 7.5 min, 192 steps per day */
#elif PLASIC_NLAT == 256 /* T170 */
#define PLASIC_DEFAULT_STEP_SECONDS 150 /* 2.5 min, 576 steps per day */
#else
#error "Unsupported resolution: PLASIC_NLAT must be 32 (T21), 48 (T31), 64 (T42), 128 (T85) or 256 (T170)"
#endif

#define PLASIC_NLON (PLASIC_NLAT + PLASIC_NLAT)
#define PLASIC_NTRU ((PLASIC_NLON - 1) / 3)
#define PLASIC_NLPP (PLASIC_NLAT / PLASIC_NPRO)
#define PLASIC_NHOR (PLASIC_NLON * PLASIC_NLPP)
#define PLASIC_NUGP (PLASIC_NLON * PLASIC_NLAT)
#define PLASIC_NLEM (PLASIC_NLEV - 1)
#define PLASIC_NLEP (PLASIC_NLEV + 1)
#define PLASIC_NLSQ (PLASIC_NLEV * PLASIC_NLEV)
#define PLASIC_NTP1 (PLASIC_NTRU + 1)
#define PLASIC_NRSP ((PLASIC_NTRU + 1) * (PLASIC_NTRU + 2))
#define PLASIC_NSPP ((PLASIC_NRSP + PLASIC_NPRO - 1) / PLASIC_NPRO)
/*
 * 每个 rank 的谱块按 NSPP 个槽位预留。NPRO 整除 NRSP 时 NESP 恰好等于
 * NRSP；否则最后一个 rank 的尾部槽位只是对齐用的 padding，既不参与
 * 谱通信，也不含任何物理意义。
 *
 * Every rank reserves NSPP spectral slots. When NPRO divides NRSP, NESP is
 * exactly NRSP; otherwise the trailing slots of the last rank are alignment
 * padding that neither participates in spectral communication nor carries
 * any physical meaning.
 */
#define PLASIC_NESP ((PLASIC_NSPP) * (PLASIC_NPRO))

#if PLASIC_USE_MPI && (PLASIC_NLAT % PLASIC_NPRO) != 0
#error "MPI builds require PLASIC_NPRO to divide PLASIC_NLAT exactly"
#endif

/*
 * 当前 rank 实际拥有的有效谱槽数：通常等于 PLASIC_NSPP，只有最后一个
 * rank 在 PLASIC_NRSP 不能被 PLASIC_NPRO 整除时少于 PLASIC_NSPP。
 *
 * Number of valid spectral slots owned by this rank: normally PLASIC_NSPP,
 * and only the last rank owns fewer when PLASIC_NPRO does not divide
 * PLASIC_NRSP.
 */
static inline int32_t plasic_local_spectral_count(int32_t rank)
{
    const int32_t remaining = PLASIC_NRSP - rank * PLASIC_NSPP;

    if (remaining < PLASIC_NSPP)
    {
        return remaining > 0 ? remaining : 0;
    }
    return PLASIC_NSPP;
}

/* Converts the normalized Y_(1,0) mode into 2 sin(latitude). */
#define PLASIC_Y10_PLANETARY_VORTICITY_FACTOR 1.63299310207


extern int32_t plasic_nstep;
extern int32_t plasic_calendar;
extern float plasic_mpstep;
extern int32_t plasic_ntspd;
extern int32_t plasic_ndatim[(7)];
extern int32_t plasic_nafter;
extern int32_t plasic_nkits;
extern int32_t plasic_nrad;
extern int32_t plasic_nadv;
extern int32_t plasic_neqsig;
extern int32_t plasic_naccuout;
extern int32_t plasic_nhdiff;
extern int32_t plasic_ndheat;
extern int32_t plasic_nsponge;
extern float plasic_ls;
extern float plasic_lv;
extern float plasic_plavor;
extern float plasic_deltsec;
extern float plasic_deltsec2;
extern float plasic_delt;
extern float plasic_delt2;
extern float plasic_tgr;
extern float plasic_psurf;
extern float plasic_tmelt;
extern float plasic_dampsp;
extern float plasic_sd[(PLASIC_NRSP) * (PLASIC_NLEV)];
extern float plasic_st[(PLASIC_NRSP) * (PLASIC_NLEV)];
extern float plasic_sz[(PLASIC_NRSP) * (PLASIC_NLEV)];
extern float plasic_sq[(PLASIC_NRSP) * (PLASIC_NLEV)];
extern float plasic_sp[(PLASIC_NRSP)];
extern float plasic_so[(PLASIC_NRSP)];
extern float plasic_sr[(PLASIC_NRSP) * (PLASIC_NLEV)];
extern float plasic_sdp[(PLASIC_NSPP) * (PLASIC_NLEV)];
extern float plasic_stp[(PLASIC_NSPP) * (PLASIC_NLEV)];
extern float plasic_szp[(PLASIC_NSPP) * (PLASIC_NLEV)];
extern float plasic_sqp[(PLASIC_NSPP) * (PLASIC_NLEV)];
extern float plasic_spp[(PLASIC_NSPP)];
extern float plasic_sop[(PLASIC_NSPP)];
extern float plasic_srp[(PLASIC_NSPP) * (PLASIC_NLEV)];
extern float plasic_sdt[(PLASIC_NSPP) * (PLASIC_NLEV)];
extern float plasic_stt[(PLASIC_NSPP) * (PLASIC_NLEV)];
extern float plasic_szt[(PLASIC_NSPP) * (PLASIC_NLEV)];
extern float plasic_sqt[(PLASIC_NSPP) * (PLASIC_NLEV)];
extern float plasic_spt[(PLASIC_NSPP)];
extern float plasic_sdm[(PLASIC_NSPP) * (PLASIC_NLEV)];
extern float plasic_stm[(PLASIC_NSPP) * (PLASIC_NLEV)];
extern float plasic_szm[(PLASIC_NSPP) * (PLASIC_NLEV)];
extern float plasic_sqm[(PLASIC_NSPP) * (PLASIC_NLEV)];
extern float plasic_spm[(PLASIC_NSPP)];
extern float plasic_Lnk[(PLASIC_NRSP) * (PLASIC_NLEV)];
extern float plasic_Lnkpp[(PLASIC_NSPP) * (PLASIC_NLEV)];
extern float plasic_spnorm[(PLASIC_NRSP)];
extern int32_t plasic_nindex[(PLASIC_NESP)];
extern int32_t plasic_ndel[(PLASIC_NLEV)];
extern float plasic_gd[(PLASIC_NHOR) * (PLASIC_NLEV)];
extern float plasic_gt[(PLASIC_NHOR) * (PLASIC_NLEV)];
extern float plasic_gz[(PLASIC_NHOR) * (PLASIC_NLEV)];
extern float plasic_gq[(PLASIC_NHOR) * (PLASIC_NLEV)];
extern float plasic_gu[(PLASIC_NHOR) * (PLASIC_NLEV)];
extern float plasic_gv[(PLASIC_NHOR) * (PLASIC_NLEV)];
extern float plasic_gtdt[(PLASIC_NHOR) * (PLASIC_NLEV)];
extern float plasic_gqdt[(PLASIC_NHOR) * (PLASIC_NLEV)];
extern float plasic_gudt[(PLASIC_NHOR) * (PLASIC_NLEV)];
extern float plasic_gvdt[(PLASIC_NHOR) * (PLASIC_NLEV)];
extern float plasic_gp[(PLASIC_NHOR)];
extern float plasic_gpj[(PLASIC_NHOR)];
extern float plasic_rcsq[(PLASIC_NHOR)];
extern float plasic_dt[(PLASIC_NHOR) * (PLASIC_NLEP)];
extern float plasic_dq[(PLASIC_NHOR) * (PLASIC_NLEP)];
extern float plasic_du[(PLASIC_NHOR) * (PLASIC_NLEP)];
extern float plasic_dv[(PLASIC_NHOR) * (PLASIC_NLEP)];
extern float plasic_dp[(PLASIC_NHOR)];
extern float plasic_dqsat[(PLASIC_NHOR) * (PLASIC_NLEP)];
extern float plasic_dqt[(PLASIC_NHOR) * (PLASIC_NLEP)];
extern float plasic_dcc[(PLASIC_NHOR) * (PLASIC_NLEP)];
extern float plasic_dql[(PLASIC_NHOR) * (PLASIC_NLEP)];
extern float plasic_dw[(PLASIC_NHOR) * (PLASIC_NLEV)];
extern float plasic_dtdt[(PLASIC_NHOR) * (PLASIC_NLEP)];
extern float plasic_dqdt[(PLASIC_NHOR) * (PLASIC_NLEP)];
extern float plasic_dudt[(PLASIC_NHOR) * (PLASIC_NLEP)];
extern float plasic_dvdt[(PLASIC_NHOR) * (PLASIC_NLEP)];
extern float plasic_dp0[(PLASIC_NHOR)];
extern float plasic_du0[(PLASIC_NHOR) * (PLASIC_NLEP)];
extern float plasic_dv0[(PLASIC_NHOR) * (PLASIC_NLEP)];
extern float plasic_dalb[(PLASIC_NHOR)];
extern float plasic_dswfl[(PLASIC_NHOR) * (PLASIC_NLEP)];
extern float plasic_dlwfl[(PLASIC_NHOR) * (PLASIC_NLEP)];
extern float plasic_dflux[(PLASIC_NHOR) * (PLASIC_NLEP)];
extern float plasic_dfu[(PLASIC_NHOR) * (PLASIC_NLEP)];
extern float plasic_dfd[(PLASIC_NHOR) * (PLASIC_NLEP)];
extern float plasic_dftu[(PLASIC_NHOR) * (PLASIC_NLEP)];
extern float plasic_dftd[(PLASIC_NHOR) * (PLASIC_NLEP)];
extern float plasic_dwetfac[(PLASIC_NHOR)];
extern float plasic_dls[(PLASIC_NHOR)];
extern float plasic_dz0[(PLASIC_NHOR)];
extern float plasic_diced[(PLASIC_NHOR)];
extern float plasic_dicec[(PLASIC_NHOR)];
extern float plasic_dtaux[(PLASIC_NHOR)];
extern float plasic_dtauy[(PLASIC_NHOR)];
extern float plasic_dust3[(PLASIC_NHOR)];
extern float plasic_dshfl[(PLASIC_NHOR)];
extern float plasic_dlhfl[(PLASIC_NHOR)];
extern float plasic_devap[(PLASIC_NHOR)];
extern float plasic_dtsa[(PLASIC_NHOR)];
extern float plasic_dmld[(PLASIC_NHOR)];
extern float plasic_dshdt[(PLASIC_NHOR)];
extern float plasic_dlhdt[(PLASIC_NHOR)];
extern float plasic_dprc[(PLASIC_NHOR)];
extern float plasic_dprl[(PLASIC_NHOR)];
extern float plasic_dprs[(PLASIC_NHOR)];
extern float plasic_dqvi[(PLASIC_NHOR)];
extern float plasic_dwmax[(PLASIC_NHOR)];
extern float plasic_dwatc[(PLASIC_NHOR)];
extern float plasic_dsnow[(PLASIC_NHOR)];
extern float plasic_dsmelt[(PLASIC_NHOR)];
extern float plasic_dsndch[(PLASIC_NHOR)];
extern float plasic_dtsoil[(PLASIC_NHOR)];
extern float plasic_dtd2[(PLASIC_NHOR)];
extern float plasic_dtd3[(PLASIC_NHOR)];
extern float plasic_dtd4[(PLASIC_NHOR)];
extern float plasic_dtd5[(PLASIC_NHOR)];
extern float plasic_dglac[(PLASIC_NHOR)];
extern float plasic_aevap[(PLASIC_NHOR)];
extern float plasic_aprl[(PLASIC_NHOR)];
extern float plasic_aprc[(PLASIC_NHOR)];
extern float plasic_aprs[(PLASIC_NHOR)];
extern float plasic_ashfl[(PLASIC_NHOR)];
extern float plasic_alhfl[(PLASIC_NHOR)];
extern float plasic_asmelt[(PLASIC_NHOR)];
extern float plasic_asndch[(PLASIC_NHOR)];
extern float plasic_acc[(PLASIC_NHOR)];
extern float plasic_assol[(PLASIC_NHOR)];
extern float plasic_asthr[(PLASIC_NHOR)];
extern float plasic_atsol[(PLASIC_NHOR)];
extern float plasic_atthr[(PLASIC_NHOR)];
extern float plasic_assolu[(PLASIC_NHOR)];
extern float plasic_asthru[(PLASIC_NHOR)];
extern float plasic_atsolu[(PLASIC_NHOR)];
extern float plasic_ataux[(PLASIC_NHOR)];
extern float plasic_atauy[(PLASIC_NHOR)];
extern float plasic_aqvi[(PLASIC_NHOR)];
extern float plasic_atsa[(PLASIC_NHOR)];
extern float plasic_atsama[(PLASIC_NHOR)];
extern float plasic_atsami[(PLASIC_NHOR)];
extern float plasic_ats0[(PLASIC_NHOR)];
extern double plasic_sid[(PLASIC_NLAT)];
extern double plasic_gwd[(PLASIC_NLAT)];
extern float plasic_csq[(PLASIC_NLAT)];
extern float plasic_cola[(PLASIC_NLAT)];
extern float plasic_rcs[(PLASIC_NLAT)];
extern float plasic_tdissd[(PLASIC_NLEV)];
extern float plasic_tdissz[(PLASIC_NLEV)];
extern float plasic_tdisst[(PLASIC_NLEV)];
extern float plasic_tdissq[(PLASIC_NLEV)];
extern float plasic_restim[(PLASIC_NLEV)];
extern float plasic_t0[(PLASIC_NLEV)];
extern float plasic_tfrc[(PLASIC_NLEV)];
extern float plasic_sigmah_configured[(PLASIC_NLEV)];
extern float plasic_damp[(PLASIC_NLEV)];
extern float plasic_dsigma[(PLASIC_NLEV)];
extern float plasic_rdsig[(PLASIC_NLEV)];
extern float plasic_sigma[(PLASIC_NLEV)];
extern float plasic_sigmah[(PLASIC_NLEV)];
extern float plasic_t01s2[(PLASIC_NLEV)];
extern float plasic_tkp[(PLASIC_NLEV)];
extern float plasic_c[(PLASIC_NLEV) * (PLASIC_NLEV)];
extern float plasic_g[(PLASIC_NLEV) * (PLASIC_NLEV)];
extern float plasic_tau[(PLASIC_NLEV) * (PLASIC_NLEV)];
extern float plasic_bm1[(PLASIC_NLEV) * (PLASIC_NLEV) * (PLASIC_NTRU)];
extern int32_t plasic_mypid;
extern int32_t plasic_seed[(8)];
extern float plasic_kap;
extern float plasic_ga;
extern float plasic_gascon;
extern float plasic_plarad;
extern float plasic_pnu;
extern float plasic_sidereal_day;
extern float plasic_solar_day;
extern float plasic_ww;
extern float plasic_ra1;
extern float plasic_ra2;
extern float plasic_ra4;
extern float plasic_cpd;
extern float plasic_cpv_cpd_minus1;
extern float plasic_cv;
extern float plasic_ct;
extern float plasic_pnu21;
extern float plasic_rdbrv;

extern const float plasic_two_pi;
extern const float plasic_water_gas_constant;
extern const float plasic_water_heat_capacity;

void plasic_state_init(void);

#endif
