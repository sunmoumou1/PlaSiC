#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "math_tools/sht_kernels.h"
#include "plasicmod_state.h"

static const double k_plarad = 6371220.0;
static const double k_sidereal_day = 86164.0916;
static const double k_gascon = 287.0;
static const double k_tgr = 288.0;
static const double k_t0 = 250.0;
static const double k_psurf0 = 101100.0;
static const double k_two_pi = 6.283185307179586476925286766559;
static double planet_cv(void)
{
    return k_plarad * (k_two_pi / k_sidereal_day);
}
/*
 * Compute the temperature scale ct = cv^2/Rd. The dimensional combination of cv^2 and Rd yields kelvin.
 */
static double planet_ct(void)
{
    const double cv = planet_cv();
    return cv * cv / k_gascon;
}

/*
 * Read exactly count float values from a binary file. Return 0 on success and 1 if the file cannot be opened or contains fewer values than expected. The data are read using the host float representation; no byte swapping is performed.
 */
static int read_bin(const char *path, float *values, size_t count)
{
    /*
     * 以二进制只读方式打开文件；values 必须已经指向至少 count*sizeof(float) 字节的可写缓冲区。
     * Open the file in binary read mode; values must already point to a writable buffer of at least count*sizeof(float) bytes.
     */
    FILE *stream = fopen(path, "rb");
    size_t got;
    if (stream == NULL)
    {
        fprintf(stderr, "cannot open '%s'\n", path);
        return 1;
    }
    /*
     * fread 的返回值是成功读取的元素个数，而不是字节数；这里要求必须恰好读取 count 个 float。
     * fread returns the number of successfully read elements, not bytes; this code requires exactly count float elements.
     */
    got = fread(values, sizeof(float), count, stream);
    fclose(stream);
    if (got != count)
    {
        fprintf(stderr, "'%s' holds %zu values, expected %zu\n", path, got, count);
        // %zu 是 C 语言 printf/fprintf 中专门用于输出 size_t 类型的 无符号十进制整数 的格式说明符。
        return 1;
    }
    return 0;
}

/*
 * 将 count 个 float 原样写入二进制文件。成功返回 0；创建文件失败或发生 short write 返回 1。
 * Write count float values verbatim to a binary file. Return 0 on success and 1 if file creation fails or a short write occurs.
 */
static int write_bin(const char *path, const float *values, size_t count)
{
    FILE *stream = fopen(path, "wb");
    if (stream == NULL)
    {
        fprintf(stderr, "cannot create '%s'\n", path);
        return 1;
    }
    /*
     * 检查 fwrite 是否写满全部元素；若未写满，先关闭文件再报告错误。
     * Verify that fwrite writes all requested elements; on a short write, close the file before reporting failure.
     */
    if (fwrite(values, sizeof(float), count, stream) != count)
    {
        fclose(stream);
        fprintf(stderr, "short write to '%s'\n", path);
        return 1;
    }
    fclose(stream);
    return 0;
}

/*
 * 读取 PlaSiC 的 sequential-unformatted restart 记录。每个变量先有一个名称记录 [u32=16][16-byte name][u32=16]，随后是数据记录 [u32 payload_bytes][payload][u32 payload_bytes]。代码通过前后 record marker 检查记录边界，并按变量名选择目标记录。
 * English: Read PlaSiC sequential-unformatted restart records. Each variable has a name record [u32=16][16-byte name][u32=16], followed by a data record [u32 payload_bytes][payload][u32 payload_bytes]. The code validates the leading/trailing record markers and selects the requested record by name.
 */
static int restart_read_named(const char *path, const char *wanted, float *out,
                              size_t expected_count)
{
    /*
     * found 用作目标记录是否已找到的标志；函数最后据此区分正常找到与扫描到文件末尾但未找到。
     * found records whether the requested variable has been located; it lets the function distinguish success from reaching EOF without finding the target.
     */
    FILE *stream = fopen(path, "rb");
    int found = 0;
    if (stream == NULL)
    {
        fprintf(stderr, "cannot open restart '%s'\n", path);
        return 1;
    }

    /*
     * 检查 restart 文件中的每个记录，直到找到 wanted 或无法再读出完整的 24 字节名称记录。
     * Scan the restart file record by record until wanted is found or a complete 24-byte name record can no longer be read.
     */
    for (;;)
    {
        /*
         * name_header 保存完整的 Fortran 名称记录：4 字节前导 marker + 16 字节名称 + 4 字节尾随 marker。payload_bytes 保存后续数据记录的字节数。
         * name_header stores the complete Fortran name record: a 4-byte leading marker, 16-byte name, and 4-byte trailing marker. payload_bytes stores the byte length of the following data record.
         */
        unsigned char name_header[24];
        uint32_t marker, payload_bytes, trailing;
        char name[17];
        // 为什么是 17 字节？后面有解释
        if (fread(name_header, 1, sizeof(name_header), stream) !=
            sizeof(name_header))
        {
            break;
        }
        /*
         * 用 memcpy 取出 marker，避免直接把未对齐的 unsigned char* 强制转换为 uint32_t*。这仍然假设 restart 文件的整数字节序与当前主机一致。
         * Use memcpy to extract the markers, avoiding an unaligned cast from unsigned char* to uint32_t*. This still assumes the restart integer byte order matches the host.
         */
        memcpy(&marker, name_header, 4);
        memcpy(&trailing, name_header + 20, 4);
        /*
         * 名称字段固定为 16 字节，因此其前后 Fortran record marker 都应等于 16；不满足时认为文件格式损坏或不兼容。
         * The name field is fixed at 16 bytes, so both Fortran record markers must equal 16; otherwise the file is treated as corrupt or incompatible.
         */
        if (marker != 16 || trailing != 16)
        {
            fprintf(stderr, "bad name record in '%s'\n", path);
            fclose(stream);
            return 1;
        }
        /*
         * 复制固定宽度 16 字节名称，并在第 17 个字节手动补 '\0'，使其可作为 C 字符串处理。
         * Copy the fixed-width 16-byte name and explicitly append '\0' in byte 17 so it can be handled as a C string.
         */
        memcpy(name, name_header + 4, 16);
        name[16] = '\0';

        /*
         * 名称记录之后紧接数据记录的前导 marker；该值就是 payload 的字节数。
         * Immediately after the name record comes the leading marker of the data record; its value is the payload size in bytes.
         */
        if (fread(&payload_bytes, 1, 4, stream) != 4)
        {
            fprintf(stderr, "truncated data marker in '%s'\n", path);
            fclose(stream);
            return 1;
        }

        /*
         * Fortran 固定宽度字符串通常用空格填充。这里构造临时字符串 trimmed，并从尾部删除填充空格后再与 wanted 比较。
         * Fortran fixed-width strings are commonly space padded. A temporary string, trimmed, is built and trailing spaces are removed before comparison with wanted.
         */
        {
            char trimmed[17];
            int i;
            memcpy(trimmed, name, 16);
            trimmed[16] = '\0';
            for (i = 15; i >= 0 && trimmed[i] == ' '; --i)
            {
                trimmed[i] = '\0';
            }
            if (strcmp(trimmed, wanted) == 0)
            {
                /*
                 * 若名称匹配，先验证 payload 大小是否正好等于 expected_count 个 float；这可避免把错误尺寸的数据读入 out。
                 * If the name matches, first verify that the payload is exactly expected_count floats; this prevents reading an incorrectly sized record into out.
                 */
                if (payload_bytes != expected_count * sizeof(float))
                {
                    fprintf(stderr,
                            "restart record '%s' has %u bytes, expected %zu\n",
                            wanted, payload_bytes,
                            expected_count * sizeof(float));
                    fclose(stream);
                    return 1;
                }

                /*
                 * 以字节为单位读取整个 payload，因为 payload_bytes 已经是字节数。
                 * Read the entire payload in bytes because payload_bytes is already expressed in bytes.
                 */
                if (fread(out, 1, payload_bytes, stream) != payload_bytes)
                {
                    fclose(stream);
                    return 1;
                }

                /*
                 * 数据记录尾随 marker 必须与前导 payload_bytes 相等，这是 Fortran sequential-unformatted 记录完整性的基本检查。
                 * The trailing data-record marker must equal the leading payload_bytes value; this is the basic integrity check for a Fortran sequential-unformatted record.
                 */
                if (fread(&trailing, 1, 4, stream) != 4 ||
                    trailing != payload_bytes)
                {
                    fclose(stream);
                    return 1;
                }
                found = 1;
                break;
            }
        }
        /*
         * 当前记录不是目标变量时，跳过 payload 本体以及 4 字节尾随 marker，继续读取下一个名称记录。
         * If the current record is not the requested variable, skip the payload plus its 4-byte trailing marker and continue with the next name record.
         */
        if (fseek(stream, (long)payload_bytes + 4, SEEK_CUR) != 0)
        {
            fclose(stream);
            return 1;
        }
    }
    fclose(stream);
    if (!found)
    {
        fprintf(stderr, "restart '%s' has no record '%s'\n", path, wanted);
        return 1;
    }
    return 0;
}

/*
 * 安全拼接输出目录与文件名
 * Safely join an output directory and filename
 */
static int make_out_path(char *buffer, size_t size, const char *dir, const char *name)
{
    const int n = snprintf(buffer, size, "%s/%s", dir, name);
    return n < 0 || (size_t)n >= size;
}

/*
 * 在相对涡度谱 sz 中加入行星涡度，使其成为模式需要的绝对涡度表示。该修正对每个垂直层执行一次。
 * Add planetary vorticity to the relative-vorticity spectrum sz so that it matches the absolute-vorticity representation required by the model. The correction is applied once per vertical level.
 */
static void add_planetary_vorticity_to_mode_zero_one(float *sz)
{
    const int32_t level_count = PLASIC_NLEV;
    for (int32_t level = 0; level < level_count; ++level)
    {
        /*
         * PlaSiC 的谱数组使用三角截断并按实部/虚部交错存储。对每个垂直层，(m=0,n=1) 的实部位于该层谱块的索引 2；这里向该系数加入行星涡度对应的 Y10 分量。
         * PlaSiC spectra use triangular truncation with interleaved real/imaginary components. For each vertical level, the real part of mode (m=0,n=1) is at index 2 of that level's spectral block; the planetary-vorticity Y10 contribution is added to that coefficient here.
         */
        sz[2 + (size_t)level * (size_t)PLASIC_NRSP] +=
            (float)PLASIC_Y10_PLANETARY_VORTICITY_FACTOR;
    }
}
/*
 * info 子命令只打印编译时确定的网格/谱维度以及 cv、ct，不读写外部数据。
 * The info subcommand only prints compile-time grid/spectral dimensions and the derived cv and ct constants; it does not read or write external data.
 */
static int command_info(void)
{
    printf("NLAT            = %d\n", PLASIC_NLAT);
    printf("NLON            = %d\n", PLASIC_NLON);
    printf("NLEV            = %d\n", PLASIC_NLEV);
    printf("NPRO            = %d\n", PLASIC_NPRO);
    printf("NTRU            = %d\n", PLASIC_NTRU);
    printf("NRSP            = %d\n", PLASIC_NRSP);
    printf("NSPP            = %d\n", PLASIC_NSPP);
    printf("NUGP            = %d\n", PLASIC_NUGP);
    printf("NHOR            = %d\n", PLASIC_NHOR);
    printf("cv              = %.10f\n", planet_cv());
    printf("ct              = %.10f\n", planet_ct());
    return 0;
}

/*
 * orography 子命令：从 restart 中读取谱空间地形 so，将其综合到高斯网格并恢复物理量纲，同时计算与模式运行时一致的参考表面气压。
 * The orography subcommand reads spectral orography so from a restart file, synthesizes it to the Gaussian grid, restores physical dimensions, and computes the same reference surface pressure used by the model runtime.
 */
static int command_orography(const char *restart_path, const char *out_grid,
                             const char *out_meta)
{
    /*
     * so 的长度为每层谱系数总数 PLASIC_NRSP；grid 的长度为全球高斯网格点数 PLASIC_NUGP。so 使用 static 存储，避免把较大的谱数组放在栈上。
     * so has PLASIC_NRSP spectral coefficients for one field; grid has PLASIC_NUGP global Gaussian-grid points. so uses static storage to avoid placing the spectral array on the stack.
     */
    static float so[PLASIC_NRSP];
    float *grid = malloc((size_t)PLASIC_NUGP * sizeof(float));
    double psurf;
    FILE *meta;
    int status;
    int point;
    /*
     * 检查网格缓冲区分配是否成功。当前代码在后续某些错误返回路径上不会 free(grid)，这里保持原逻辑不变。
     * Check whether the grid buffer allocation succeeded. In the current code, some later error-return paths do not free(grid); the original behavior is preserved here.
     */
    if (grid == NULL)
    {
        return 1;
    }
    /*
     * 初始化球谐变换内核：使用编译期三角截断 NTRU、纬向/经向网格维度以及本地纬度数。最后一个参数 0 保持原程序配置。
     * Initialize the spherical-harmonic transform kernel with compile-time triangular truncation NTRU, latitude/longitude dimensions, and local latitude count. The final argument remains 0 as in the original program.
     */
    status = sht_init(PLASIC_NTRU, PLASIC_NLAT, PLASIC_NLON, PLASIC_NLAT, 0);
    // int sht_init(int ntru, int nlat, int nlon, int nlpp, int latitude_offset)
    if (status != PLASIC_SHT_OK)
    {
        fprintf(stderr, "sht_init failed (%d)\n", status);
        return 1;
    }
    /*
     * 从 restart 中查找名为 "so" 的地形谱记录，并要求其恰好包含 PLASIC_NRSP 个 float。
     * Find the orography spectrum named "so" in the restart file and require exactly PLASIC_NRSP float values.
     */
    if (restart_read_named(restart_path, "so", so, PLASIC_NRSP) != 0)
    {
        return 1;
    }
    /*
     * 把单个标量谱场 so 综合到高斯网格。第一个参数 1 表示这里只处理一个二维标量场，而不是多个垂直层。
     * Synthesize the single scalar spectral field so to the Gaussian grid. The first argument 1 indicates one 2-D scalar field rather than multiple vertical levels.
     */
    if (sht_scalar_to_grid(1, so, grid) != PLASIC_SHT_OK)
    {
        return 1;
    }

    /*
     * 谱空间地形在模式内部按 cv^2 无量纲化；这里逐网格点乘以 cv^2，将其恢复到物理尺度。
     * The model stores orography nondimensionalized by cv^2; multiply every grid point by cv^2 here to restore the physical scale.
     */
    for (point = 0; point < PLASIC_NUGP; ++point)
    {
        grid[point] = (float)(grid[point] * planet_cv() * planet_cv());
    }

    /*
     * 将恢复量纲后的全球地形网格写到 out_grid。
     * Write the dimensionalized global orography grid to out_grid.
     */
    if (write_bin(out_grid, grid, (size_t)PLASIC_NUGP) != 0)
    {
        return 1;
    }

    /*
     * 复现 runtime/plasic_runtime.c 中的参考表面气压修正：psurf = exp(log(101100) - so[0]*cv^2/(sqrt(2)*Rd*Tgr))。so[0] 是地形谱的最低阶系数，cv^2 将模式无量纲地形恢复为物理尺度。
     * Reproduce the reference surface-pressure correction used in runtime/plasic_runtime.c: psurf = exp(log(101100) - so[0]*cv^2/(sqrt(2)*Rd*Tgr)). so[0] is the lowest-order orography spectral coefficient, and cv^2 restores the physical scale of the model orography.
     */
    /*
     * 使用最低阶地形谱系数 so[0] 修正参考表面气压，结果 psurf 的单位为 Pa。
     * Use the lowest-order orography coefficient so[0] to correct the reference surface pressure; psurf is in pascals.
     */
    psurf = exp(log(k_psurf0) -
                (double)so[0] * planet_cv() * planet_cv() /
                    sqrt(2.0) / (k_gascon * k_tgr));
    /*
     * 元数据文本文件记录 psurf、cv 和 ct，便于后续 analyse 与模式运行时使用同一组尺度参数。
     * The metadata text file records psurf, cv, and ct so that later analyse steps and the model runtime can use the same scaling parameters.
     */
    meta = fopen(out_meta, "w");
    if (meta == NULL)
    {
        return 1;
    }
    fprintf(meta, "psurf %.10f\n", psurf);
    fprintf(meta, "cv %.10f\n", planet_cv());
    fprintf(meta, "ct %.10f\n", planet_ct());
    fclose(meta);
    printf("psurf = %.6f Pa, cv = %.6f, ct = %.6f\n", psurf, planet_cv(),
           planet_ct());
    free(grid);

    return 0;
}
/*
 * analyse 子命令：读取物理网格场 T、U、V、Q、PS，构造 PlaSiC 使用的无量纲，再执行球谐分析；风场直接转换为散度谱 sd 与涡度谱 sz。
 * The analyse subcommand reads physical gridded T, U, V, Q, and PS fields, constructs the nondimensional used by PlaSiC, and performs spherical-harmonic analysis; the wind fields are converted directly to divergence spectrum sd and vorticity spectrum sz.
 */
static int command_analyse(int nlev, const char *t_path, const char *u_path,
                           const char *v_path, const char *q_path,
                           const char *ps_path, const char *out_dir)
{
    /*
     * plane 是单层全球网格点数；volume 是 nlev 层三维场元素数；spectral 是编译期 NLEV 个垂直层的谱系数总数。
     * plane is the number of global grid points in one level; volume is the number of elements in an nlev-level 3-D field; spectral is the total number of spectral coefficients for the compile-time NLEV levels.
     */
    const size_t plane = (size_t)PLASIC_NUGP;
    const size_t volume = plane * (size_t)nlev;
    const size_t spectral = (size_t)PLASIC_NRSP * (size_t)PLASIC_NLEV;

    float *t = malloc(volume * sizeof(float));
    float *u = malloc(volume * sizeof(float));
    float *v = malloc(volume * sizeof(float));
    float *q = malloc(volume * sizeof(float));
    float *ps = malloc(plane * sizeof(float));
    /*
     * 中间网格变量：gt=(T-T0)/ct；gq=q*(ps/psurf)；gp=ln(ps/psurf)。这些正是随后送入标量球谐分析的场。
     * Intermediate grid variables: gt=(T-T0)/ct; gq=q*(ps/psurf); gp=ln(ps/psurf). These are the fields passed to scalar spherical-harmonic analysis.
     */
    float *gt = malloc(volume * sizeof(float));
    float *gq = malloc(volume * sizeof(float));
    float *gp = malloc(plane * sizeof(float));

    /*
     * 输出谱变量：st 为温度谱，sd 为散度谱，sz 为涡度谱，sq 为湿度谱，sp 为对数表面气压谱。sp 只有一个二维层，因此只分配 PLASIC_NRSP。
     * Output spectral variables: st is temperature spectrum, sd divergence spectrum, sz vorticity spectrum, sq humidity spectrum, and sp log-surface-pressure spectrum. sp is a single 2-D field, so only PLASIC_NRSP elements are allocated.
     */
    float *st = malloc(spectral * sizeof(float));
    float *sd = malloc(spectral * sizeof(float));
    float *sz = malloc(spectral * sizeof(float));
    float *sq = malloc(spectral * sizeof(float));
    float *sp = malloc((size_t)PLASIC_NRSP * sizeof(float));
    const double ct = planet_ct();
    char path[4096];
    int status;
    size_t index;

    if (t == NULL || u == NULL || v == NULL || q == NULL || ps == NULL ||
        gt == NULL || gq == NULL || gp == NULL || st == NULL || sd == NULL ||
        sz == NULL || sq == NULL || sp == NULL)
    {
        return 1;
    }

    if (nlev != PLASIC_NLEV)
    {
        fprintf(stderr, "grid file levels %d != compiled %d\n", nlev,
                PLASIC_NLEV);
        return 1;
    }

    if (read_bin(t_path, t, volume) || read_bin(u_path, u, volume) ||
        read_bin(v_path, v, volume) || read_bin(q_path, q, volume) ||
        read_bin(ps_path, ps, plane))
    {
        return 1;
    }

    /*
     * 在执行任何球谐分析之前初始化 SHT 内核。
     * Initialize the spherical-harmonic-transform kernel before any spectral analysis.
     */
    status = sht_init(PLASIC_NTRU, PLASIC_NLAT, PLASIC_NLON, PLASIC_NLAT, 0);
    if (status != PLASIC_SHT_OK)
    {
        fprintf(stderr, "sht_init failed (%d)\n", status);
        return 1;
    }

    /*
     * 运行时参考气压通过环境变量 PLASIC_PSURF 传入；通常可使用 orography 子命令输出的 psurf。若该环境变量不存在，则退回未修正的 101100 Pa。
     * The runtime reference pressure is supplied through the PLASIC_PSURF environment variable; the psurf produced by the orography subcommand can normally be used. If the environment variable is absent, the code falls back to the uncorrected 101100 Pa.
     */
    {
        /*
         * 读取可选环境变量 PLASIC_PSURF；atof 将字符串转成 double。若变量不存在则使用 k_psurf0。
         * Read the optional PLASIC_PSURF environment variable; atof converts it to double. If absent, use k_psurf0.
         */
        const char *psurf_env = getenv("PLASIC_PSURF");
        const double psurf =
            psurf_env != NULL ? atof(psurf_env) : k_psurf0;
        /*
         * 三维数组按 level-major 连续存储，每个垂直层包含 plane 个点，因此 index % plane 可得到对应的水平网格点索引，并据此为每一层复用同一个二维 ps[point]。
         * The 3-D arrays are contiguous and level-major, with plane points per level. Therefore index % plane gives the corresponding horizontal point index, allowing the same 2-D ps[point] field to be reused at every vertical level.
         */
        for (index = 0; index < volume; ++index)
        {
            const size_t point = index % plane;
            gt[index] = (float)(((double)t[index] - k_t0) / ct);
            gq[index] = (float)((double)q[index] * (double)ps[point] / psurf);
        }

        for (index = 0; index < plane; ++index)
        {
            gp[index] = (float)log((double)ps[index] / psurf);
        }
    }
    /*
     * 中文：对 gt、gq、gp 做标量球谐分析得到 st、sq、sp；对 u/v 做矢量球谐分析，直接得到散度谱 sd 和相对涡度谱 sz。任一内核失败都统一报告 spherical analysis failed。
     * English: Perform scalar spherical-harmonic analysis on gt, gq, and gp to obtain st, sq, and sp; perform vector analysis on u/v to obtain divergence spectrum sd and relative-vorticity spectrum sz. Any kernel failure is reported as spherical analysis failed.
     */
    if (sht_grid_to_scalar(nlev, gt, st) != PLASIC_SHT_OK ||
        sht_grid_to_scalar(nlev, gq, sq) != PLASIC_SHT_OK ||
        sht_grid_to_scalar(1, gp, sp) != PLASIC_SHT_OK ||
        sht_wind_to_vortdiv(nlev, u, v, sd, sz) != PLASIC_SHT_OK)
    {
        fprintf(stderr, "spherical analysis failed\n");
        return 1;
    }

    add_planetary_vorticity_to_mode_zero_one(sz);

    if (make_out_path(path, sizeof(path), out_dir, "st.bin") ||
        write_bin(path, st, spectral))
    {
        return 1;
    }
    if (make_out_path(path, sizeof(path), out_dir, "sd.bin") ||
        write_bin(path, sd, spectral))
    {
        return 1;
    }
    if (make_out_path(path, sizeof(path), out_dir, "sz.bin") ||
        write_bin(path, sz, spectral))
    {
        return 1;
    }
    if (make_out_path(path, sizeof(path), out_dir, "sq.bin") ||
        write_bin(path, sq, spectral))
    {
        return 1;
    }
    if (make_out_path(path, sizeof(path), out_dir, "sp.bin") ||
        write_bin(path, sp, (size_t)PLASIC_NRSP))
    {
        return 1;
    }
    printf("analysis written to %s\n", out_dir);
    return 0;
}

static int command_synth(const char *dir)
{
    const size_t plane = (size_t)PLASIC_NUGP;
    const size_t volume = plane * (size_t)PLASIC_NLEV;
    const size_t spectral = (size_t)PLASIC_NRSP * (size_t)PLASIC_NLEV;

    float *st = malloc(spectral * sizeof(float));
    float *sd = malloc(spectral * sizeof(float));
    float *sz = malloc(spectral * sizeof(float));
    float *sq = malloc(spectral * sizeof(float));
    float *sp = malloc((size_t)PLASIC_NRSP * sizeof(float));

    float *gt = malloc(volume * sizeof(float));
    float *gq = malloc(volume * sizeof(float));
    float *gp = malloc(plane * sizeof(float));
    float *zero = calloc(volume, sizeof(float));
    float *u = malloc(volume * sizeof(float));
    float *v = malloc(volume * sizeof(float));
    float *t = malloc(volume * sizeof(float));
    float *q = malloc(volume * sizeof(float));
    float *ps = malloc(plane * sizeof(float));
    const double ct = planet_ct();
    char path[4096];
    size_t index;

    if (st == NULL || sd == NULL || sz == NULL || sq == NULL || sp == NULL ||
        gt == NULL || gq == NULL || gp == NULL || zero == NULL || u == NULL ||
        v == NULL || t == NULL || q == NULL || ps == NULL)
    {
        return 1;
    }

    if (make_out_path(path, sizeof(path), dir, "st.bin") ||
        read_bin(path, st, spectral))
    {
        return 1;
    }
    if (make_out_path(path, sizeof(path), dir, "sd.bin") ||
        read_bin(path, sd, spectral))
    {
        return 1;
    }
    if (make_out_path(path, sizeof(path), dir, "sz.bin") ||
        read_bin(path, sz, spectral))
    {
        return 1;
    }
    if (make_out_path(path, sizeof(path), dir, "sq.bin") ||
        read_bin(path, sq, spectral))
    {
        return 1;
    }
    if (make_out_path(path, sizeof(path), dir, "sp.bin") ||
        read_bin(path, sp, (size_t)PLASIC_NRSP))
    {
        return 1;
    }

    if (sht_init(PLASIC_NTRU, PLASIC_NLAT, PLASIC_NLON, PLASIC_NLAT, 0) !=
        PLASIC_SHT_OK)
    {
        return 1;
    }

    if (sht_scalar_to_grid(PLASIC_NLEV, st, gt) != PLASIC_SHT_OK ||
        sht_scalar_to_grid(PLASIC_NLEV, sq, gq) != PLASIC_SHT_OK ||
        sht_scalar_to_grid(1, sp, gp) != PLASIC_SHT_OK ||
        sht_vortdiv_to_wind(PLASIC_NLEV,
                            (float)PLASIC_Y10_PLANETARY_VORTICITY_FACTOR, sd,
                            sz, u, v) != PLASIC_SHT_OK)
    {
        return 1;
    }

    for (index = 0; index < volume; ++index)
    {
        const size_t point = index % plane;
        const double ps_value = exp((double)gp[point]);
        t[index] = (float)((double)gt[index] * ct + k_t0);
        q[index] = (float)((double)gq[index] / ps_value);
    }

    if (make_out_path(path, sizeof(path), dir, "synth_t.bin") ||
        write_bin(path, t, volume))
    {
        return 1;
    }
    if (make_out_path(path, sizeof(path), dir, "synth_u.bin") ||
        write_bin(path, u, volume))
    {
        return 1;
    }
    if (make_out_path(path, sizeof(path), dir, "synth_v.bin") ||
        write_bin(path, v, volume))
    {
        return 1;
    }
    if (make_out_path(path, sizeof(path), dir, "synth_q.bin") ||
        write_bin(path, q, volume))
    {
        return 1;
    }

    for (index = 0; index < plane; ++index)
    {
        ps[index] = (float)exp((double)gp[index]);
    }
    if (make_out_path(path, sizeof(path), dir, "synth_ps.bin") ||
        write_bin(path, ps, plane))
    {
        return 1;
    }

    (void)zero;
    printf("synthesis written to %s\n", dir);
    return 0;
}

int main(int argc, char **argv)
{
    if (argc < 2)
    {
        fprintf(stderr,
                "usage: %s info|orography|analyse|synth ...\n",
                argv[0]);
        return 2;
    }

    /*
     * 中文：info 不需要额外参数。
     * English: info takes no additional arguments.
     */
    if (strcmp(argv[1], "info") == 0)
    {
        return command_info();
    }

    /*
     * orography 需要 3 个参数：restart 输入、地形网格输出、元数据文本输出。
     * orography requires three arguments: restart input, orography-grid output, and metadata-text output.
     */
    if (strcmp(argv[1], "orography") == 0 && argc == 5)
    {
        return command_orography(argv[2], argv[3], argv[4]);
    }

    /*
     * analyse 需要 nlev、T/U/V/Q/PS 五个输入文件和一个输出目录；atoi 将 nlev 从字符串转为 int。
     * analyse requires nlev, five T/U/V/Q/PS input files, and one output directory; atoi converts nlev from text to int.
     */
    if (strcmp(argv[1], "analyse") == 0 && argc == 9)
    {
        return command_analyse(atoi(argv[2]), argv[3], argv[4], argv[5],
                               argv[6], argv[7], argv[8]);
    }

    /*
     * synth 只需要包含谱文件的目录。
     * synth only requires the directory containing the spectral files.
     */
    if (strcmp(argv[1], "synth") == 0 && argc == 3)
    {
        return command_synth(argv[2]);
    }

    /*
     * 如果子命令名称不匹配，或参数个数不满足对应分支，则统一报告 bad arguments。
     * If the subcommand name does not match, or its argument count is incorrect, report bad arguments.
     */
    fprintf(stderr, "bad arguments\n");
    return 2;
}
