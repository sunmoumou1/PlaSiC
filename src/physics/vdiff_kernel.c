#include "vdiff_kernel.h"

#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

/*
 * 本文件实现 PlaSiC 垂直扩散过程的 C 内核。
 *
 * 主要计算内容：
 *   1. 根据垂直风切变、温度层结和水汽修正计算动量/热量扩散系数；
 *   2. 使用半隐式离散和三对角追赶法求解 u、v 的垂直扩散；
 *   3. 在 ndheat > 0 时，将动量耗散造成的机械能损失转化为热量；
 *   4. 使用相同的隐式框架求解比湿 q 和温度 T 的垂直扩散。
 *
 * 数组布局：
 *   所有二维场均按 [level][horizontal] 的逻辑顺序存储，
 *   但在一维内存中 horizontal 是连续维，即：
 *
 *       index = horizontal + level * nhor
 *
 *   level = 0       表示模式顶层；
 *   level = nlev-1  表示最低模式层。
 */

/*
 * 工作区中按完整三维场大小分配的临时数组数量。
 * 每个“场”包含 nhor * nlev 个 float 元素。
 */
enum
{
    VDIFF_FIELD_COUNT = 18
};

/*
 * 将二维逻辑坐标 (horizontal, level) 转换为一维数组下标。
 *
 * 参数：
 *   horizontal：当前水平格点编号，范围 [0, nhor-1]；
 *   level      ：当前垂直层编号，范围 [0, nlev-1]；
 *   nhor       ：水平格点总数。
 *
 * 返回值：
 *   对应的一维数组下标。
 *
 * 显式转换为 size_t，可以避免在数组寻址时使用有符号整数。
 */
static size_t field_index(int horizontal, int level, int nhor)
{
    return (size_t)horizontal + (size_t)level * (size_t)nhor;
}

/*
 * 统一处理工作区内存分配失败。
 * 输出错误信息后直接终止程序，避免后续解引用空指针。
 */
static void allocation_failure(void)
{
    fputs("vdiff: unable to allocate C-kernel workspace\n", stderr);
    abort();
}

/*
 * 计算温度、比湿和水平风的垂直湍流扩散倾向。
 *
 * 参数说明：
 *   vdiff_lamm  ：动量混合长度的渐近上限；
 *   kap        ：κ = R_d / c_pd；
 *   rdbrv       ：R_d / R_v，用于虚温/虚位温的水汽修正；
 *   cpd        ：干空气定压比热 c_pd；
 *   cpv_cpd_minus1         ：水汽对空气定压比热的修正常数；
 *   sigma       ：模式全层的 σ 坐标；
 *   sigmah      ：模式半层（层界面）的 σ 坐标；
 *   dsigma      ：各模式层的 σ 厚度；
 *   dt, dq      ：输入温度和比湿；
 *   du, dv      ：输入纬向风和经向风；
 *   dtdt, dqdt  ：输入/输出温度和比湿总倾向；
 *   dudt, dvdt  ：输入/输出纬向风和经向风总倾向。
 */
void vdiff(
    int nhor,
    int nlev,
    int ndheat,
    float zumin,
    float vdiff_lamm,
    float vdiff_b,
    float vdiff_c,
    float vdiff_d,
    float deltsec2,
    float ga,
    float gascon,
    float kap,
    float rdbrv,
    float cpd,
    float cpv_cpd_minus1,
    const float *sigma,
    const float *sigmah,
    const float *dsigma,
    const float *dt,
    const float *dq,
    const float *du,
    const float *dv,
    float *dtdt,
    float *dqdt,
    float *dudt,
    float *dvdt)
{
    /*
     * ztscal：用于由 σ 坐标估算几何高度和层厚的参考温度，单位 K。
     * 这里采用常数 250 K，属于扩散参数化中的尺度温度。
     */
    const float ztscal = 250.0f;

    /* von Kármán 常数，用于构造近地层混合长度。 */
    const float vonkarman = 0.4f;

    /*
     * 相邻全层之间共有 nlev-1 个内部界面，因此 nlem = nlev-1。
     * 后续扩散系数只在这些相邻层界面上计算。
     */
    const int nlem = nlev - 1;

    /* 一个完整二维场包含的 float 元素数量。 */
    const size_t field_size = (size_t)nhor * (size_t)nlev;

    /*
     * 一次性工作区总元素数：
     *   - 18 个大小为 nhor*nlev 的二维临时场；
     *   - zskap 和 zskaph 两个长度为 nlev 的垂直数组。
     *
     * 使用单次 malloc，可减少多次分配开销，并使临时数据连续存储。
     */
    const size_t work_count = (size_t)VDIFF_FIELD_COUNT * field_size + 2U * (size_t)nlev;
    // u 表示 unsigned int

    /* 整块工作区的起始地址。 */
    float *work;

    /* 在整块工作区内依次切分各临时数组时使用的游标。 */
    float *cursor;

    /* 以下四个数组保存本过程单独产生的 T、u、v、q 倾向。 */
    float *zdtdt;
    float *zdudt;
    float *zdvdt;
    float *zdqdt;

    /*
     * zkdiffm：动量垂直扩散系数，供 u、v 和动能扩散使用；
     * zkdiffh：标量垂直扩散系数，供热量和水汽扩散使用。
     *
     * 二者在 level 位置表示 level 与 level+1 之间界面的耦合强度。
     */
    float *zkdiffm;
    float *zkdiffh;

    /* 水平风速模长 sqrt(u^2+v^2)。 */
    float *zabsu;

    /* 隐式扩散求解后的新温度、新风场和新比湿。 */
    float *ztn;
    float *zun;
    float *zvn;
    float *zqn;

    /*
     * 在进入垂直扩散前，将已有倾向推进 deltsec2 后得到的工作场。
     * 它们是本函数进行半隐式扩散求解时所使用的初始状态。
     */
    float *zt;
    float *zu;
    float *zv;
    float *zq;

    /*
     * 三对角追赶法中的消元系数。
     * 同一块数组会在风、动能、水汽和温度求解阶段重复使用。
     */
    float *zebs;

    /* 扩散前和扩散后的单位质量水平风动能。 */
    float *zke;
    float *zken;

    /*
     * σ^κ：
     *   zskap[level]  = sigma[level]^kap；
     *   zskaph[level] = sigmah[level]^kap。
     */
    float *zskap;
    float *zskaph;

    /* 热量/水汽所用混合长度的渐近上限。 */
    float zlamh;

    /* 将时间步长、重力和气体常数组合起来的预计算常数。 */
    float zkonst1;
    float zkonst2;
    float zkonst3;

    /* 垂直层循环和水平格点循环下标。 */
    int level;
    int horizontal;

    /* 为全部临时变量一次性分配连续工作区。 */
    work = malloc(work_count * sizeof(*work));
    if (work == NULL)
    {
        allocation_failure();
    }

    /*
     * 从 work 起始位置开始，依次为每个二维临时场分配一段地址。
     * TAKE_FIELD 只移动指针，不执行额外的动态内存分配。
     */
    cursor = work;

#define TAKE_FIELD(name)      \
    do                        \
    {                         \
        name = cursor;        \
        cursor += field_size; \
    } while (0)
    TAKE_FIELD(zdtdt);
    TAKE_FIELD(zdudt);
    TAKE_FIELD(zdvdt);
    TAKE_FIELD(zdqdt);
    TAKE_FIELD(zkdiffm);
    TAKE_FIELD(zkdiffh);
    TAKE_FIELD(zabsu);
    TAKE_FIELD(ztn);
    TAKE_FIELD(zun);
    TAKE_FIELD(zvn);
    TAKE_FIELD(zqn);
    TAKE_FIELD(zt);
    TAKE_FIELD(zu);
    TAKE_FIELD(zv);
    TAKE_FIELD(zq);
    TAKE_FIELD(zebs);
    TAKE_FIELD(zke);
    TAKE_FIELD(zken);
#undef TAKE_FIELD

    /*
     * 18 个完整二维场之后，剩余空间依次分给两个长度为 nlev 的数组。
     */
    zskap = cursor;
    cursor += nlev;
    zskaph = cursor;

    /*
     * 热量/水汽混合长度上限。
     */
    zlamh = vdiff_lamm * sqrtf(3.0f * vdiff_d * 0.5f);

    /*
     * 预先合并重复出现的常数因子：
     *
     *   zkonst1 = g * (2Δt) / R_d
     *   zkonst2 = g² * (2Δt) / R_d²
     *   zkonst3 = g³ * (2Δt) / R_d³
     *
     * 其中 deltsec2 = 2Δt。
     */
    zkonst1 = ga * deltsec2 / gascon;
    zkonst2 = zkonst1 * ga / gascon;
    zkonst3 = zkonst2 * ga / gascon;

    /*
     * 预计算全层和半层上的 σ^κ，避免在水平格点循环中重复调用 powf。
     */
    for (level = 0; level < nlev; ++level)
    {
        zskap[level] = powf(sigma[level], kap);
        zskaph[level] = powf(sigmah[level], kap);
    }

    /*
     * ------------------------------------------------------------------
     * 第 1 步：构造进入垂直扩散求解的工作状态
     * ------------------------------------------------------------------
     *
     * 将调用本函数之前已经累积的物理倾向推进 deltsec2：
     *
     *   X_work = X_input + (dX/dt)_existing * deltsec2
     *
     * 随后使用更新后的 u、v 计算水平风速模长。
     */
    for (level = 0; level < nlev; ++level)
    {
        for (horizontal = 0; horizontal < nhor; ++horizontal)
        {
            const size_t index = field_index(horizontal, level, nhor);
            zt[index] = dt[index] + dtdt[index] * deltsec2;
            zq[index] = dq[index] + dqdt[index] * deltsec2;
            zu[index] = du[index] + dudt[index] * deltsec2;
            zv[index] = dv[index] + dvdt[index] * deltsec2;
            zabsu[index] = sqrtf(zu[index] * zu[index] + zv[index] * zv[index]);
        }
    }

    /*
     * ------------------------------------------------------------------
     * 第 2 步：计算各相邻模式层之间的动量和标量扩散系数
     * ------------------------------------------------------------------
     *
     * level 对应 level 与 level+1 之间的内部界面。
     */
    for (level = 0; level < nlem; ++level)
    {
        const int next_level = level + 1;

        /*
         * 使用等温静力近似，由半层 σ 值估算该界面的几何高度：
         *
         *   z = -(R_d T_scale / g) ln(σ_h)
         *
         * σ_h 越小，表示位置越高，zzlev 越大。
         */
        const float zzlev = -gascon * ztscal * logf(sigmah[level]) / ga;

        /*
         * 使用相邻全层 σ 比值估算两层之间的几何厚度：
         *
         *   Δz = (R_d T_scale / g) ln(σ_{k+1}/σ_k)
         */
        const float zdz = gascon * ztscal * logf(sigma[next_level] / sigma[level]) / ga;

        /*
         * Blackadar 型混合长度：
         *
         *   l = l_inf κz / (l_inf + κz)
         *
         * 近地面时 l≈κz；高度较大时逐渐趋近渐近上限 l_inf。
         * zmixh 用于热量/水汽，zmixm 用于动量。
         */
        const float zmixh = zlamh * vonkarman * zzlev / (zlamh + vonkarman * zzlev);
        const float zmixm = vdiff_lamm * vonkarman * zzlev / (vdiff_lamm + vonkarman * zzlev);

        /* 扩散系数中使用混合长度平方，提前计算以避免重复乘法。 */
        const float zmixh2 = zmixh * zmixh;
        const float zmixm2 = zmixm * zmixm;

        /* 相邻全层 σ 间距的倒数。 */
        const float zrdsig = 1.0f / (sigma[next_level] - sigma[level]);

        /*
         * 扩散方程离散后出现的公共系数，包含：
         *   时间步长、重力/气体常数比例、半层 σ³ 和 1/Δσ。
         */
        const float zfac = zkonst3 * sigmah[level] * sigmah[level] * sigmah[level] * zrdsig;

        /*
         * 以下量用于不稳定层结的稳定度函数。
         * 它们只依赖当前垂直界面的几何尺度，因此放在水平循环外计算。
         */
        const float zratio = (zzlev + zdz) / zzlev;
        const float zroot = powf(zratio, 1.0f / 3.0f);
        const float zcube_base = (zroot - 1.0f) / zdz;
        const float zcube = zcube_base * zcube_base * zcube_base;
        const float zzkfac = sqrtf(zcube / zzlev);

        for (horizontal = 0; horizontal < nhor; ++horizontal)
        {
            const size_t index = field_index(horizontal, level, nhor);
            const size_t next_index = field_index(horizontal, next_level,
                                                  nhor);

            /*
             * 两层温度的 dsigma 加权倒数。
             * 可理解为界面附近 1/T 的离散近似，后续以三次方进入扩散系数。
             */
            const float zrthl = (dsigma[next_level] + dsigma[level]) / (zt[index] * dsigma[next_level] + zt[next_index] * dsigma[level]);

            /*
             * 相邻两层“风速模长”的差，而不是 u、v 矢量差。
             * fabsf 保证差值非负。
             */
            const float wind_difference = fabsf(zabsu[index] - zabsu[next_index]);

            /*
             * σ 坐标中的垂直风切变尺度：
             *
             *   zdvds = max(zumin, | |V_k|-|V_{k+1}| |) / Δσ
             *
             * zumin 防止风切变严格为零，避免 Richardson 数分母为零。
             */
            const float zdvds = fmaxf(zumin, wind_difference) * zrdsig;

            /*
             * 虚温修正因子：
             *
             *   1 + (R_v/R_d - 1) q
             *     = 1 + (1/rdbrv - 1) q
             *
             * 用来考虑水汽对空气密度和静力稳定度的影响。
             */
            const float zqf1 = 1.0f + (1.0f / rdbrv - 1.0f) * zq[index];
            const float zqf2 = 1.0f + (1.0f / rdbrv - 1.0f) * zq[next_index];

            /*
             * 虚位温在 σ 坐标中的层间梯度近似：
             *
             *   θ_v ≈ T [1 + (R_v/R_d-1)q] / σ^κ
             *
             *   zdthds ≈ (θ_{v,k} - θ_{v,k+1}) / Δσ
             */
            const float zdthds = (zt[index] * zqf1 / zskap[level] - zt[next_index] * zqf2 / zskap[next_level]) * zrdsig;

            /*
             * 梯度 Richardson 数的 σ 坐标形式。
             *
             *   Ri ∝ 静力稳定度 / 风切变²
             *
             * Ri <= 0：不稳定或中性层结，湍流混合被增强；
             * Ri >  0：稳定层结，湍流混合被抑制。
             */
            const float zri = zskaph[level] * gascon * zdthds / (sigmah[level] * zdvds * zdvds);

            /* 热量/水汽和动量对应的稳定度修正函数。 */
            float zrifh;
            float zrifm;

            if (zri <= 0.0f)
            {
                const float unstable_zri = fminf(zri, 0.0f);

                /* sqrt(-Ri)，仅用于不稳定层结修正。 */
                const float zsqrt = sqrtf(-1.0f * unstable_zri);

                /*
                 * 不稳定条件下，热量/水汽和动量使用不同的经验分母。
                 * 混合长度平方越大、层结越不稳定，修正幅度通常越强。
                 */
                const float zdenomh = 1.0f + 3.0f * vdiff_c * vdiff_b * zsqrt * zzkfac * zmixh2;
                const float zdenomm = 1.0f + 3.0f * vdiff_c * vdiff_b * zsqrt * zzkfac * zmixm2;

                /* Ri 为负，因此下面两式通常使稳定度函数大于 1。 */
                zrifh = 1.0f - 3.0f * vdiff_b * unstable_zri / zdenomh;
                zrifm = 1.0f - 2.0f * vdiff_b * unstable_zri / zdenomm;
            }
            else
            {
                const float stable_zri = fmaxf(zri, 0.0f);

                /* 稳定层结经验函数中的 sqrt(1 + d Ri)。 */
                const float zdenom = sqrtf(1.0f + vdiff_d * stable_zri);

                /*
                 * Ri 增大时，zrifh 和 zrifm 均减小，从而抑制垂直湍流交换。
                 * 热量/水汽和动量使用不同的系数形式。
                 */
                zrifh = 1.0f / (1.0f + 3.0f * vdiff_b * stable_zri * zdenom);
                zrifm = 1.0f / (1.0f + 2.0f * vdiff_b * stable_zri / zdenom);
            }

            /*
             * 最终标量和动量扩散耦合系数：
             *
             *   K_discrete ∝ zfac × l² × f(Ri) × |垂直风切变| × (1/T)³
             *
             * 这些系数已经包含时间离散和 σ 坐标变换所需的因子，
             * 后面直接作为相邻层之间的三对角耦合系数使用。
             */
            zkdiffh[index] = zfac * zmixh2 * zrifh * zdvds * zrthl * zrthl * zrthl;
            zkdiffm[index] = zfac * zmixm2 * zrifm * zdvds * zrthl * zrthl * zrthl;
        }
    }

    /*
     * ------------------------------------------------------------------
     * 第 3 步：隐式求解 u、v 的垂直扩散
     * ------------------------------------------------------------------
     *
     * 对每个水平格点，垂直方向构成一个三对角线性系统。
     * 以下代码采用 Thomas 追赶法解三对角矩阵：
     *   1. 顶层初始化；
     *   2. 从上到下前向消元；
     *   3. 最低层闭合；
     *   4. 从下到上回代。
     *
     * 由于 u 和 v 的扩散系数完全相同，因此我们在同一轮消元中同时求解。
     */

    /* 
     * 【1. 顶层方程初始化 (level = 0)】
     * 因为顶层之上没有大气层，所以来自上方的大气动量通量项 $K_{-1} = 0$
     * 我们需要将其化为标准递推形式：$u_0^{new} = E_0 u_1^{new} + F_0$
     */
    for (horizontal = 0; horizontal < nhor; ++horizontal)
    {
        const size_t top = field_index(horizontal, 0, nhor);
        
        // 提取方程左侧合并后的公分母：Denominator = $\Delta\sigma_0 + K_0$
        const float denominator = dsigma[0] + zkdiffm[top];
        
        // 计算消元系数 E_0 = $K_0 / Denominator$
        zebs[top] = zkdiffm[top] / denominator;
        
        // 计算临时常数项 F_0 = $(\Delta\sigma_0 u_0^{old}) / Denominator$
        // 注意：此处 zun 和 zvn 数组被临时借用来存储 F 变量，它们此时还不是最终的风速。
        zun[top] = dsigma[0] * zu[top] / denominator;
        zvn[top] = dsigma[0] * zv[top] / denominator;
    }

    /* 
     * 【2. 从第 2 层向下进行前向消元 (level = 1 到 nlem - 1)】
     */
    for (level = 1; level < nlem; ++level)
    {
        const int previous_level = level - 1;
        for (horizontal = 0; horizontal < nhor; ++horizontal)
        {
            const size_t index = field_index(horizontal, level, nhor);
            const size_t previous = field_index(horizontal, previous_level, nhor);

            // 分母 Denominator = $\Delta\sigma_k + K_k + K_{k-1} * (1 - E_{k-1})$
            // 这里的 zebs[previous] 就是上一层算好的 $E_{k-1}$
            const float denominator = dsigma[level] + zkdiffm[index] + zkdiffm[previous] * (1.0f - zebs[previous]);

            // 计算当前层的消元系数 E_k = $K_k / Denominator$
            zebs[index] = zkdiffm[index] / denominator;

            // 计算当前层的临时常数项 F_k = $(\Delta\sigma_k u_k^{old} + K_{k-1} F_{k-1}) / Denominator$
            // 此时上一层的 zun[previous] 里恰好存储着 $F_{k-1}$，直接拿来用。
            // 算出的新 F_k 继续存放在 zun[index] 中供下一层消元使用。
            zun[index] = (zu[index] * dsigma[level] + zkdiffm[previous] * zun[previous]) / denominator;
            zvn[index] = (zv[index] * dsigma[level] + zkdiffm[previous] * zvn[previous]) / denominator;
        }
    }

    /*
     * 【3. 最低层方程闭合 (level = nlev - 1)】
     */
    for (horizontal = 0; horizontal < nhor; ++horizontal)
    {
        const size_t bottom = field_index(horizontal, nlev - 1, nhor);
        const size_t previous = field_index(horizontal, nlem - 1, nhor);
        
        // 分母 Denominator = $\Delta\sigma_{bottom} + K_{k-1} * (1 - E_{k-1})$ (因为 $K_k=0$，被去掉了)
        const float denominator = dsigma[nlev - 1] + zkdiffm[previous] * (1.0f - zebs[previous]);
        
        // $u_{bottom}^{new} = (\Delta\sigma_{bottom} u_{bottom}^{old} + K_{k-1} F_{k-1}) / Denominator$
        // 此时存入 zun[bottom] 的不再是临时变量 F，而是真真实实的底层隐式新风速了。
        zun[bottom] = (zu[bottom] * dsigma[nlev - 1] + zkdiffm[previous] * zun[previous]) / denominator;
        zvn[bottom] = (zv[bottom] * dsigma[nlev - 1] + zkdiffm[previous] * zvn[previous]) / denominator;
    }

    /*
     * 【4. 自底向上回代 (Backward substitution)】
     * 我们已经求出了最底层确定的新风速。
     * 根据我们之前建立并保存的关系：$u_k^{new} = E_k u_{k+1}^{new} + F_k$
     * 我们从倒数第二层向顶层依次倒推，利用已知的下一层风速，求解当前层的风速。
     */
    for (level = nlem - 1; level >= 0; --level)
    {
        const int next_level = level + 1;
        for (horizontal = 0; horizontal < nhor; ++horizontal)
        {
            const size_t index = field_index(horizontal, level, nhor);
            const size_t next = field_index(horizontal, next_level, nhor);
            
            // 公式：$u_k^{new} = F_k + E_k * u_{k+1}^{new}$
            // 此时，zun[index] 里存的正是之前计算好的常数 $F_k$；
            // zebs[index] 是 $E_k$；
            // zun[next] 则是刚刚在上一轮循环(或底层闭合)中算好的确定风速 $u_{k+1}^{new}$。
            // 赋值完成后，zun[index] 正式从“临时常数 F”被刷新为了“物理上的隐式新风速”。
            zun[index] = zun[index] + zebs[index] * zun[next];
            zvn[index] = zvn[index] + zebs[index] * zvn[next];
        }
    }

    /*
     * 【5. 将隐式求得的新风场转换为垂直扩散倾向】
     */
    for (level = 0; level < nlev; ++level)
    {
        for (horizontal = 0; horizontal < nhor; ++horizontal)
        {
            const size_t index = field_index(horizontal, level, nhor);
            
            // 计算由于垂直扩散引起的独立风速变化率
            zdudt[index] = (zun[index] - zu[index]) / deltsec2;
            zdvdt[index] = (zvn[index] - zv[index]) / deltsec2;
            
            // 将该变化率累加到输入输出参数(总倾向)中
            dudt[index] = dudt[index] + zdudt[index];
            dvdt[index] = dvdt[index] + zdvdt[index];
        }
    }

    /*
     * ------------------------------------------------------------------
     * 第 4 步：可选的动量耗散加热
     * ------------------------------------------------------------------
     *
     * 垂直扩散会削弱风切变并耗散机械能。
     * 当 ndheat > 0 时，将不能由动能垂直输送解释的风动能损失
     * 转换成温度倾向，以改善能量守恒。
     */
    if (ndheat > 0)
    {
        /* 扩散前单位质量水平风动能：KE = 0.5*(u^2+v^2)。 */
        for (level = 0; level < nlev; ++level)
        {
            for (horizontal = 0; horizontal < nhor; ++horizontal)
            {
                const size_t index = field_index(horizontal, level, nhor);
                zke[index] = 0.5f * zabsu[index] * zabsu[index];
            }
        }

        /*
         * 使用与动量相同的 zkdiffm，对动能标量执行一次隐式垂直扩散。
         * 这样可以区分“动能在层间的重新分配”和“真正的机械能耗散”。
         */

        /* 顶层初始化。 */
        for (horizontal = 0; horizontal < nhor; ++horizontal)
        {
            const size_t top = field_index(horizontal, 0, nhor);
            const float denominator = dsigma[0] + zkdiffm[top];
            zebs[top] = zkdiffm[top] / denominator;
            zken[top] = dsigma[0] * zke[top] / denominator;
        }

        /* 内部层前向消元。 */
        for (level = 1; level < nlem; ++level)
        {
            const int previous_level = level - 1;
            for (horizontal = 0; horizontal < nhor; ++horizontal)
            {
                const size_t index = field_index(horizontal, level, nhor);
                const size_t previous = field_index(horizontal,
                                                    previous_level, nhor);
                const float denominator = dsigma[level] + zkdiffm[index] + zkdiffm[previous] * (1.0f - zebs[previous]);
                zebs[index] = zkdiffm[index] / denominator;
                zken[index] = (zke[index] * dsigma[level] + zkdiffm[previous] * zken[previous]) / denominator;
            }
        }

        /* 最低层闭合。 */
        for (horizontal = 0; horizontal < nhor; ++horizontal)
        {
            const size_t bottom = field_index(horizontal, nlev - 1, nhor);
            const size_t previous = field_index(horizontal, nlem - 1, nhor);
            const float denominator = dsigma[nlev - 1] + zkdiffm[previous] * (1.0f - zebs[previous]);
            zken[bottom] = (zke[bottom] * dsigma[nlev - 1] + zkdiffm[previous] * zken[previous]) / denominator;
        }

        /* 自底向上回代，得到纯“动能扩散”后的 zken。 */
        for (level = nlem - 1; level >= 0; --level)
        {
            const int next_level = level + 1;
            for (horizontal = 0; horizontal < nhor; ++horizontal)
            {
                const size_t index = field_index(horizontal, level, nhor);
                const size_t next = field_index(horizontal, next_level,
                                                nhor);
                zken[index] = zken[index] + zebs[index] * zken[next];
            }
        }

        /* 耗散加热来自动能输送倾向与风场动能倾向之差。 */
        for (level = 0; level < nlev; ++level)
        {
            for (horizontal = 0; horizontal < nhor; ++horizontal)
            {
                const size_t index = field_index(horizontal, level, nhor);
                const float wind_energy_tendency =
                    0.5f * (zun[index] * zun[index] - zu[index] * zu[index] +
                            zvn[index] * zvn[index] - zv[index] * zv[index]) /
                    deltsec2;
                const float transported_energy_tendency =
                    (zken[index] - zke[index]) / deltsec2;

                zdtdt[index] =
                    (transported_energy_tendency - wind_energy_tendency) /
                    cpd / (1.0f + cpv_cpd_minus1 * dq[index]);

                /* 将耗散加热累加到总温度倾向。 */
                dtdt[index] = dtdt[index] + zdtdt[index];

                /*
                 * 清零临时数组当前元素。
                 * 该值稍后会被温度扩散倾向重新覆盖。
                 */
                zdtdt[index] = 0.0f;
            }
        }
    }

    /*
     * ------------------------------------------------------------------
     * 第 5 步：隐式求解比湿 q 的垂直扩散
     * ------------------------------------------------------------------
     *
     * q 是普通标量，直接使用 zkdiffh 构造与动量类似的三对角系统。
     */

    /* 顶层初始化。 */
    for (horizontal = 0; horizontal < nhor; ++horizontal)
    {
        const size_t top = field_index(horizontal, 0, nhor);
        const float denominator = dsigma[0] + zkdiffh[top];
        zebs[top] = zkdiffh[top] / denominator;
        zqn[top] = dsigma[0] * zq[top] / denominator;
    }

    /* 内部层前向消元。 */
    for (level = 1; level < nlem; ++level)
    {
        const int previous_level = level - 1;
        for (horizontal = 0; horizontal < nhor; ++horizontal)
        {
            const size_t index = field_index(horizontal, level, nhor);
            const size_t previous = field_index(horizontal, previous_level,
                                                nhor);
            const float denominator = dsigma[level] + zkdiffh[index] + zkdiffh[previous] * (1.0f - zebs[previous]);
            zebs[index] = zkdiffh[index] / denominator;
            zqn[index] = (zq[index] * dsigma[level] + zkdiffh[previous] * zqn[previous]) / denominator;
        }
    }

    /* 最低层闭合。 */
    for (horizontal = 0; horizontal < nhor; ++horizontal)
    {
        const size_t bottom = field_index(horizontal, nlev - 1, nhor);
        const size_t previous = field_index(horizontal, nlem - 1, nhor);
        const float denominator = dsigma[nlev - 1] + zkdiffh[previous] * (1.0f - zebs[previous]);
        zqn[bottom] = (zq[bottom] * dsigma[nlev - 1] + zkdiffh[previous] * zqn[previous]) / denominator;
    }

    /* 自底向上回代，得到各层扩散后的新比湿。 */
    for (level = nlem - 1; level >= 0; --level)
    {
        const int next_level = level + 1;
        for (horizontal = 0; horizontal < nhor; ++horizontal)
        {
            const size_t index = field_index(horizontal, level, nhor);
            const size_t next = field_index(horizontal, next_level, nhor);
            zqn[index] = zqn[index] + zebs[index] * zqn[next];
        }
    }

    /* 将新旧比湿之差转换成垂直扩散倾向并累加。 */
    for (level = 0; level < nlev; ++level)
    {
        for (horizontal = 0; horizontal < nhor; ++horizontal)
        {
            const size_t index = field_index(horizontal, level, nhor);
            zdqdt[index] = (zqn[index] - zq[index]) / deltsec2;
            dqdt[index] = dqdt[index] + zdqdt[index];
        }
    }

    /*
     * ------------------------------------------------------------------
     * 第 6 步：隐式求解温度 T 的垂直扩散
     * ------------------------------------------------------------------
     *
     * 模式保存和预报的是完整层温度 T_k，但湍流热通量按照位温
     *
     *     theta_k = T_k / Pi_k,       Pi_k = sigma_k^kappa
     *
     * 的垂直差计算。这里有两套位置不同、不能相互约掉的转换因子：
     *
     *   zskap[k]  = Pi_k   = sigma_k^kappa
     *       定义在完整模式层 k，用于把该层温度 T_k 转成位温 theta_k；
     *
     *   zskaph[j] = Pi_h,j = sigma_h,j^kappa
     *       定义在界面 j（完整层 j 和 j+1 之间），用于构造界面热通量系数。
     *
     * 本步骤开始前，zkdiffh[j] 是已经用于比湿扩散的原始界面扩散系数 K_j。
     * 比湿求解已经结束，因此可以原地把它改成温度方程使用的界面耦合系数：
     *
     *     H_j = K_j * Pi_h,j
     *
     * 后续看到的 zkdiffh[j] 都表示 H_j，而不再表示原始 K_j。
     * 相邻两层之间的热扩散项于是具有下面的结构：
     *
     *     H_j * (T_{j+1}/Pi_{j+1} - T_j/Pi_j)
     *
     * 因此代码中既会出现“乘 zskaph[j]”构造 H_j，也会出现
     * “除 zskap[k]”把完整层温度转换为位温。这两个操作分别发生在界面和
     * 完整层，物理位置与用途都不同。
     *
     * 对内部完整层 k，尚未进行前向消元的三对角方程可以写成：
     *
     *   -H_{k-1}/Pi_{k-1} * T_{k-1}^{new}
     *   +[dsigma_k + (H_{k-1}+H_k)/Pi_k] * T_k^{new}
     *   -H_k/Pi_{k+1} * T_{k+1}^{new}
     *   = dsigma_k * T_k^{old}.
     *
     * 前向消元把每一层保存成下面的递推形式：
     *
     *     T_k^{new} = ztn_k + zebs_k * T_{k+1}^{new}/Pi_{k+1}.
     *
     * 此时 ztn_k 只是消元后的常数项，并非最终温度；最终温度要等最低层
     * 求出之后，再通过自底向上的回代得到。
     */

    /*
     * 6.1 构造温度方程的界面耦合系数 H_j = K_j * Pi_h,j。
     *
     * level 在这一循环中表示“界面编号”。界面 level 位于完整层 level 和
     * level+1 之间。只存在 nlem=nlev-1 个内部界面，所以循环不处理 nlev-1。
     */
    for (level = 0; level < nlem; ++level)
    {
        for (horizontal = 0; horizontal < nhor; ++horizontal)
        {
            const size_t interface_index =
                field_index(horizontal, level, nhor);

            zkdiffh[interface_index] =
                zkdiffh[interface_index] * zskaph[level];
        }
    }

    /*
     * 6.2 顶层方程初始化。
     *
     * 顶层上方没有界面，只与其下方的界面 0 相连。顶层方程为：
     *
     *   [dsigma_0 + H_0/Pi_0] * T_0^{new}
     *     = dsigma_0*T_0^{old} + H_0*T_1^{new}/Pi_1.
     *
     * 除以主对角系数 D_0 后：
     *
     *   T_0^{new}
     *     = dsigma_0*T_0^{old}/D_0
     *       + (H_0/D_0)*T_1^{new}/Pi_1.
     *
     * 因此 ztn[0] 保存第一项，zebs[0] 保存 H_0/D_0；Pi_1 留到回代时
     * 与下一层温度一起处理。
     */
    for (horizontal = 0; horizontal < nhor; ++horizontal)
    {
        const size_t top_level_index = field_index(horizontal, 0, nhor);

        /* 界面 0 的存储偏移与完整层 0 相同，但二者代表不同物理位置。 */
        const size_t interface_below_top_index = top_level_index;
        const float top_sigma_kappa = zskap[0];
        const float interface_heat_coupling =
            zkdiffh[interface_below_top_index];
        const float denominator =
            dsigma[0] + interface_heat_coupling / top_sigma_kappa;

        zebs[top_level_index] = interface_heat_coupling / denominator;
        ztn[top_level_index] =
            dsigma[0] * zt[top_level_index] / denominator;
    }

    /*
     * 6.3 内部层前向消元。
     *
     * 对当前完整层 k：
     *
     *   H_{k-1} 是上方界面耦合系数；
     *   H_k     是下方界面耦合系数；
     *   Pi_{k-1} 和 Pi_k 分别属于上方完整层和当前完整层。
     *
     * 上一层已经写成：
     *
     *   T_{k-1}^{new}
     *     = ztn_{k-1} + zebs_{k-1}*T_k^{new}/Pi_k.
     *
     * 将它代入当前层三对角方程后，当前层主对角系数变为：
     *
     *   D_k = dsigma_k
     *         + [H_k + H_{k-1}(1-zebs_{k-1}/Pi_{k-1})]/Pi_k.
     *
     * 消元后的右端常数项为：
     *
     *   R_k = dsigma_k*T_k^{old}
     *         + H_{k-1}*ztn_{k-1}/Pi_{k-1}.
     *
     * 最后保存 ztn_k=R_k/D_k、zebs_k=H_k/D_k，供下一层继续消元。
     */
    for (level = 1; level < nlem; ++level)
    {
        const int level_above = level - 1;

        for (horizontal = 0; horizontal < nhor; ++horizontal)
        {
            const size_t current_level_index =
                field_index(horizontal, level, nhor);
            const size_t level_above_index =
                field_index(horizontal, level_above, nhor);

            /*
             * 界面数组按其上方完整层编号存储，因此：
             *   下方界面 H_k     的偏移等于 current_level_index；
             *   上方界面 H_{k-1} 的偏移等于 level_above_index。
             */
            const size_t interface_below_index = current_level_index;
            const size_t interface_above_index = level_above_index;
            const float current_sigma_kappa = zskap[level];
            const float above_sigma_kappa = zskap[level_above];
            const float coupling_below = zkdiffh[interface_below_index];
            const float coupling_above = zkdiffh[interface_above_index];
            const float above_elimination_correction =
                1.0f - zebs[level_above_index] / above_sigma_kappa;
            const float denominator =
                dsigma[level]
                + (coupling_below
                   + coupling_above * above_elimination_correction)
                  / current_sigma_kappa;
            const float eliminated_right_hand_side =
                zt[current_level_index] * dsigma[level]
                + coupling_above / above_sigma_kappa
                  * ztn[level_above_index];

            zebs[current_level_index] = coupling_below / denominator;
            ztn[current_level_index] =
                eliminated_right_hand_side / denominator;
        }
    }

    {
        /*
         * 6.4 最低层方程闭合。
         *
         * 最低层 b=nlev-1 的下方没有界面，只保留上方界面 H_{b-1}。
         * 上一层消元结果代入后：
         *
         *   D_b = dsigma_b
         *         + H_{b-1}/Pi_b*(1-zebs_{b-1}/Pi_{b-1}),
         *
         *   R_b = dsigma_b*T_b^{old}
         *         + H_{b-1}*ztn_{b-1}/Pi_{b-1},
         *
         *   T_b^{new} = R_b/D_b.
         *
         * 注意 Pi_b 和 Pi_{b-1} 属于两个不同完整层，不能使用同一个下标。
         */
        const int bottom_level = nlev - 1;
        const int level_above_bottom = bottom_level - 1;

        /* 最低层上方的界面编号等于其上方完整层编号。 */
        const int interface_above_bottom = nlem - 1;

        for (horizontal = 0; horizontal < nhor; ++horizontal)
        {
            const size_t bottom_level_index =
                field_index(horizontal, bottom_level, nhor);
            const size_t level_above_index =
                field_index(horizontal, level_above_bottom, nhor);
            const size_t interface_above_index =
                field_index(horizontal, interface_above_bottom, nhor);
            const float bottom_sigma_kappa = zskap[bottom_level];
            const float above_sigma_kappa = zskap[level_above_bottom];
            const float coupling_above = zkdiffh[interface_above_index];
            const float above_elimination_correction =
                1.0f - zebs[level_above_index] / above_sigma_kappa;
            const float denominator =
                dsigma[bottom_level]
                + coupling_above / bottom_sigma_kappa
                  * above_elimination_correction;
            const float eliminated_right_hand_side =
                zt[bottom_level_index] * dsigma[bottom_level]
                + coupling_above * ztn[level_above_index]
                  / above_sigma_kappa;

            ztn[bottom_level_index] =
                eliminated_right_hand_side / denominator;
        }
    }

    /*
     * 6.5 自底向上回代。
     *
     * 前向消元保存的是：
     *
     *   T_k^{new} = ztn_k + zebs_k*T_{k+1}^{new}/Pi_{k+1}.
     *
     * 所以必须从已经求出的最低层开始向上恢复各层温度。这里除以
     * lower_sigma_kappa，是把下一完整层的温度 T_{k+1} 转换成该层位温；
     * 它不是对界面耦合系数的再次修正。
     */
    for (level = nlem - 1; level >= 0; --level)
    {
        const int level_below = level + 1;

        for (horizontal = 0; horizontal < nhor; ++horizontal)
        {
            const size_t current_level_index =
                field_index(horizontal, level, nhor);
            const size_t level_below_index =
                field_index(horizontal, level_below, nhor);
            const float below_sigma_kappa = zskap[level_below];

            ztn[current_level_index] =
                ztn[current_level_index]
                + zebs[current_level_index] * ztn[level_below_index]
                  / below_sigma_kappa;
        }
    }

    /*
     * 6.6 将隐式求得的新温度转换成垂直扩散倾向。
     *
     * zt 是进入 VDIFF 时、已经包含此前物理过程倾向的临时温度；ztn 是完成
     * 位温扩散后对应的新温度。二者之差除以 deltsec2，得到本过程温度倾向，
     * 再累加到模式总温度倾向 dtdt。这里没有再次进行 sigma^kappa 转换，
     * 因为 ztn 和 zt 都已经是温度 T，而不是位温 theta。
     */
    for (level = 0; level < nlev; ++level)
    {
        for (horizontal = 0; horizontal < nhor; ++horizontal)
        {
            const size_t level_index =
                field_index(horizontal, level, nhor);

            zdtdt[level_index] =
                (ztn[level_index] - zt[level_index]) / deltsec2;
            dtdt[level_index] = dtdt[level_index] + zdtdt[level_index];
        }
    }

    /* 释放函数入口处分配的整块临时工作区。 */
    free(work);
}
