# ivreghdfe2 设计文档

日期：2026-07-30

## 1. 背景与问题

`ivreghdfef`（本作者早前的命令）用手工两阶段实现 IV：

```stata
reghdfe x  iv control, absorb(...)      // 第一阶段
predict xhat
reghdfe y  xhat control, absorb(...)    // 第二阶段
```

第二阶段的**点估计是正确的**（含 `_cons`），因为 2SLS 的系数确实等于用 `xhat` 做 OLS 的系数。但**方差矩阵是错的**：`reghdfe` 用的残差是 `y - xhat*b`，而 2SLS 的正确残差是 `y - x*b`。

- 同方差情形下两者只差一个标量 `s2` 的比值，尚可缩放修正；
- **cluster / robust 情形下 meat 矩阵结构不同，不存在标量修正**。

而现有工具都不能直接给出想要的输出：

- `ivreghdfe` 的 VCE 是对的，但在 `absorb()` 存在时会强制 `noconstant`（`ivreghdfe.ado` 674–681 行），**不报常数项**，也不保存可用于出表的两阶段估计。
- `reghdfe` 的 `_cons` 是用"回归过均值 + 分块矩阵求逆"加回去的（`reghdfe.mata` 3830–3870 行），口径明确，但它不是 IV 估计。

## 2. 目标

新命令 `ivreghdfe2`：

1. 语法贴近 `ivreghdfe`（括号语法 + 选项直写），不再用 `ivreghdfef` 的 `iv()` / `opt()` / `regopt()` 包装。
2. 输出 `ivreghdfef` 那套额外统计量（CD-Wald F、KP rk Wald F、KP rk LM、Hansen J）。
3. **报告常数项**，并把**两个阶段**都保存成可 `esttab` / `reg2docx` 的估计结果。
4. 第二阶段的标准误必须是正确的 2SLS 标准误，`_cons` 那一行也必须正确。

## 3. 方案选择

### 方案 A（采用）：去均值 + `ivreg2` 增广回归

用 `reghdfe v, absorb() residuals()` 把每个变量的固定效应 partial out，再把各自的样本均值加回去，然后交给 `ivreg2` 跑**带常数项**的 2SLS，传入 `sdofminus(df_a - 1)` 对齐自由度。

正确性依据：

- **斜率与残差不变**。由 FWL 定理，对 `y~ = M_D y + ybar` 关于 `[X~, 1]`（`X~ = M_D X + Xbar`）用 `[Z~, 1]` 做 2SLS，工具变量矩阵张成的列空间等于 `span{M_D Z, 1}`，因此斜率与残差和 `ivreghdfe`（`M_D y` 关于 `M_D X`，无常数）逐位相同。
- **常数项口径与 `reghdfe` 一致**：`_cons = ybar - Xbar'b`，正是 `reghdfe_extend_b_and_xx()` 用的"回归过均值"公式。
- **VCE 全部由 `ivreg2` 计算**：robust / 单向与双向 cluster / HAC kernel / GMM / LIML 全部走它自己的成熟代码，包括 `_cons` 那一行，无需我们重写任何方差公式。

#### 自由度的处理：一条走错的路，和最终的做法

最初的设计是给增广 run 传 `sdofminus(df_a - 1)`，用来抵消 `rankxx` 因为多了常数项而 +1。**这是错的**，理由在 `ivreghdfe.ado` 644–700 行：

```
loc partialcons `absorb_ct'                              // 656 行，absorb_ct 通常远大于 1
if ("`partial1'" != "" | `partialcons'==1) { ... }       // 659 行
else { local partial_ct 0 ; local partialcons 0 }        // 688-691 行
if "`nopartialsmall'"=="" { local sdofminus = `sdofminus'+`partial_ct' }   // 694-695 行
```

`partialcons` 被赋成 `absorb_ct`（固定效应个数，例）而不是 0/1 布尔量，所以 `partialcons==1` 这个条件在正常情况下**根本不成立**，走的是 else 分支：`partial_ct = 0`，`sdofminus` 保持用户传入的值不变。也就是说固定效应个数根本没有进入 `sdofminus`。任何基于"`sdofminus` 里含 `absorb_ct`"的推导都是空中楼阁。

最终做法**完全绕开自由度约定**：

1. 增广 run 不传 `small`。由 `ivreg2.ado` 1180–1192 行，不传 `small` 时 V 完全不缩放。
2. 从参照 run 反解出那个标量缩放因子 `f = V_ref[i,i] / V_aug[i,i]`。
3. 因为 `small` 修正是乘在整个 V 上的**纯标量**（1184 行与 1186–1187 行都是 `matrix V = V * <标量>`），`f` 对每一个斜率必须完全相同。命令逐个算出 `f`，取 `(min+max)/2`，并检查离散度；离散度超阈值就 `error 498`。
4. 把 `f` 乘到整个增广 V 上（含 `_cons` 行），再把斜率块用参照值逐元素覆盖。

这个做法同时是**结构正确性的硬校验**：如果增广 run 的残差、cluster 划分或设计矩阵与参照 run 有任何实质差异，`f` 就不会是常数，命令会报错而不是静默输出错结果。校验结果留在 `e(chk_b)`、`e(chk_v)`、`e(smallfac)` 里，可从外部断言。

`_cons` 行采用参照 run 的标量 `f`（而不是"增广模型自己该用的" `rankxx=K+1` 版本），是刻意的选择：这样斜率与 `ivreghdfe` 逐位相同，而 `_cons` 与斜率共享同一套自由度约定，内部一致。这也正是 `reghdfe` 的口径（`df_m = K`，常数项不计入，见 `reghdfe.mata` 3625 行 `sol.df_m = sol.rank = sol.K`）。

`e(df_m)`、`e(df_r)` 等标量一律用参照 run 的值覆盖。

### 方案 B（否决）：手工给 `e(b)`/`e(V)` 追加 `_cons` 行

照抄 `reghdfe.mata` 的分块求逆公式。需要 `(Xhat'Xhat)^{-1}` 与残差 meat 矩阵，`ivreg2` 都不 post，等于要按 robust / cluster / kernel 各重写一遍 VCE，极易与 `ivreg2` 对不上。

### 方案 C（否决）：只修斜率

斜率 SE 取 `ivreghdfe`，`_cons` 那一行沿用手工二阶段的。代码最少，但 `_cons` 的 SE 依然错误，而它会被印进表格。

## 4. 语法

```stata
ivreghdfe2 depvar [exogvars] (endog = instruments) [if] [in] [weight] ,
    Absorb(absvars) [ FIRST STOre(name) Format(%fmt) NOSUMmary  其他 ivreghdfe 选项 ]
```

- 括号前后都可以放外生控制变量，两处都会被并入 `inexog`。
- `absorb()` 可选。留空时跳过去均值与增广回归，直接使用参照 run 的结果（此时 `ivreghdfe` 本来就报常数项）。
- `weight` 支持 `aweight` / `fweight` / `pweight`。不支持 `iweight`（`reghdfe` 不接受）。

自有选项：

| 选项 | 作用 |
|---|---|
| `first` | 额外显示第一阶段回归表 |
| `store(name)` | 保存 `name_second`；第一阶段为 `name_first`（单内生）或 `name_first1`、`name_first2`…（多内生） |
| `format(%fmt)` | 汇总区数字格式，默认 `%10.4f` |
| `nosummary` | 不显示统计量汇总区 |

其余选项（`cluster()`、`robust`、`small`、`gmm2s`、`liml`、`bw()`、`kernel()`、`dofminus()`、`ffirst`…）透传给 `ivreghdfe` 与 `ivreg2`。`vce(robust)` / `vce(cluster ...)` 会被翻译成 `robust` / `cluster()`。

## 5. 执行流程

1. **解析**。按 `(`、`=`、`)` 把命令行切成 depvar / exog / endog / iv 四组。检测并拒绝 factor variable；ts 算子经 `tsrevar` 解析成临时变量（1:1 映射）。`marksample touse`。
2. **参照回归**。`qui ivreghdfe <原样输入> if touse [w], absorb() <options> nooutput`。生成 `esample`，把 `e()` 的全部 scalars / macros / matrices 抓进 tempname（用 `scalar` 保精度，不用 `local =`，后者会截断有效数字），取 `e(b)`、`e(V)` 作为斜率的权威参照。
3. **第一阶段**。对每个内生变量 j：`reghdfe endog_j iv exog if esample [w], absorb() <vce> <hdfeopts>` → `est store`。这是标准 OLS，标准误本来就对，且天然带 `_cons`。记录每个（内生变量 × 排除性工具）的系数与 p 值，存成 `e(fsb)` / `e(fsp)` 供汇总区与 replay 使用。
   - **注意 `reghdfe` 与 `ivreg2` 的 vce 拼写不同**：`reghdfe` 只有 `vce()` 和 `cluster()`，**没有裸 `robust`**（`reghdfe.ado` 103 行起的 syntax）。所以要维护两套写法：给 `ivreghdfe` / `ivreg2` 的 `robust` / `cluster()`，和给 `reghdfe` 的 `vce(robust)` / `vce(cluster ...)`。
4. **去均值**。对 `depvar + endog + exog + iv` 每个变量 v：
   - `qui reghdfe v if esample [w], absorb() residuals(dm_v) <hdfeopts>`
   - `qui summarize v if esample [aw/fw]` → `replace dm_v = dm_v + r(mean)`（`pweight` 在求均值时按 `aweight` 处理，均值相同）
   - 指定 `noconstant` 时跳过加回均值这一步。
5. **增广二阶段**。`qui ivreg2 dm_y dm_exog (dm_endog = dm_iv) if esample [w], <ivopts 去掉 small> nooutput`。
   - `nooutput` 不只是为了输出干净：当被 absorb 的固定效应嵌套在 cluster 变量里时，残差在每个 cluster 内精确求和为零，于是新加的常数项那一列的 cluster score 恒为 0，`ivreg2` 会警告 moment covariance 矩阵秩亏（`ivreg2.ado` 1793 行判定、2577 行打印）。这条警告在这里是预期且无关的——所有对外报告的统计量都来自参照 run，警告只涉及增广 run 自己那个被丢弃的 J 统计量。`quietly` 拦不住它（`di in r`），`nooutput` 会整块跳过显示。
6. **校验与标量缩放**。把 `e(b)` 的列名从临时变量映射回真实变量名，**按名字**逐一比对第 5 步与第 2 步的斜率，同时算出每个斜率的方差比 `f`。斜率相对误差 > 1e-5 直接 `error 498`；`f` 的离散度 > 1e-4 报错、> 1e-6 提示。然后 `V = f * V`（含 `_cons` 行），斜率块用参照值逐元素覆盖。绝不静默出错。
   - 阈值定在这里是量出来的：`f` 的离散度是 meat 矩阵条件数带来的浮点噪声，随 cluster 数减少而变大（3410 个 cluster 约 5e-11，12 个约 6.5e-9），且**与 `tolerance()` 无关**（1e-10 / 1e-12 / 1e-13 都是 6.503e-09）。结构性错误会给出 0.1–1 量级的离散度，所以 1e-4 的报错线留了充足余量。
7. **Post**。`ereturn post b V, esample() depname() obs()`；`small` 时额外 `dof(df_r)`（`ereturn display` 靠 `e(df_r)` 决定显示 t 还是 z，这与 `ivreg2` 的 `small` 语义一致）。回填第 2 步抓下来的全部 `e()`，再补 `e(cdf)` `e(rkf)` `e(rklm)` `e(rklmp)` `e(hansenj)` `e(hansenjp)`（与 `ivreghdfef` 同名，`reg2docx ... scalars()` 可直接复用）。`est store name_second`。
8. **显示**。紧凑表头（N、F、R²、absorbed FE）+ `ereturn display`（含 `_cons`）+ 汇总区。带 `first` 时先打第一阶段表。

## 6. 汇总区格式

```
=================================    Results Summary   =================================
==> First stage coefficient
beta for ttl_exp <- hours:  0.1234***
==> Second stage coefficient
beta for ln_wage <- ttl_exp_hat:  0.5678**
==> Weak-identification
Cragg Donald Wald F Test:  123.4567
Kleibergen Paap rk Wald F Test:  98.7654
==> Under-identification
Kleibergen Paap rk LM Test:  45.6789***
==> Over-identification
Hansen J Test:  1.2345
```

星号：`*** p<0.01`、`** p<0.05`、`* p<0.1`。多内生变量时第一阶段那段按（内生变量 × 工具）逐行展开。

## 7. e() 返回

- 沿用参照 run 的全部 scalars / macros / matrices（`N`、`df_m`、`rmse`、`rss`、`r2`、`r2_a`、`j`、`jdf`、`jp`、`idstat`、`idp`、`cdf`、`rkf`、`widstat`、`df_a`、`N_hdfe`、`N_clust`、`clustvar`、`depvar`、`instd`、`insts`、`exexog`、`inexog`、`absvar` 等）。
- 新增：`e(rklm)`、`e(rklmp)`、`e(hansenj)`、`e(hansenjp)`（`e(cdf)`、`e(rkf)` 参照 run 本来就有）、`e(smallfac)`、`e(chk_b)`、`e(chk_v)`、`e(df_r_ref)`、`e(fsb)`、`e(fsp)`。
- `e(cmd) = "ivreghdfe2"`、`e(cmdline)`、`e(endovars)`、`e(exexvars)`、`e(inexogvars)`、`e(absvarlist)`、`e(fsstored)`、`e(ssstored)`。
- `e(b)`、`e(V)` 含 `_cons`。
- `e(rklmp)` / `e(hansenjp)` 直接取 `ivreg2` 的 `e(idp)` / `e(jp)`，不自己重算。旧 `ivreghdfef` 用 `chi2tail(e(exexog_ct), ...)`，在单内生变量下恰好等于正确的 `L-K+1` 与 `L-K`，但多内生变量时会错。

## 8. 测试方案

`tests/test_ivreghdfe2.do`，数据用 `webuse nlswork`，共 68 项断言。

1. **对齐 `ivreghdfe`**：unadjusted / robust / cluster / 双向 cluster，各带与不带 `small`，共 9 组。
   - 注意：斜率块是从参照 run 复制过来的，所以直接比对 `e(b)` / SE 必然相等，**这一项只验证管道（列名映射、顺序、样本）**。真正的数值校验靠命令内部暴露出来的 `e(chk_b)`（原始增广估计与参照的最大相对系数差）与 `e(chk_v)`（方差比的离散度）。断言 `chk_b < 1e-10`、`chk_v < 1e-7`。
2. **`_cons` 等于 `ybar - Xbar'b`**（手工算）。
3. **退化检验（对 `_cons` 的 VCE 唯一的硬检验）**：拿内生变量的副本当工具变量，此时 2SLS ≡ OLS，`e(b)` 与 `e(V)` 应与 `reghdfe y x controls, absorb() <同 vce>` 逐元素相等。注意不能写 `(x = x)`，`ivreg2` 不接受内生变量出现在自己的工具变量列表里，要先 `gen` 一个副本。
   - **嵌套 cluster 单列一项**：cluster 变量同时被 absorb 时，`reghdfe` 会多减一个自由度（`reghdfe.mata` 4026 行的 `nested_adj`）而 `ivreghdfe` 不减，两者差一个 1/N 量级的标量。断言"**所有系数（含 `_cons`）的 SE 比值是同一个标量**"——这仍然是对 `_cons` 方差的真实检验：如果 `_cons` 那一行算坏了，它就不会跟着斜率走。
4. **与旧命令对比**：与 `ivreghdfef` 点估计相同、SE 不同，记录差异幅度作为"修了什么"的证据。旧命令会 `drop` 再重建内生变量，测试里要 `preserve` / `restore` 包住。
5. **出表与 replay**：`esttab` / `reg2docx` 跑通、`ivreghdfe2` 空调用可重放。
6. **边界**：多内生变量、`if` / `in`、三种权重、ts 算子、`absorb()` 留空、factor variable 与 `partial()` 被拒。

### 实测结果（68/68 通过）

| 检验 | 结果 |
|---|---|
| `chk_b`（系数） | 4.0e-13（全部 9 组） |
| `chk_v`（方差比离散度） | 1.1e-11 ～ 6.6e-09（随 cluster 数减少而变大） |
| `_cons` = `ybar - Xbar'b` | 6.4e-12 |
| 退化检验 SE vs `reghdfe`：unadjusted / robust / cluster | 4.3e-13 / 1.7e-12 / 1.5e-12 |
| 嵌套 cluster 的公共标量 | 六个系数含 `_cons` 一致到 8.6e-12（标量 = 1.00002808633） |
| 旧命令低估 `ttl_exp` 的 SE | 10.68%（0.03870 vs 0.04333，p 从 0.006 变 0.014） |

## 9. 刻意划出去的范围（写进 help）

- **变量列表不支持 factor variable**（`i.x`、`i.a#i.b`、`c.x`）。`fvrevar` 展开后基准组与共线组会让列名和设计矩阵错位，可靠对齐需要另写一层解析。检测到 `i.` / `#` / `c.` 直接报错，提示把分类变量放进 `absorb()`。ts 算子（`L.x`、`D.x`）支持，因为它是 1:1 映射。
- **不支持 `partial()`**。它与 `absorb()` 叠加后的小样本自由度算术依赖 `ivreg2` 的内部 `partialcons`，外部取不到。
- **`bw()` / `kernel()` 下的第一阶段**：`reghdfe` 不支持 HAC，第一阶段退化为 `vce(robust)` 并给出提示；第二阶段不受影响。
- **不支持 `iweight`**：`reghdfe` 的 syntax 只接受 `fw aw pw`。
- **t 还是 z**：跟随 `ivreg2` 的语义——只有指定 `small` 时才报 t。实现上通过 `ereturn post ... dof()` 是否传入来控制（`ereturn display` 靠 `e(df_r)` 存在与否决定）。
- **性能**：比 `ivreghdfe` 多 1 次参照回归 + k 次去均值（k = 变量个数），大面板上约 2–3 倍耗时。曾考虑用 `hdfe` 一次 partial 掉所有变量，但本机未安装、且会新增依赖，故用 `reghdfe ..., residuals()` 逐个跑。
- **`_cons` 的含义**：与 `reghdfe` 口径一致，等于固定效应的均值。含固定效应的模型里它不是结构参数，仅供出表。

## 10. 文件清单

```
ivreghdfe2.ado
ivreghdfe2.sthlp
ivreghdfe2.pkg
stata.toc
README.md
LICENSE
tests/test_ivreghdfe2.do
docs/superpowers/specs/2026-07-30-ivreghdfe2-design.md
```
