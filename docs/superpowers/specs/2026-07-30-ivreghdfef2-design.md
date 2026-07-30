# ivreghdfef2 设计文档

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

新命令 `ivreghdfef2`：

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

自由度算术（对照 `ivreg2.ado` 1140–1200 行）：

| 量 | 原 `ivreghdfe`（`noconstant`, `rankxx=K`, `sdofminus=a`） | 增广 run（含常数, `rankxx=K+1`, `sdofminus=a-1`） |
|---|---|---|
| `small` 无 cluster 的 V 缩放分母 | `N-K-a-dofminus` | `N-(K+1)-(a-1)-dofminus` = 同 |
| `small` + cluster 的 V 缩放分母 | `N-K-a` | `N-(K+1)-(a-1)` = 同 |
| `Fdf1 = rankxx - cons` | `K-0` | `(K+1)-1` = 同 |
| 无 `small` 时 | V 完全不缩放 | 同 |

即 `sdofminus(df_a - 1)` 在所有分支精确抵消 `rankxx` 的 +1。其中 `a = max(e(df_a), 1)`，与 `ivreghdfe.ado` 654 行的 `if (absorb_ct==0) loc absorb_ct 1` 一致。

唯一对不上的是 `e(df_m) = rankxx - cons + (sdofminus - partialcons)`：`partialcons` 是 `ivreg2` 的内部量，外部拿不到，因此增广 run 的 `e(df_m)` 会多出 `a-1`。重新 post 时用参照 run 的值覆盖。

### 方案 B（否决）：手工给 `e(b)`/`e(V)` 追加 `_cons` 行

照抄 `reghdfe.mata` 的分块求逆公式。需要 `(Xhat'Xhat)^{-1}` 与残差 meat 矩阵，`ivreg2` 都不 post，等于要按 robust / cluster / kernel 各重写一遍 VCE，极易与 `ivreg2` 对不上。

### 方案 C（否决）：只修斜率

斜率 SE 取 `ivreghdfe`，`_cons` 那一行沿用手工二阶段的。代码最少，但 `_cons` 的 SE 依然错误，而它会被印进表格。

## 4. 语法

```stata
ivreghdfef2 depvar [exogvars] (endog = instruments) [if] [in] [weight] ,
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
2. **参照回归**。`qui ivreghdfe <原样输入> if touse [w], absorb() <options>`。生成 `esample`，把 `e()` 的全部 scalars / macros / matrices 抓进 tempname（用 `scalar` 保精度，不用 `local =`），取 `e(b)`、`e(V)` 作为斜率的权威参照，取 `a = max(e(df_a), 1)`。
3. **第一阶段**。对每个内生变量 j：`reghdfe endog_j iv exog if esample [w], absorb() <vce> <hdfeopts>` → `est store`。这是标准 OLS，标准误本来就对，且天然带 `_cons`。记录每个（内生变量 × 排除性工具）的系数与 p 值，存成 `e(fsb)` / `e(fsp)` 供汇总区与 replay 使用。
4. **去均值**。对 `depvar + endog + exog + iv` 每个变量 v：
   - `qui reghdfe v if esample [w], absorb() residuals(dm_v) <hdfeopts>`
   - `qui summarize v if esample [aw/fw]` → `replace dm_v = dm_v + r(mean)`（`pweight` 在求均值时按 `aweight` 处理，均值相同）
   - 指定 `noconstant` 时跳过加回均值这一步。
5. **增广二阶段**。`qui ivreg2 dm_y dm_exog (dm_endog = dm_iv) if esample [w], <ivopts> sdofminus(sdofminus + a - 1)`。
6. **校验**。把 `e(b)` 的列名从临时变量映射回真实变量名，**按名字**逐一比对第 5 步与第 2 步的斜率与 SE。相对误差 > 1e-6 打警告，> 1e-4 直接 `error 498` 并打印对照表。绝不静默出错。
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
- 新增：`e(cdf)`、`e(rkf)`、`e(rklm)`、`e(rklmp)`、`e(hansenj)`、`e(hansenjp)`、`e(fsb)`、`e(fsp)`。
- `e(cmd) = "ivreghdfef2"`、`e(cmdline)`。
- `e(b)`、`e(V)` 含 `_cons`。

## 8. 测试方案

`tests/test_ivreghdfef2.do`，数据用 `webuse nlswork`：

1. **斜率与 SE 等于 `ivreghdfe`**：unadjusted / robust / cluster / 双向 cluster / `small` 各一组，逐元素比对，容差 1e-8 相对误差。
2. **`_cons` 等于 `ybar - Xbar'b`**（手工算）。
3. **退化检验（对 `_cons` 的 VCE 唯一的硬检验）**：拿内生变量自己当工具变量 `(x = x)`，此时 2SLS ≡ OLS，`e(b)` 与 `e(V)` 应与 `reghdfe y x controls, absorb() <同 vce>` 逐元素相等。
4. **与旧命令对比**：与 `ivreghdfef` 点估计相同、SE 不同，记录差异幅度作为"修了什么"的证据。
5. **出表**：`esttab` / `reg2docx` 跑通。
6. **边界**：多内生变量、`if` / `in`、三种权重、缺失值下的样本一致性、`absorb()` 留空。

## 9. 刻意划出去的范围（写进 help）

- **变量列表不支持 factor variable**（`i.x`、`i.a#i.b`）。`fvrevar` 展开后基准组与共线组会让列名和设计矩阵错位，可靠对齐需要另写一层解析。检测到 `i.` / `#` 直接报错，提示把分类变量放进 `absorb()`。ts 算子（`L.x`、`D.x`）支持。
- **不支持 `partial()`**。它与 `absorb()` 叠加后 `sdofminus` 的算术需要 `ivreg2` 的内部 `partialcons`，外部取不到。
- **`bw()` / `kernel()` 下的第一阶段**：`reghdfe` 不支持 HAC，第一阶段退化为 robust 并给出提示；第二阶段不受影响。
- **性能**：比 `ivreghdfe` 多 1 次参照回归 + k 次去均值（k = 变量个数），大面板上约 2–3 倍耗时。
- **`_cons` 的含义**：与 `reghdfe` 口径一致，等于固定效应的均值。含固定效应的模型里它不是结构参数，仅供出表。

## 10. 文件清单

```
ivreghdfef2.ado
ivreghdfef2.sthlp
ivreghdfef2.pkg
stata.toc
README.md
LICENSE
tests/test_ivreghdfef2.do
docs/superpowers/specs/2026-07-30-ivreghdfef2-design.md
```
