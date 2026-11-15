# ivreghdfe2

[English](README.md) | 简体中文

在 `ivreghdfe` 的基础上修正二阶段方差矩阵，报告常数项，补充工具变量诊断，并把两个阶段都保存为估计结果。

## 用途

使用 `absorb()` 时，`ivreghdfe` 不报告常数项，也不会把两个阶段都保存为估计结果。常见的替代做法是手工跑两阶段最小二乘：

```stata
reghdfe x  iv controls, absorb(...)   // 一阶段
predict xhat
reghdfe y  xhat controls, absorb(...) // 二阶段
```

这样得到的系数是对的：2SLS 的斜率确实等于对拟合值做 OLS 回归的结果。**但报告的标准误是错的。** 手工二阶段用 `y - xhat*b` 构造残差，而 2SLS 要求的是 `y - x*b`。同方差假设下两者只差一个标量；在 `robust` 或 `cluster()` 下，三明治估计量中间那一项的结构本身就不同，无法靠缩放修正。

在下面的示例中，手工做法把内生解释变量的标准误低估了 **10.7%**，*p* 值从 0.006 变成 0.014。

`ivreghdfe2` 的斜率和方差矩阵直接取自 `ivreghdfe`，因此在构造上就是正确的；常数项另行推导，不改动这两者。

## 安装

```stata
* 最新版本
cap ado uninstall ivreghdfe2
net install ivreghdfe2, from("https://raw.githubusercontent.com/codefoxs/ivreghdfe2/main/") replace
```

依赖 [`ivreghdfe`](https://github.com/sergiocorreia/ivreghdfe)、[`ivreg2`](http://fmwww.bc.edu/RePEc/bocode/i/) 和 [`reghdfe`](https://github.com/sergiocorreia/reghdfe)。

## 语法

模型写法与 `ivreghdfe` 完全相同：

```stata
ivreghdfe2 depvar [exogvars] (endog = instruments) [if] [in] [weight] ,
    absorb(absvars) [first store(name) format(%fmt) nosummary  ...ivreghdfe options]
```

`cluster()`、`robust`、`small`、`gmm2s`、`liml`、`bw()`、`kernel()`、`tolerance()` 等选项直接写，原样传给 `ivreghdfe`。

统计推断与 `ivreghdfe` 一致。使用 `absorb()` 时，`ivreghdfe` 总会做小样本调整，所以无论是否写 `small`，二阶段都按 `e(df_r)` 报告 t 统计量。不使用 `absorb()` 时，它与 `ivreg2` 相同：写 `small` 报告 t，否则报告 z。

## 示例

```stata
webuse nlswork

ivreghdfe2 ln_wage age tenure not_smsa south (ttl_exp = hours union), ///
    absorb(idcode year) cluster(idcode) first store(IV)

esttab IV_first IV_second, se scalars(N r2_a cdf rkf rklm hansenj)
```

`store(IV)` 会留下 `IV_first` 和 `IV_second` 两个估计结果，都带有 `_cons` 行。有多个内生解释变量时，一阶段依次保存为 `IV_first1`、`IV_first2` 等。

## 实现方法

`absorb()` 会让 `ivreghdfe` 把常数项也一并剔除（`ivreghdfe.ado` 第 674–681 行强制 `noconstant`），所以本来就没有 `_cons` 行可报告。`ivreghdfe2` 按以下步骤恢复常数项：

1. **参考运行。** 按用户写法运行 `ivreghdfe`。其 `e(sample)`、诊断统计量、斜率和方差矩阵是其余所有结果的基准。
2. **一阶段。** 对每个内生解释变量，用 `reghdfe` 对工具变量和控制变量回归。这是普通 OLS，标准误无需修正，并自带常数项。
3. **剔除固定效应。** 每个变量先对吸收的固定效应取残差，再加回样本均值。加回均值使截距重新可识别；根据 Frisch–Waugh–Lovell 定理，斜率和残差都不受影响。这与 `reghdfe` 采用的“过均值回归”约定相同（`reghdfe.mata` 第 3830–3870 行）。
4. **增广二阶段。** 对变换后的变量运行*带*常数项的 `ivreg2`。因此各种方差设定（`robust`、单向和双向 `cluster()`、HAC 核、GMM）包括 `_cons` 行在内，全部由 `ivreg2` 自身的代码计算，没有重新实现任何方差公式。
5. **校准与核验。** 增广运行*不带* `small`，`ivreg2` 不对方差矩阵做缩放；小样本调整因子由参考运行与增广运行的方差之比得到。`ivreg2` 的小样本调整是给整个矩阵乘一个数，所以这个比值对每个系数都必须相同；命令会检查这一点，不满足就报错退出。这样既不必复现 `ivreghdfe` 内部的 `sdofminus` 计算，又顺带核验了增广运行确实复现了参考拟合。核验结果保存在 `e(chk_b)`、`e(chk_v)` 和 `e(smallfac)`。
6. **报告。** 斜率部分原样复制自参考运行，所以报告的斜率及其标准误与 `ivreghdfe` 逐位相同，只有 `_cons` 行是新增的。

## 正确性证据

`tests/test_ivreghdfe2.do` 在 `nlswork` 数据上运行 92 项检查（92/92 通过），除 `ivreghdfe`、`ivreg2` 和 `reghdfe` 外不依赖其他命令。其中关键的几项：

- **与 `ivreghdfe` 比对**，覆盖九种方差设定（未调整、`robust`、单向和双向 `cluster()`，各自带或不带 `small`）：命令内部把*未经替换*的增广估计与参考运行比较，系数差 `e(chk_b)` 约为 4e-13，方差比值的一致程度在 1e-11 到 7e-9 之间。`e(df_r)` 和二阶段 p 值也与 `ivreghdfe` 一致，所以无论是否写 `small`，都使用 t 推断。
- **常数项的数值**等于 `ybar - Xbar'b`，误差 6e-12。
- **常数项的方差。** 用内生解释变量本身的副本作工具变量时，2SLS 在数值上等同于 OLS，因此完整的 `e(b)` 和 `e(V)`（包括 `_cons` 行）都必须等于 `reghdfe` 的结果。实测误差为 4e-13（未调整）、1.7e-12（`robust`）和 1.5e-12（`cluster()`）。
- **手工两阶段的问题**：用 `reghdfe`、`predict`、`reghdfe` 重现手工做法，其内生解释变量和所有控制变量的系数与 `ivreghdfe2` 相差 3e-15，而内生解释变量的标准误小了 10.7%。
- **嵌套聚类。** 按同时被吸收的面板个体聚类时，`reghdfe` 和 `ivreghdfe` 确实会出现差异：`reghdfe` 多扣一个自由度，两者相差一个量级为 1/N 的标量。全部六个系数（含 `_cons`）共享同一个标量（1.00002808633 对 1.00002808633），12 位有效数字一致，说明 `_cons` 行与斜率同步，而不是单独计算的。
- 此外还检查了权重（`aweight`、`fweight`、`pweight`）、`if`/`in`、时间序列算子、多个内生解释变量、不使用 `absorb()`、`esttab` 和 `reg2docx` 导出，以及结果重显（replay）。

## 局限

- **变量列表中不支持因子变量**（`i.x`、`i.a#i.b`、`c.x`）。要在增广回归中让系数名与设计矩阵保持对应，需要追踪基准组和被省略的水平。类别变量请放进 `absorb()`。支持时间序列算子（`L.x`、`D.x`）。
- **不支持 `partial()`。** 它与 `absorb()` 同用时，小样本自由度的计算依赖 `ivreg2` 内部的一个量。
- **比 `ivreghdfe` 慢**：一次参考运行，加上每个变量一次剔除固定效应，在大型面板上耗时约为 2–3 倍。
- **使用 `bw()`/`kernel()` 时，一阶段**退回 `vce(robust)`，因为 `reghdfe` 不能计算 HAC 标准误。二阶段不受影响。
- **不支持 `iweight`**，因为 `reghdfe` 不接受。
- `_cons` 沿用 `reghdfe` 的约定（吸收的固定效应的均值）。有固定效应时它不是结构参数，报告它只是因为表格需要这一行。

## 参考

- https://github.com/sergiocorreia/ivreghdfe
- https://github.com/sergiocorreia/reghdfe
- http://fmwww.bc.edu/RePEc/bocode/i/ivreg2.html

## 更新记录

### 1.0.0（2026-11-15）

- **修复：二阶段在 `ivreghdfe` 报告 t 的情形下报告了 z 统计量。** 使用 `absorb()` 时，`ivreghdfe` 强制启用 `small`：方差矩阵带有小样本调整，并写入 `e(df_r)`，因此系数表显示 t 统计量。0.1.0 版照搬了这个方差矩阵，却只在用户写了 `small` 时才传入自由度，结果不写 `small` 时，二阶段对一个已经做过小样本调整的方差矩阵报告了 z 统计量和正态分布 p 值。现在 `ivreghdfe2` 直接采用参考运行的 `e(df_r)`。
- **影响范围：** 系数和标准误不变。二阶段的 p 值、置信区间和显著性星号改为按 `e(df_r)` 自由度的 t 分布计算，与 `ivreghdfe` 一致。聚类数多时差别可以忽略，聚类数少时差别明显：12 个聚类时 `e(df_r)` = 11，z 检验的 p 值明显偏小。用 0.1.0 版在 `absorb()` 且未写 `small` 的情形下得到的结果，建议重新运行。
- 不使用 `absorb()` 时行为不变：与 `ivreg2` 相同，默认报告 z，写 `small` 报告 t。
- 测试新增了在全部九种方差设定下与 `ivreghdfe` 比对 `e(df_r)` 和二阶段 p 值的检查（92/92 通过）。

### 0.1.0（2026-07-30）

- 首次发布。
