# 基线、峰型比较与峰数稳定性

[中文首页](../README.zh-CN.md) · [自动选峰](automatic.zh-CN.md) · [实测验证](real-validation.zh-CN.md) · [公式与实现](profiles-and-baselines.md)

## 从哪一步开始

先确认原始谱的横坐标单位、表头和有效测量范围，再比较基线和峰型。程序不读取材料名称来猜峰位、峰数或压力关系。
自动模式默认仍只使用 Gaussian，保证旧清单的行为；要进行三种峰型比较，明确写入 `shapes`：

```json
{
  "data_kind": "measurement",
  "input": {"path": "data/spectrum.xlsx", "header": null, "x": 0, "y": 1},
  "baseline": {"method": "arpls", "lam": 1000000, "iterations": 100, "tolerance": 0.001, "resample": true},
  "fit": {
    "auto": {"shapes": ["gaussian", "lorentzian", "voigt"], "max_peaks": 5, "criterion": "bic", "delta_threshold": 6, "prominence_snr": 4, "multistarts": 3},
    "max_nfev": 3000,
    "residual_correlation_limit": 0.35
  }
}
```

保存为 `config.json` 后运行：

```bash
raman-analyze --lang zh run config.json --output runs/three-shapes
```

以上数值说明配置格式，不能当作所有仪器通用的推荐参数。批处理时把 `baseline` 和 `fit` 放进清单的 `defaults`。

## 三种峰型的参数分别是什么

所有峰型的 `amplitude` 均表示峰顶高度，`center` 表示中心。一个候选模型中的各个组分使用同一种峰型；当前不搜索 Gaussian 和 Lorentzian 混合组分组合。

| 峰型 | 每个组分的拟合参数 | FWHM | 完整实轴上的积分面积 |
|---|---|---|---|
| Gaussian | 高度 A、中心 μ、标准差 σ | `2 sqrt(2 ln 2) σ` | `A σ sqrt(2π)` |
| Lorentzian | 高度 A、中心 μ、半高半宽 γ | `2γ` | `π A γ` |
| Voigt | 高度 A、中心 μ、Gaussian σ、Lorentzian γ | 数值解精确的半高方程 | `A / V(0;σ,γ)` |

Voigt 使用 Normal 与 Cauchy 的卷积，不是两条曲线的线性混合，也不是 pseudo-Voigt。
采用 [SciPy 的 Voigt 定义](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.voigt_profile.html)。
`peaks.csv` 的 `sigma_cm1` 对 Lorentzian 留空，`gamma_cm1` 对 Gaussian 留空；不能将两者解释成同一种宽度。
面积是整个理论峰的积分，窗口截断时不能把它当成窗口内实际测得的积分。
面积与 FWHM 的局部误差传播使用完整的参数协方差交叉项；Voigt 的 σ/γ 可能难以分别辨识。

手动模式使用 `fit.shape`，并填写对应的宽度参数：

```json
{"fit": {"shape": "voigt", "peaks": [{"center": 500, "sigma": 8, "gamma": 5, "sigma_bounds": [0.1, 50], "gamma_bounds": [0.1, 50]}]}}
```

Lorentzian 手动初值使用 `gamma`/`gamma_bounds`；Gaussian 使用 `sigma`/`sigma_bounds`。手动与自动设置不能混用。

## 如何公平比较峰型和峰数

在一个固定基线和窗口内，三种峰型拟合完全相同的原始点；零峰线性背景只比较一次。
Gaussian/Lorentzian 的均值参数数目是 `3N+2`，Voigt 是 `4N+2`；未提供测量误差时，还计入一个估计噪声方差参数。
不会只按 RSS 选更灵活的 Voigt，也不会排除诊断不佳的收敛模型来制造更好看的排序。

自动初始化使用等效 Gaussian 宽度。`min_sigma_cm1`/`max_sigma_cm1` 定义的分辨率区间对应总 FWHM 为这两个值乘 `2 sqrt(2 ln 2)`。
Lorentzian 的 γ 边界按相同 FWHM 区间换算。Voigt 的两种内禀宽度允许接近单一峰型极限，但拟合后总 FWHM 超出该区间仍标为 INCONCLUSIVE。
允许卷积分量接近零，不等于放宽可分辨峰的宽度要求。

查看逐谱 `model-selection.html`：它同时列出各峰型内部的最佳峰数，以及全部“峰型 × 峰数”候选的 BIC、AICc、代价差和尝试次数。
`selected_shape` 是总体最佳峰型，`plausible_models` 保存代价差阈值内的候选组合。
`shape_count_stability:SENSITIVE` 表示各峰型内部的最佳峰数不同，并使最终判断保留 INCONCLUSIVE。
只有一种峰型被请求时，`AGREEMENT` 没有跨峰型证据。

`shape-comparison.csv` 保存每种峰型各自最佳模型的曲线、残差和组分，供复核；不是只保留最终胜出的曲线。
组分列及协方差按输入参数顺序保存，峰参数表按中心排序，两种顺序不一定相同。

## AsLS 和 arPLS 如何选择

AsLS 使用固定的正、负残差权重，适用于明确假设正峰位于背景之上的场景。它可能把对称噪声的下沿当作基线。
arPLS 从负残差的均值和标准差调整 logistic 权重，改善这类噪声偏移；它仍可能吸收宽峰、受负向异常影响，或估计过高的小峰背景。
方法依据 [Baek 等人的 arPLS 论文](https://doi.org/10.1039/C4AN01061B)；权重方程也见 [pybaselines 方法文档](https://pybaselines.readthedocs.io/en/latest/algorithms/whittaker.html#arpls-asymmetrically-reweighted-penalized-least-squares)。实现独立编写，没有增加 pybaselines 依赖。

arPLS 参数：`lam` 默认 10^6，`iterations` 默认 100、范围 1–500，`tolerance` 默认 0.001、范围 (0,1)。
相邻迭代权重的相对变化低于 tolerance 才报告收敛；耗尽预算或残差统计退化时保留基线曲线、原因和 INCONCLUSIVE。
`report.json.baseline_diagnostics` 保存实际迭代数、终止原因及是否使用插值网格。
AsLS 延续固定迭代预算，不声称每次已经收敛。

两者默认要求近似均匀采样。显式 `resample:true` 只对基线估计构造同点数的均匀网格，再把基线映射回原轴。
拟合仍使用原始未平滑点；原始强度、尖峰和负向异常均保留。λ 随采样密度和背景尺度变化，没有自动寻找最高 PASS 率的 λ。

## 跨基线与窗口的稳定性

先按事前确定的基线/窗口设置分别完成批处理，保留每次输出。不要在不同校正后数据间直接比较 BIC/AICc 数值，选择对自己有利的基线。
可以复用保存结果生成稳定性总览：

```bash
raman-analyze --lang zh stability runs/asls/batch-report.json runs/arpls-low/batch-report.json runs/arpls-mid/batch-report.json runs/arpls-high/batch-report.json --output runs/stability
```

所有批次必须有相同的光谱 ID、源文件哈希、列定义与原始点；自动设置、误差处理和拟合预算须一致，仅基线及明确声明的窗口可变。
程序先验证保存结果的文件哈希，随后比较每次总体最佳峰数和每种峰型内部的最佳峰数，输出 `stability.csv`、`stability.json` 与中文离线页面。

| 稳定性字段 | 含义 |
|---|---|
| `SUPPORTED_STABLE` | 各变体及各峰型最佳峰数一致，全部相关诊断通过，搜索完整且竞争峰数唯一 |
| `AGREEMENT_ONLY` | 峰数相同，但诊断、边界、竞争峰数或搜索仍有问题 |
| `SENSITIVE` | 至少一组设置给出不同的最佳峰数 |
| `UNRESOLVED` | 存在 BLOCKED 记录，不能完整比较 |

稳定性逐条光谱判断，不要求不同压力下的峰数保持相同。它只覆盖提供的设置，不能证明在所有合理参数下稳定，更不等于真实物理峰数的准确率。
局部误差仍不包含基线、校准和峰型/峰数选择不确定性。

## 独立核对保存结果

```bash
python tools/check_saved_batch.py runs/arpls-mid/batch-report.json --output runs/arpls-mid-check.json
```

检查器独立从保存的点重新计算每种峰型最佳模型的 RSS、参数惩罚与 BIC/AICc，用分布函数重建理论曲线，并核对输出哈希和恒等式。
这检验数据与软件输出的一致性，不替代真实峰数的独立参考标注。


下图是公开合成 Voigt 单峰的对照，用于展示错误峰型如何用多个组分吸收失配；它不代表实测性能。

![Synthetic profile comparison](synthetic-profile-comparison.png)
