# 自动选峰与峰数判断

[中文首页](../README.zh-CN.md) · [完整使用指南](usage.zh-CN.md) · [方法公式与参考](automatic.md) · [实测验证](real-validation.zh-CN.md)

## 使用方式

在单谱配置或批处理的 `defaults.fit` 中，使用 `auto` 替换 `peaks`：

```json
{
  "input": {"path": "data/spectrum.xlsx", "header": null, "x": 0, "y": 1},
  "baseline": {"method": "asls", "lam": 1000000, "p": 0.01, "iterations": 30, "resample": true},
  "fit": {
    "auto": {"max_peaks": 5, "criterion": "bic", "delta_threshold": 6, "multistarts": 3},
    "max_nfev": 3000
  },
  "data_kind": "measurement"
}
```

这是格式示例。文件必须使用 cm^-1 横坐标；基线、宽度范围和拟合窗口需要结合数据选择。
波长数据需要明确的物理换算，不能直接标成 cm^-1。

```bash
raman-analyze --lang zh run config.json --output runs/auto-single
raman-analyze --lang zh validate manifest.json --check-files
raman-analyze --lang zh batch manifest.json --output runs/auto-batch
raman-analyze --lang zh demo-batch --automatic --scenario crossing --output runs/auto-demo
```

`fit.auto` 和 `fit.peaks` 不能同时存在。逐谱切换模式时，会清除共享配置中另一种模式的设置；`fit.window` 等共同参数仍可继承。

## 程序实际做什么

1. 使用声明的基线和拟合范围。零峰模型是线性残余背景，也是必须比较的候选。
2. 由一阶差分的 MAD 估计噪声尺度；它只是启发式估计，不会冒充用户提供的测量误差。
3. 在用于提出初值的平滑信号上寻找具有突出度的局部峰，并把宽度转换到实际横坐标单位。
4. 对 1–N 个正 Gaussian 组分分别进行有界拟合；初值来自检测、前一模型未解释的残差或宽组分拆分。默认最多保留 32 个检测提议，每个峰数进行 3 次局部初始化尝试。
5. 所有模型都拟合相同的原始采样点，并含有线性残余背景。平滑谱不参与最终拟合，不增加样本量。
6. 对成功收敛的候选比较 BIC 或 AICc。保存代价、竞争峰数、初值、每次优化的状态和拟合诊断。

这是一种有资源边界的局部搜索，不保证找到每个峰数的全局最优解。`search_complete:true` 仅表示请求的每个峰数至少得到一个收敛结果；精确平坦/线性输入会明确跳过无峰对比度的正峰模型。

## 参数

| 字段 | 默认值 | 含义 |
|---|---:|---|
| `max_peaks` | 5 | 比较 0 到此上限；范围 1–8 |
| `criterion` | `bic` | `bic` 或 `aicc`，越小越优 |
| `delta_threshold` | 6 | 距最小代价在此范围内的峰数均保留；不是置信概率 |
| `prominence_snr` | 4 | 检测突出度相对噪声尺度的阈值，同时检查选中组分的幅度 |
| `smooth_sigma_points` | 1.5 | 提议阶段的 Gaussian 平滑尺度，单位为采样点；0 关闭 |
| `min_sigma_cm1` | 中位横坐标间距 | 组分 sigma 下限；可明确指定正值 |
| `max_sigma_cm1` | 拟合跨度 / 3 | 组分 sigma 上限；可明确指定正值 |
| `min_separation_cm1` | 中位间距的两倍 | 检测与选中组分的最小间距检查，不能替代仪器分辨率 |
| `multistarts` | 3 | 每个峰数的初始化尝试次数，1–12 |

`fit.max_nfev` 控制每次自动拟合的计算上限，默认 3000。`fit.residual_correlation_limit` 延续单谱诊断，默认 0.35；关闭告警不会使模型获得额外的物理证据。

## 读取结果

先打开每条谱的 `model-selection.html`，直接比较各峰数的代价和诊断；批处理总览中也有链接。

- `report.json` 的 `model_selection.selected_peak_count` 是当前代价最小的候选。
- `plausible_peak_counts` 保存代价差在阈值内的候选峰数，不是概率区间。
- `count_status` 同时考虑选中模型的拟合诊断和峰数判断问题。为 INCONCLUSIVE 时，不能把最佳候选当作已确认峰数。
- `model-selection.csv` 方便比较各峰数的 BIC、AICc、参数数目、代价差和拟合状态。
- `model-selection.json` 还包含提议、初始化参数和每次优化的诊断。
- `fit.csv` 和拟合图用于检查原始点、模型分量与残差。零峰结果仍输出可读取的 `peaks.csv` 表头。
- 批处理的离线总览显示“最佳候选 / 竞争峰数”，汇总表保留峰数状态。非 PASS 拟合不会更新压力轨迹。

达到峰数上限、竞争模型接近、部分峰数未收敛、局部协方差不可靠、参数触界、低幅度组分、窗口截断或相关残差，都会使结果保留不确定性。

## 在真实数据中必须复核什么

荧光背景、仪器尖峰、负向异常、非 Gaussian 线型和重叠峰都可能被额外组分吸收。较低的 BIC 不证明发现了更多物理振动模式。
建议比较合理的基线、窗口和宽度范围，检查判断是否稳定，并保留所有运行。

AsLS 默认要求近似均匀采样。显式 `resample:true` 只在基线估计阶段构造同样点数的均匀网格，再映射回原始坐标。原始谱的拟合和导出不被重采样；基线不确定性仍不包含在局部协方差内。

目前已用真实压力光谱开展本地适用性检查；没有独立人工或仪器参考峰数标签，因此不能报告真实峰数准确率。
合成真值检查、真实数据可运行、模型诊断通过、物理峰数正确是不同层次的证据。
