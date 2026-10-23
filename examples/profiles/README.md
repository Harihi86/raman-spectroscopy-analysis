# Synthetic family-comparison example / 合成峰型比较

```bash
raman-analyze --lang zh batch examples/profiles/manifest.json --output runs/profiles
```

Open `runs/profiles/report.html` and each spectrum's `model-selection.html`.
Each input contains one generated peak, using a different true family. Seeds and
counts are selected from the observations; `truth.json` is not read by the fitter.
The baseline is deliberately fixed to `none` to isolate the profile comparison.
A mismatched family can use extra components; disagreement remains INCONCLUSIVE,
so an exit code of 1 can be an expected completed result.

每条合成谱包含一个峰；三条谱分别由 Gaussian、Lorentzian、Voigt 生成。
用于验证软件的真值与实测验证分开保存。`shape-comparison.csv` 保留所有峰型的最佳曲线。
基线和峰型的完整设置见[中文说明](../../docs/profiles-and-baselines.zh-CN.md)；
真实 184 条光谱的结果另见[实测记录](../../docs/real-validation.zh-CN.md)，不会从这些合成示例推断真实准确率。
