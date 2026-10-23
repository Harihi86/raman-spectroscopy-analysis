# Automatic selection / 自动选峰

This manifest uses the existing synthetic crossing CSV files but supplies no peak positions, widths or counts to the fitter.

这份清单引用公开的合成交叉峰数据；拟合器不接收生成器的峰位、宽度或峰数。

```bash
raman-analyze --lang zh batch examples/automatic/manifest.json --output runs/automatic
```

Read each spectrum's `model-selection.csv` and `model-selection.json`, alongside the fitted curve and residuals.
See [自动选峰说明](../../docs/automatic.zh-CN.md). Synthetic truth is available in `../batch/truth.json`.
