# Synthetic batch example / 合成批处理示例

Six generated CSV spectra sample both sides of a frequency-order crossing.
`manifest.json` contains explicit model seeds and tracking settings;
`truth.json` stores the independently generated branch centers.
These are synthetic software fixtures, not experimental measurements or a blind peak-discovery benchmark.

六条 CSV 均为程序生成的合成光谱。清单中保留明确的初始峰参数和追踪阈值；
真值文件用于检查频率排序交叉后的轨迹连续性，不代表真实实验结果。

From the repository root / 在仓库根目录运行：

```bash
raman-analyze --lang zh validate examples/batch/manifest.json --check-files
raman-analyze --lang zh batch examples/batch/manifest.json --output runs/batch-example
```

Open `runs/batch-example/report.html` for the bilingual offline overview.
See [中文使用指南](../../docs/usage.zh-CN.md) for adapting the manifest to your own data.
