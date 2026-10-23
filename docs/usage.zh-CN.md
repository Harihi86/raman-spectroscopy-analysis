# 中文使用指南：从文件到批处理和压力峰追踪

[返回中文首页](../README.zh-CN.md) · [常见问题](troubleshooting.zh-CN.md) · [参数参考](configuration.md)

## 1. 先确定光谱模型

每条光谱至少需要横坐标与强度。横坐标单位是 Raman shift 的 `cm^-1`。
如果输入为波长，请依据激发波长完成明确的物理换算；当前导入器不会猜测单位。

你需要根据原始谱、分辨率和物理背景指定 Gaussian 峰数、初始中心、宽度及必要的范围。
它们可以在整个批次共享，也可以逐谱覆盖。一个程序不应把某次实验的峰位置当成所有材料的通用常数。

## 2. 准备文件和清单

示例目录：

```text
project/
  manifest.json
  data/
    spectrum-01.csv
    spectrum-02.csv
```

CSV 第一行写列名，例如 `shift_cm1,intensity,noise_sd`。TXT 也是带表头的分隔文本；可显式指定 `delimiter`。
XLSX 使用第一行表头，`sheet` 可以是工作表名称或从 0 开始的索引。

完整清单示例：

```json
{
  "schema_version": 1,
  "data_kind": "measurement",
  "metadata": {"description": "自行填写样品和采集说明"},
  "defaults": {
    "input": {
      "x": "shift_cm1",
      "y": "intensity",
      "sigma": "noise_sd",
      "units": "cm^-1",
      "missing": "error"
    },
    "baseline": {"method": "asls", "lam": 1000000, "p": 0.01, "iterations": 30},
    "fit": {
      "window": [100, 1100],
      "peaks": [
        {"center": 350, "sigma": 20, "center_bounds": [300, 400]},
        {"center": 720, "sigma": 30, "center_bounds": [650, 800]}
      ]
    }
  },
  "outputs": {"per_spectrum_figures": true, "tracking_figures": true},
  "tracking": {
    "enabled": true,
    "prediction": "linear",
    "max_distance_cm1": 20,
    "ambiguity_margin_cm1": 1,
    "max_gap": 1,
    "width_weight": 0
  },
  "spectra": [
    {
      "id": "s01", "label": "初始压力", "series": "compression",
      "pressure_gpa": 0, "input": {"path": "data/spectrum-01.csv"}
    },
    {
      "id": "s02", "label": "第二次采集", "series": "compression",
      "pressure_gpa": 1, "input": {"path": "data/spectrum-02.csv"},
      "analysis": {"fit": {"peaks": [{"center": 353, "sigma": 20}, {"center": 718, "sigma": 30}]}}
    }
  ]
}
```

这组数值仅用于解释格式。请按实际谱段和峰形修改。没有可靠误差列时，删除 `sigma` 字段；不要编造误差。
缺少 `sigma` 时使用残差缩放的局部误差模型，而不是绝对测量误差模型。

### 同一个 Excel 文件中的多条谱

每条记录可以引用同一文件的不同强度列；横坐标列也可以不同。
例如第一条 `input` 使用 `{"path":"data/book.xlsx","sheet":"Sheet1","x":"shift","y":"P1"}`，
第二条使用相同文件和横坐标、`y` 改为 `P2`，并各自填写实际压力。
列也可以使用从 0 开始的非负索引。当前接口要求有表头，不自动判断共享 X 或交替 X/Y 布局。

## 3. 覆盖规则和数据身份

`defaults` 与每条记录的 `analysis` 按对象合并，列表整体替换。
因此，逐谱提供 `fit.peaks` 会替换整个默认峰列表，可以改变某一条谱的峰数。
`input` 的逐谱字段覆盖共享输入选项；每条记录都必须明确提供自己的 `input.path`。

切换基线方法时，例如由 `asls` 改为 `none`，新基线对象整体替换旧对象，避免遗留另一种方法的参数。
元数据可以包含中文说明，并保存在输出中。

`id` 在批次内必须唯一，忽略大小写后也不能重复；避免 Windows 保留名称，例如 CON。
显示名称用 `label`，可写中文。`series` 是明确的独立序列名称；省略时归入 `default`。
不同样品、需要分开的压缩/减压段或加热状态，请使用不同的 `series`。
程序不会从文件名或实验代号猜测这些关系。

清单顺序是采集顺序。重复压力的谱保留为不同记录，不合并、不平均。
减压序列也按采集顺序处理，压力不必单调增加。

## 4. 校验和执行

```bash
raman-analyze --lang zh validate manifest.json --check-files
raman-analyze --lang zh batch manifest.json --output runs/batch-01
```

校验检查版本、字段、必需参数、稳定 ID、追踪设置；加 `--check-files` 会读取各输入文件并检查列和数据。
拟合能否成立、峰是否可辨识，仍需在实际分析中判断。

输出目录必须是新目录。批量处理时可增加 `--no-figures` 节省绘图时间和空间；HTML 总览与数值结果仍保存。

```bash
raman-analyze --lang zh batch manifest.json --output runs/batch-02 --no-figures --json
```

`--json` 将终端摘要保持为机器可读 JSON；逐条处理进度写入 stderr。
`--lang zh` 控制帮助和自然语言提示，CSV/JSON 的字段名与状态码保持一致。

## 5. 先看离线总览，再查具体数值

用浏览器打开输出目录中的 `report.html`。它直接链接每条谱的数值文件和诊断。

- `summary.csv`：光谱身份、名称、数据类型、序列、压力、状态、峰数和结果目录。
- `spectra/<id>/processed.csv`：读取的原始强度、基线与校正后强度；输入文件不被改写。
- `spectra/<id>/fit.csv`：拟合曲线、背景、各峰分量和残差。
- `spectra/<id>/peaks.csv`：峰位、宽度、面积及局部误差。
- `tracking/assignments.csv`：所有拟合峰、配对结果、排除原因及显式缺口。
- `tracking/tracks.csv`：轨迹起止光谱和已采用点数。
- `tracking/series-*.png/.pdf`：各独立序列的压力轨迹图。

图中的圆点是采用的拟合峰位。灰色叉号表示歧义或被排除的拟合结果，缺失帧处断线。
连接线帮助阅读相邻采集点，不代表增加了测量值。未知的局部误差不会被填成零误差。

## 6. 追踪参数的含义

| 参数 | 含义 |
|---|---|
| `prediction` | `linear` 用最近两个不同压力的已采用点预测；`nearest` 使用最近采用峰位 |
| `max_distance_cm1` | 观测峰与预测中心的最大允许距离；也定义出生/缺口代价的尺度 |
| `ambiguity_margin_cm1` | 与全局最优代价相比，仍视为竞争配对的代价差范围；不是概率 |
| `max_gap` | 同一序列允许跨过的缺失采集帧数，0–50；0 表示下一帧前不可跨缺口 |
| `width_weight` | 可选的 FWHM 连续性代价权重，默认 0；需要有正的有效 FWHM |
| `enabled` | `false` 表示只做批量光谱分析，追踪明确记为 NOT_RUN |

线性预测需要局部变化近似线性。它不保证跨突变、分支变化或噪声很大的区间仍能识别同一物理峰。
对于压力重复点，程序保留每次采集，并避免用近乎相同的压力计算不稳定斜率。

程序会比较全局竞争方案，并把“属于已有轨迹”和“新峰”都纳入候选。
有多个竞争身份时记为 `AMBIGUOUS`，不赋予确定轨迹编号，也不拿该点更新预测。

## 7. 配对事件词典

| 事件 | 含义 |
|---|---|
| `NEW` | 当前采用点建立新轨迹；新编号不意味着发现了新物相 |
| `MATCHED` | 在配置的模型及阈值下唯一配对到已有轨迹 |
| `AMBIGUOUS` | 存在多个竞争身份，候选保存在 `candidate_track_ids` |
| `GAP` | 某条已有轨迹在此采集帧没有被采用的峰；`center_cm1` 留空，预测另列保存 |
| `EXCLUDED_FIT` | 拟合状态不是 PASS，因此不更新轨迹 |
| `MISSING_PRESSURE` | 缺少压力，不推造压力值，也不参与压力预测 |

编号不能自动证明化学或振动模式归属。峰完全重合、分辨率不足或前史不足时，保留歧义是正确的结果。

## 8. 调整追踪而不重拟合

```bash
raman-analyze --lang zh track runs/batch-01/batch-report.json --settings tracking.json --output runs/tracking-02
```

`tracking.json` 可以只填写需要使用的追踪参数；其余字段使用已定义的默认值。
传入该文件时它替换整套保存的追踪设置。省略 `--settings` 则复用原批次的设置。
新结果记录源批处理报告的哈希，原批次不被改写。

## 9. 部分失败和中断

单个文件找不到或列无效时，该记录记为 BLOCKED，其余记录继续。
拟合失败或不确定的数值可以保留供复核，但不会进入可信轨迹状态。

如果按 Ctrl+C，中断前已完成的记录保留在 `spectra.jsonl`，`run.json` 标为 INTERRUPTED。
进程被强制结束时，状态也可能停留在 RUNNING；这同样不代表完成。
只有 COMPLETE 的 `batch-report.json` 才可用于重新追踪。

当前不自动续跑。可以根据 journal 中已完成的 ID，把尚未完成的记录另存清单并选新输出目录继续分析；
跨两个独立批次的轨迹不会被程序默认为已经连通，需要完整序列的结果再作统一评估。

## 10. 科学与数据边界

支持的是指定的 Gaussian 模型，不是自动解释一切 Raman 光谱。局部协方差只在所选模型和固定基线条件下成立。
AsLS 参数随采样密度与噪声而变；不均匀采样应使用合适的参考区间基线或另行明确重采样。

合成示例的 truth 文件用于软件验证，不属于实验结果。真实谱的校准、物理解释和数据发布权限需由实际研究背景支持。
