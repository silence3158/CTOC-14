> **历史归档说明**：本文记录的是早期 `CTOC13.xml` / `ZPY15.xml` 离线数据提取工作，仅用于资料核对，不属于当前 CTOC14 乙题首版 MATLAB 筛选管线。当前正式约定见 `../../AGENTS.md`，模型见 `../MODEL.md`，筛选方案见 `../BATCH_SCREENING.md`；当前目标入口为 `../../data/ctoc14b_targets.csv`。

# CTOC13.xml 数据提取方法记录（历史归档）

## 1. 任务与结论

- 源文件：`D:\CTOC-14\往届题目\CTOC13-丙题\CTOC13.xml`
- 源文件编码：XML 声明为 `GB2312`，实际可用 `gb18030` 解码。
- 提取脚本：`D:\CTOC-14\simulation\tools\old\extract_ctoc13_xml.py`
- 输出目录：`D:\CTOC-14\simulation\data\old`

**重要结论：`CTOC13.xml` 中包含的不是卫星轨道数据，而是“20 个地面固定目标 + 1 个海上移动目标状态”。场景里没有 `Satellite` 对象，也没有开普勒根数 `a/e/i/RAAN/omega/M`。因此无法按 `starlink卫星轨道数据.xlsx` 的卫星轨道格式填充轨道根数列。**

辅助验证：
- 场景根节点为 `Scenario Name="CTOC13"`；
- 根下包含 20 个 `Facility`，名称为 `Target1` ~ `Target20`，每个只有 `Position/Lat`、`Position/Lon`、`Position/Alt`；
- 根下包含 1 个 `Ship`，其中有 `Route`（2 个航路点）和 `HistoryData`（289 个状态点）；
- `Scenario1.xml`、`Scenario2.xml` 都是只有空 `Group` 的场景，也没有轨道对象；
- 丙题说明文档说明该题是地面目标覆盖/重访问题，不是星间多目标交会。

## 2. 源文件校验

```text
SHA256: D22D24D2C92AD7D1B379AAD26732AD2AC91B2444BF9B5FF7A7730993A5B634DB
大小:   194215 bytes
```

## 3. 提取流程

1. 以二进制读取 `CTOC13.xml`；
2. 按 `utf-8-sig`、`gb18030`、`gbk`、`big5`、`utf-16`、`utf-16-le`、`utf-16-be` 依次尝试解码；
3. 用 Python `xml.etree.ElementTree` 解析 XML；
4. 读取 `Scenario` 下的元数据字段；
5. 遍历所有 `Facility`，提取 `Name` 和 `Position` 下的 `Lat`、`Lon`、`Alt`；
6. 读取 `Ship/Route/WayStartUTC` 和所有 `WayPt`，提取航路点；
7. 读取 `Ship/HistoryData`：
   - 第 1 行为基准 UTC 时间；
   - 后续每行为一个状态点：`Time_s` + 10 个数值；
   - 已知含义为：`Time_s, X_m, Y_m, Z_m, Vx_mps, Vy_mps, Vz_mps`；
   - 末尾 4 个数值暂记为 `Extra1` ~ `Extra4`，推测为姿态四元数，原样保留；
8. 所有 CSV 使用 `utf-8-sig` 编码写出，便于 Excel 直接打开。

## 4. 输出文件

### 4.1 `CTOC13_scenario_meta.csv`

| 列名 | 说明 |
|---|---|
| Key | 字段名 |
| Value | 字段值 |

包含 `StartTime`、`StopTime`、`StepSize`、`Mode` 等。

### 4.2 `CTOC13_fixed_targets.csv`

20 个地面固定目标。

| 列名 | 说明 |
|---|---|
| Index | 序号 1~20 |
| Name | `Target1` ~ `Target20` |
| ObjectType | 固定为 `Facility` |
| Lat_deg | 纬度，度 |
| Lon_deg | 经度，度 |
| Alt_m | 高度，m |
| SourceFile | `CTOC13.xml` |
| SourceXPath | 对应的 XML 路径 |

### 4.3 `CTOC13_ship_route.csv`

海上移动目标的航路点。

| 列名 | 说明 |
|---|---|
| Index | 航路点序号 |
| Name | `Ship` |
| WayStartUTC | 航路起始 UTC |
| Time_s | 相对航路起始时刻的秒数 |
| Lat_deg | 纬度，度 |
| Lon_deg | 经度，度 |
| Alt_m | 高度，m |
| SourceFile | `CTOC13.xml` |
| SourceXPath | 对应的 XML 路径 |

### 4.4 `CTOC13_ship_history_ecef.csv`

海上移动目标的 289 个历史状态点，时间步长 600 s。

| 列名 | 说明 |
|---|---|
| Index | 状态点序号 |
| Time_s | 相对基准时刻的秒数 |
| UTC_ISO | 换算后的 UTC 时间 |
| X_m, Y_m, Z_m | 位置，m（从数值看为 ECEF/地固系） |
| Vx_mps, Vy_mps, Vz_mps | 速度，m/s |
| Extra1~Extra4 | 源文件中未使用的 4 个附加数值，推测为姿态四元数 |
| SourceFile | `CTOC13.xml` |

## 5. 当前 `starlink卫星轨道数据.xlsx` 的真实格式

对 `D:\CTOC-14\simulation\data\old\starlink卫星轨道数据.xlsx` 解析后，工作表 `Sheet1` 的表头为：

| 列 | 表头 | 含义 |
|---|---|---|
| A | 序号 | 目标编号 |
| B | `name` | 卫星名称 |
| C | `catalogNumber` | 编目号 |
| D | `epoch[MJD2000]` | 历元 |
| E | `a` | 半长轴，m |
| F | `e` | 偏心率 |
| G | `i（deg）` | 倾角，deg |
| H | `Omega（deg）` | 升交点赤经，deg |
| I | `omega（deg）` | 近地点幅角，deg |
| J | `MM（deg）` | 平近点角，deg |

`Main.m` 中 `load_starlink_targets` 按 B~J 列读取，并根据 a 的数量级自动把 m 转成 km。  
当前 xlsx 另有 `Sheet2` 用于变量变化范围；但 `parameters1.m` 期望的服务航天器工作表名为 `Service`，当前 xlsx 中没有该工作表，这一点在直接运行参考代码时需要另行处理。

## 6. 为什么 CTOC13 不能生成同格式轨道 CSV

`CTOC13.xml` 中不存在卫星对象，也没有上述 E~J 列所需的半长轴、偏心率、倾角、升交点赤经、近地点幅角、平近点角。  
场景里的 20 个目标是地球表面固定点，海上目标也是地固系/经纬高航迹。强行把经纬度填进轨道根数列会产生物理上错误的输入。因此本次提取没有伪造轨道六根数，而是按实际对象类型输出地面目标与移动目标数据。

## 7. 复现命令

```powershell
python D:\CTOC-14\simulation\tools\old\extract_ctoc13_xml.py
```

脚本会自动从 `D:\CTOC-14` 下递归查找 `CTOC13.xml`，并写入 `D:\CTOC-14\simulation\data\old`。

## 8. 后续可选工作

- 如果后续找到包含 `Satellite` 和 `Keplerian/Classical` 的 ATK 想定 XML，可新增一个转换器，输出与 `starlink卫星轨道数据.xlsx` 完全一致的 10 列 CSV；
- 如果目标是做 CTOC13 式地面目标覆盖/重访问题，应使用本方法提取的 `CTOC13_fixed_targets.csv` 和 `CTOC13_ship_history_ecef.csv`，建立卫星/星座对地面和海上目标的可见性模型；
- 如果仍要把 CTOC13 数据接入 `Main.m` 的星间交会序列优化模型，需要先明确物理建模方式，不能直接把地面点当作轨道目标。



