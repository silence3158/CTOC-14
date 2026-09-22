> **历史归档说明**：本文记录的是早期 `CTOC13.xml` / `ZPY15.xml` 离线数据提取工作，仅用于资料核对，不属于当前 CTOC14 乙题首版 MATLAB 筛选管线。当前正式约定见 `../../AGENTS.md`，模型见 `../MODEL.md`，筛选方案见 `../BATCH_SCREENING.md`；当前目标入口为 `../../data/ctoc14b_targets.csv`。

# ZPY15.xml 卫星数据提取方法记录（历史归档）

## 1. 任务与结论

- 源文件：`D:\CTOC-14\往届题目\周培源力学竞赛15团体赛\ZPY15.xml`
- 源文件编码：XML 声明为 `GB2312`，实际可用 `gb18030` 解码。
- 源文件大小：`9243702` 字节。
- SHA256：

```text
E51B2760AA3B9047E41F0253089216BAD488E76E58F1EAA660E678674344395D
```

**结论：ZPY15.xml 确实包含卫星数据。文件中有 345 颗卫星，名称为 `Debris1` ~ `Debris345`。每颗卫星的 `Orbit` 节点给出同一历元的 ECI 直角坐标位置和速度。本文将每颗卫星的直角坐标状态转换为经典轨道六根数，并生成与 `starlink卫星轨道数据.xlsx` Sheet1 相同列结构的 CSV。**

## 2. 源文件结构

- 根节点：`Scenario Name="ZPY15"`
- 共 345 个 `Satellite` 对象，名称 `Debris1` ~ `Debris345`
- 每颗卫星的关键数据位于 `Satellite/Orbit`：

| 字段 | 说明 |
|---|---|
| `OrbEpoch` | 轨道历元，UTC |
| `StartUTC` / `StopUTC` | 仿真起止时间 |
| `PositionX/Y/Z` | ECI 位置，m |
| `VelocityX/Y/Z` | ECI 速度，m/s |
| `TleElem` | 每个卫星都有该节点，但全为占位值（`SSCNum=99999`、其余为 0），不包含有效 TLE |

所有卫星的示例历元均为：

```text
2030-11-14 08:00:00
```

## 3. 处理方法

1. 以二进制读取 `ZPY15.xml`；
2. 按 `utf-8-sig`、`gb18030`、`gbk` 等编码依次尝试解码；
3. 用 Python `xml.etree.ElementTree` 解析；
4. 遍历根节点下所有 `Satellite`；
5. 从 `Satellite/Orbit` 读取 `OrbEpoch`、`PositionX/Y/Z`、`VelocityX/Y/Z`；
6. 使用地球引力常数：

```text
mu = 3.986004418e14 m^3/s^2
```

将 ECI 直角坐标状态转换为经典轨道六根数：

- `a`：半长轴，m
- `e`：偏心率
- `i`：倾角，deg
- `Omega`：升交点赤经，deg
- `omega`：近地点幅角，deg
- `MM`：平近点角，deg

7. 将 `OrbEpoch` 转换为 `MJD2000`：

```text
MJD2000 = (OrbEpoch - 2000-01-01 12:00:00) / 86400
```

8. 按 `starlink卫星轨道数据.xlsx` Sheet1 的列结构写出 CSV。

## 4. 输出文件

### 4.1 `ZPY15_satellites_starlink_format.csv`

位置：

```text
D:\CTOC-14\simulation\data\old\ZPY15_satellites_starlink_format.csv
```

列结构与 starlink xlsx Sheet1 一致：

| 列 | 表头 | 说明 |
|---|---|---|
| A | 序号 | 1 ~ 345 |
| B | `name` | `Debris1` ~ `Debris345` |
| C | `catalogNumber` | 源文件 `TleElem/SSCNum`，当前全为 `99999` 占位值 |
| D | `epoch[MJD2000]` | 由 `OrbEpoch` 转换 |
| E | `a` | 半长轴，m |
| F | `e` | 偏心率 |
| G | `i（deg）` | 倾角，deg |
| H | `Omega（deg）` | 升交点赤经，deg |
| I | `omega（deg）` | 近地点幅角，deg |
| J | `MM（deg）` | 平近点角，deg |

共 345 行数据。

### 4.2 `ZPY15_satellites_source_states.csv`

位置：

```text
D:\CTOC-14\simulation\data\old\ZPY15_satellites_source_states.csv
```

保留源文件的原始直角坐标状态，便于核查：

| 列 | 说明 |
|---|---|
| `index` | 序号 |
| `name` | 卫星名称 |
| `catalogNumber` | `SSCNum` |
| `orb_epoch_utc` | 历元 |
| `start_utc` / `stop_utc` | 起止时间 |
| `position_x_m` / `position_y_m` / `position_z_m` | ECI 位置，m |
| `velocity_x_mps` / `velocity_y_mps` / `velocity_z_mps` | ECI 速度，m/s |
| `source_file` | `ZPY15.xml` |
| `source_xpath` | 对应的 XML 路径 |

## 5. 结果校验

对 345 颗卫星的转换结果进行统计：

| 轨道根数 | 最小值 | 最大值 | 平均值 |
|---|---:|---:|---:|
| 半长轴 a / m | 7,081,316.79 | 11,023,552.34 | 7,308,469.69 |
| 偏心率 e | 0.00328 | 0.36631 | 0.03427 |
| 倾角 i / deg | 98.0018 | 100.0033 | 98.0891 |
| 升交点赤经 Omega / deg | 73.8483 | 96.2727 | 75.9257 |

全部卫星均为有界椭圆轨道，未出现 e >= 1 或 a <= 0 的异常值。

## 6. 注意事项

- `catalogNumber` 来自源文件 `SSCNum`，当前全部为 `99999`。如果需要真实编目号，需要提供带真实 SSC 编号的数据源。
- `TleElem` 全为占位值，因此本文没有使用 TLE，而是使用 `Orbit` 中的直角坐标状态。
- 转换得到的是历元时刻的瞬时（osculating）开普勒根数；若后续用 J2 或高精度模型传播，仍应以原始直角坐标状态为初值。
- `MM（deg）` 按 `starlink卫星轨道数据.xlsx` 和 `Main.m` 中的用法表示平近点角，而不是平均运动。
- 全部卫星使用同一个历元 `2030-11-14 08:00:00`，对应 `epoch[MJD2000]` 约为 `11274.833333333334`。

## 7. 复现命令

```powershell
python D:\CTOC-14\simulation\tools\old\extract_zpy15_satellites.py
```

脚本会自动在 `D:\CTOC-14` 下递归查找 `ZPY15.xml`，读取 345 颗卫星并写入上面的两个 CSV。



