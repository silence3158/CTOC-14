# V3 必备基础预处理实施记录

日期：2026-09-23。对应用户授权的五项基础工作：输入/物理配置核验、35 个目标的摄动星历、连续插值及误差检查、缓存元数据与失效检查、MATLAB 批量查询接口。仅完成基础预处理，不代表集束蚁群、统一联合优化或 SA 已实现。

## 1. 已实现范围与尚未核实的条件

| 项目 | 实现 | 边界 |
|---|---|---|
| 输入与动力学审计 | CSV 与原 XML 35 个唯一目标、初始状态、历元核对；核对中心天体与力模型设置；读取 EGM96 常数和 C20 | XML 原生惯性轴系与竞赛引擎的未来 EOP 策略仍待核实 |
| 摄动目标星历 | 对全部 35 个目标在 `[0,864000] s` 数值积分中心引力 + J2 | 当前为明确标注的 nominal J2，不是 ATK 已对齐星历 |
| 连续插值 | 分段五次 Hermite，端点使用位置、速度、加速度；速度为同一位置多项式的解析一阶导数 | 每段检查不等于严格的全区间误差上界 |
| 缓存 | MAT、模型和配置元数据、输入/代码 SHA-256、MATLAB 版本、文件完整性摘要 | 变更任一签名依赖需重新生成，不静默复用旧缓存 |
| 查询 | 成对查询、标量广播、目标×时间笛卡尔积；可返回物理加速度 | 禁止区间外外推；查询热路径不访问文件 |

用户补充：CSV 直接从 `.atk` 提取，不能回忆是否指定过轴系，认为应为地心系。**地心原点不等于确定惯性坐标轴**。代码不旋转原始 CSV 数值，暂以 GCRF 解释，明确写入 `frame=ATK_native_assumed_GCRF`、`official_alignment_verified=false`。后续需以竞赛版 ATK 的明确说明或逐时刻对照消除此项；当前未启动 ATK，也未调用其 DLL/MEX。未核实时禁止将此缓存当作正式对齐结果。

## 2. 数据、单位与地球定向

所有原文件保持不变。入口 `data/ctoc14b_targets.csv`；原场景 `D:\CTOC-14\problem\ctoc14乙题\CTOC14B.atk`；地球常数来自同目录 `ATK-CTOC14/AstroData/Earth/EGM96.grv`。

- 内部采用 km、s、km/s、km/s²；历元 UTC `2035-01-01 12:00:00`。
- μ = 398600.4415 km³/s²，RE = 6378.1363 km。
- 全归一化 C20 = −4.841653717360e−4；J2 = −√5 C20。采用瞬时笛卡尔状态传播，不使用平均根数传播器。
- 极轴来源：随题 `AstroData/IERS-conventions/2010/IAU2006_XYS.dat`，表头标明 TT/TDT、每日采样、X/Y/s 单位角秒，覆盖 1974-12-15 至 2050-01-15；本地按表头推荐的 9 阶 Lagrange 插值读取 X、Y，并构造 `[X,Y,sqrt(1-X²-Y²)]`。
- TT = UTC + 37 + 32.184 s；37 s 来自附带闰秒表最后记录。2035 年真实/竞赛闰秒策略尚未确认。
- 极轴预存为 300 s 网格；运行时线性插值并归一化；全部中点和直接 XYS 插值比较并记录差异。这是第二层插值的检查，不独立校准原始 XYS 表。
- 附带 EOP 数据止于 2025-12-27，覆盖不到任务历元。当前取零极移、零 dCIP；不能称为 ATK 的既定回退行为。轴对称 J2 在零极移假设下不依赖绕极轴的地球自转角，故无需预存完整地固旋转矩阵或设置 DUT1 来驱动它。
- `Earth.cb`、EOP、闰秒表和 XYS 表均纳入输入签名；`Earth.cb` 仅作随题参考，不声称竞赛程序实际采用其中每项配置。

本机 MATLAB R2024a 能使用基础积分器。最初尝试 Aerospace Toolbox 的 `dcmeci2ecef`，实际调用因产品许可不可用失败；已移除该运行依赖。当前实现仅使用 MATLAB 基础能力（包括 JVM/XML/SHA-256），没有 Python、Aerospace、Parallel Computing Toolbox 的运行依赖。

## 3. 数值检查与发布规则

生产积分器 `ode113`；交叉参考 `ode89` 使用单独书写的径向/轴向 J2 公式。两者均使用同一明确记录的极轴模型，因此交叉一致不能证明该模型已与 ATK 对齐。

最终默认两者 `RelTol=3e-14`，位置 `AbsTol=1e-12 km`、速度 `AbsTol=1e-15 km/s`，`MaxStep=60 s`。插值初始步长 120 s，逐目标按需减半至最低 15 s；这些都是计算配置，不是题目限制。

在**每一段**的三个内点 `u = 0.211324865405187, 0.5, 0.788675134594813` 及全部节点检查查询结果与参考积分的差异。默认最大检查位置差≤0.1 m，速度差≤0.0001 m/s。分别记录插值与生产积分的差异、两个积分器的差异，避免将积分误差误认为只需细化插值即可消除。无随机抽查依赖。

任何目标未通过则报错，不发布有效缓存。全部通过且生成前后的源签名相同时，先写临时 MAT 再移为正式文件，并生成 `.sha256`；缺少摘要或摘要不一致的文件禁止加载。下游仍须对最终任务使用固定原脉冲独立积分，不能将星历数值检查当作 35/35 飞越或高度验证。

## 4. 入口与查询协议

从 `D:\CTOC-14\simulation` 启动 MATLAB：

```matlab
addpath('scripts');
folder = run_v3_preprocessing('my_targets_01');  % 使用新的目录名
results = test_v3_preprocessing(folder);
```

生成完整 10 天目标缓存是本次已授权的预处理，不运行优化实验。配置可作为第二参数传入，例如 `struct('initial_step_s',60)`；未知字段或非正数值报错。新目录要求避免覆写现有结果。

加载一次，再多次查询：

```matlab
addpath('src');
file = fullfile(folder, 'target_ephemeris.mat');
eph = ctocscreen.v3LoadTargetEphemeris(file);

% 一个目标，一个或多个时刻：N×3
[r,v] = ctocscreen.v3QueryTargets(eph, 12, [0 1000 864000]);

% 全部35个目标，同一时刻：35×3；a 为中心引力+J2物理加速度
[r,v,a] = ctocscreen.v3QueryTargets(eph, [], 12345);

% 成对：目标2在1000s、目标7在2000s
[r,v] = ctocscreen.v3QueryTargets(eph, [2 7], [1000 2000]);

% 每个目标在每个时刻：2×3×3
[r,v] = ctocscreen.v3QueryTargets(eph, [2 7], [0 1000 2000], 'grid');
```

目标号、时间使用 `double`。默认 `pairs` 模式：一个输入为标量可广播，否则长度必须一致。支持重复、无序目标和时刻；`[]` 目标表示全部 35 个。时间闭区间包含首尾，不允许负数、NaN/Inf 或超过 864000 s。`a` 来自物理力模型，不是五次多项式二阶导数。

每个工作进程加载一次即可；不要在每次目标函数调用中重复读取 MAT、计算签名。未来并行搜索可在 worker 初始化时加载，本次未引入并行依赖。当前查询无需访问 ATK 或即时重新积分目标。

`v3LoadTargetEphemeris(file,true)` 要求正式轴系/动力学对齐，当前缓存必须报 `alignmentPending`。默认允许显式标注的 nominal 模型供本地开发；严禁通过手改字段伪装对齐。

## 5. 文件与复现

源码位于 `src/+ctocscreen/v3*.m`；命令入口为 `scripts/run_v3_preprocessing.m`、`scripts/test_v3_preprocessing.m`；回归测试 `tests/TestV3Preprocessing.m`。

输出 `runs/v3/preprocessing/<label>/`：

- `target_ephemeris.mat`、`target_ephemeris.mat.sha256`：系数、初态、模型、签名、配置、验证与查询计时。
- `validation.csv`：35 行逐目标误差、实际步长、样本量、时间。
- `report.md`：自动生成的动力学假设、审计误差、精度和性能结果。
- `test_results.mat`：显式执行测试入口后生成的测试结果。

签名覆盖 CSV、原 XML、EGM96、Earth.cb、EOP、闰秒、XYS 以及预处理/查询和源核验函数；记录 MATLAB 完整版本与绝对路径。变更源、版本或移动整个工作区后重新生成。配置、模型参数、schema 保存在 MAT 内，MAT 文件摘要同时保护这些字段；摘要用于完整性检查，不是防恶意篡改的签名。

开发过程记录：`targets_20260923` 因 Aerospace 许可失败未发布；`targets_20260923_v2` 在 Target33 达不到 0.1 m 门槛后拒绝发布；`targets_20260923_final` 为收紧积分容差前启动的重复尝试，已主动中止；`targets_20260923_validated` 全部目标数值检查通过，但性能报告使用数字开头字段名导致发布前报错，随后修正。上述尝试均未发布有效缓存，目录保留，不能算成功批次。最终运行与测试结果见下节。

## 6. 最终执行结果

最终运行标签 `targets_20260923_release`，已成功发布并加载，完整流程进程退出码为 0。自动日志位于 `runs/v3/preprocessing/release_console.log`；测试汇总 `runs/v3/preprocessing/test_summary.log`。

| 实测项 | 结果 |
|---|---:|
| 目标 / 时间覆盖 | 35 / 全部 864000 s |
| 检查状态点总数 | 1,036,835 |
| 最大检查位置差 | 0.0238114026 m |
| 最大检查速度差 | 2.05228705×10⁻⁵ m/s |
| 插值步长 | 34 颗 120 s；Target32 为 60 s |
| 生成与数值验证时间（不含后续性能计时/写盘/测试） | 104.873 s |
| MAT 文件大小 | 32,163,217 bytes，约 30.67 MiB |
| 极轴 300 s 插值中点最大差 | 1.02302921×10⁻¹³ rad |
| CSV/XML 最大位置分量差 | 4.94517735×10⁻⁷ km |
| CSV/XML 最大速度分量差 | 4.97433206×10⁻¹⁰ km/s |
| 回归测试 | 12 passed / 0 failed / 0 incomplete |
| MATLAB `checkcode` | 8 个包函数、2 个入口、1 个测试类均 0 findings |

预热后每组 100 次、重复 3 组取中位数；只衡量已加载缓存的查询，不包含磁盘读取/源码哈希，也不作为优化器提速倍数：

| 查询 | 每次时间 |
|---|---:|
| 一颗目标一个时刻 | 0.01409 ms |
| 一颗目标 1000 个时刻 | 0.098591 ms |
| 35 颗目标同一时刻 | 0.245064 ms |
| 35 对目标/时刻 | 0.252856 ms |
| 35 颗目标同一时刻含加速度 | 0.248513 ms |

测试涵盖初态一致、首尾时刻、成对/网格对应及重复无序目标、位置/速度导数一致、节点连续、越界/非法输入、强制正式对齐拒绝、源签名失效、缓存损坏拒绝、赤道/极轴 J2 解析特例和全部目标误差预算。代码检查最初的三处多余抑制标记已经移除并复查；检查日志分别为 `codecheck.log` 与 `test_summary.log`。

直接使用本次缓存，无需重新生成：

```matlab
cd('D:\CTOC-14\simulation');
addpath('src');
folder = fullfile('runs','v3','preprocessing','targets_20260923_release');
eph = ctocscreen.v3LoadTargetEphemeris(fullfile(folder,'target_ephemeris.mat'));
[r,v,a] = ctocscreen.v3QueryTargets(eph, [], 3600);
```

当前结论为 **nominal J2 数值查询基础可用**。输入数值核验已完成，正式轴系/EOP 对齐仍待核实；没有运行搜索、生成新的巡检解或验证任何竞赛成绩。
