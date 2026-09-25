# beamsearch —— 扩大 V1 搜索（四项开关 + N 并行）

这个目录**只新增文件，不修改 `src/` 下的任何东西**。它调用既有的 `ctocscreen.*` 数值内核
（Lambert、二体传播、高度检查、SQP 精修、独立复核），把 V1 的"宽度 1 随机贪心构造"
换成**宽度 W 的束搜索**。

## 相对历史 V1 搜索的四个开关

**规则完全照抄 `ctocscreen.constructGreedy`，唯一增加的是"束宽"（保留 W 条半成品而不是 1 条）。**
目标集合、时间提案（10 档 + 0.2 对数正态抖动）、时间预算系数、分支（全部）、高度门槛
（首段后 ≥6000 km）、末端残差、fminbnd 连续时间精修——全部与 `constructGreedy.m` 逐条一致。

| # | 开关 | 历史默认 | 现在 | 代码位置 |
|---|---|---|---|---|
| 1 | **束宽**（每层保留多少条半成品）——**唯一新增的算力** | **1** | **8** | `cfg.beam_width` |
| 2 | 每步候选目标数 | 5 | **35**（全部剩余） | `cfg.construct_target_count` |
| 3 | 时间提案 | 5 档 | **10 档**（1800…86400 s） | `cfg.construct_times_s` |
| 4 | Lambert 圈数上限 | 2 | **3** | `cfg.branch_policy.max_revolutions` |

再加两项原代码本就有、此前未启用的规则：`construct_time_budget_factor=1.6`（单段时间预算）
与 `construct_refine_count=4`（对最省的 4 条腿做 fminbnd 连续时间精修）。

开关 2/3/4 在 `expanded_20260921_a` 那轮已经用过；**开关 1 是以前从未用过的**。

## 实测：束宽到底有没有用（同等算力对照）

同一批种子、同样规则、**每个候选的算力相同**：

| 构造器 | 每次耗时 | 候选 ΔV（km/s） |
|---|---|---|
| `constructGreedy`（宽度 1）× 5 次 | 4.5–4.9 s | 24.91 / 32.05 / 约束违反 / 37.07 / 32.05 → **最好 24.91** |
| **束（宽度 3）× 3 次** | 13–15 s | 全废 ｜ 25.07 / 28.29 / 31.88 ｜ **22.80 / 23.10 / 23.97** → **最好 22.80** |

即：4.9 秒/候选 vs 14.8÷3 = 4.9 秒/候选（算力一致），**束更好而且每条构造稳定给出多条不同解**。

> 历史教训：我最初把每段候选做了抽样（只取 15/35 目标、2/10 时间档、2 个分支），
> 结果束比原构造器**差一个数量级**（44–105 km/s vs 32 km/s）。**不要减少每段的候选**，
> 束宽的收益必须建立在与原构造器相同的每段强度之上。

## 怎么运行

```matlab
cd('D:\CTOC-14\simulation')
addpath('beamsearch')

% 第一步：自检（1–3 分钟，验证能跑、并给出吞吐量）
selfTest()

% 第二步：8 任务 × 8 worker × 每任务 480 s × 束宽 8（与 campaign16 同预算）
b = run_beam_v1('beam_v1_a', 8, 8, 480, 8);

% 中断后续跑（同一 run id、同参数）
b = run_beam_v1('beam_v1_a', 8, 8, 480, 8, true);

% 放大（确认性能之后）
b = run_beam_v1('beam_v1_b', 64, 8, 1800, 24);
```

### 自己指定种子

```matlab
% 方式一：指定主种子（任务 i 的种子 = master_seed + 104729*(i-1)）
b = run_beam_v1('beam_v1_s1', 4, 8, 480, 8, false, 'master_seed', 20260922);

% 方式二：逐个任务指定（不够则循环使用）
b = run_beam_v1('beam_v1_s2', 4, 8, 480, 8, false, 'seeds', [101 202 303 404]);

% 任何 bsConfig 选项都能这样覆盖
b = run_beam_v1('beam_v1_s3', 8, 8, 900, 16, false, 'seeds', 1:8, 'refine_every', 10);
```

### 出图（只做 V1 有的图）

```matlab
report_beam_v1('beam_v1_a')                  % 基线默认 16.287951653319，参考线 7
report_beam_v1('beam_v1_a', 16.287951653319, 7)
```

与历史 V1 脚本（`tmp/archive_expanded.m`、`tmp/archive_campaign16.m`）一致，只产出：

| 产物 | 内容 |
|---|---|
| `convergence.png` | X=外层迭代，Y=总 ΔV，每个任务一条线 + 基线虚线 + 参考点线 |
| `wave_results.png` | X=任务序号，Y=各任务最好**独立复核** ΔV（V1 里 X 是波次，这里批次按任务组织） |
| `trajectory.csv` | leg, target, departure_s, arrival_s, dvx/dvy/dvz, dv_km_s, independent_error_km |

**不做** V2 的那几张（`events_and_burns.png`、`altitude.png`、`double_flyby.png`、`maneuvers.csv`、`events.csv`）。

**受控停止**：在某一任务的目录里放一个名为 `STOP` 的空文件，例如
`runs\screening\beam_v1_a_job0001\STOP`，该任务会在下一次迭代边界停止并保存检查点。

MATLAB 里运行前请确认 `startup.m` 不会拖慢启动（本机 `startup.m` 里的
`satk_initialize` 会占几分钟；如遇卡顿可用 `matlab -sd <空目录>` 绕过）。

## 文件清单

| 文件 | 作用 |
|---|---|
| `bsConfig.m` | 配置：四个开关 + 算力/并行 + 束扇出参数 |
| `bsBeamConstruct.m` | **核心**：宽度 W 的束构造，返回若干条完整候选 + 诊断 |
| `bsSearchTask.m` | 单任务搜索循环（软预算 / STOP / 检查点 / 独立复核入档） |
| `bsBatch.m` | 多任务并行（parfor），合并各任务的独立复核结果 |
| `run_beam_v1.m` | 你运行的入口 |
| `selfTest.m` | 1–3 分钟自检 |

## 产物

```
runs/screening/<run_id>_jobNNNN/checkpoint.mat   断点（可用 resume 续跑）
runs/screening/<run_id>_jobNNNN/result.mat       该任务全部结果
runs/screening/<run_id>_jobNNNN/elite.mat        该任务独立复核通过的最好解
runs/screening/<run_id>/batch.mat                合并结果（含每个任务的 archive）
runs/screening/<run_id>/elite.mat                合并后独立复核通过的最好解
```

`elite.mat` 里的 `elite` 结构包含 `candidate`（初轨 + 访问顺序 + 每段飞行时长 + 分支）
与 `evaluation`（真实脉冲、到达状态、最低高度、验证层级）。

## 算力与吞吐（实测，MATLAB R2024a 本机）

每次边求值 ≈ 一次 `enumerateBranches` + 一次 `checkArc`，实测 **≈ 1400 条边/秒**（单线程）。
每层边数 ≈ `beam_width × 本节点剩余目标数 × beam_times_per_target × beam_branches_per_edge`
（默认 8 × 剩余目标 × 3 × 2；随深度递减），一整条 35 层的束构造实测 **4–45 秒**（取决于宽度）。

**已实测跑通并记录的结果**：

| 测试 | 配置 | 结果 |
|---|---|---|
| 自检 `selfTest` | 1 任务，W=2 | 束走满 35 层，762 条边 0.9 秒；搜索循环 68 次迭代全部产出完整候选 |
| 单任务生产参数 | 1 任务 × 120 s，W=8 | **35 次迭代**；最好解 **42.19 km/s**，独立复核通过（35/35，最大误差 0.030 km） |
| 并行路径 | 2 任务 × 2 worker × 60 s，W=4 | 77 秒完成；task2 独立复核 38.04 km/s；合并导出 1 条 |

| 配置 | 搜索进程时间 | 墙钟（8 worker） |
|---|---|---|
| 8 任务 × 480 s | 64 进程分钟 | ~8–9 分钟 |
| 64 任务 × 1800 s | 32 进程小时 | ~4 小时 |
| 128 任务 × 3600 s | 128 进程小时 | ~16 小时 |

## 一个重要的调参经验（已内置）

飞行时间**不能**从时间网格里随机抽样：LEO→12h 轨道的合理飞行时间约 10800 s，随机抽到
1800 s 会产生巨大 ΔV（实测把最好解从 42 km/s 拖到 105 km/s）。现在的做法是按**转移几何**
给出时间建议（当前半径与目标半径所构椭圆的半周期）再乘 `beam_time_factors = [0.55 1 1.8]`。
这条改动把单任务 120 s 的最好解从 **105 → 42 km/s**。

## 跑完看什么

1. `batch.mat` → `batch.export_archive(1).total_dv_km_s`：**独立复核通过**的最好解
   （这才是成绩，屏幕上的"筛选值"不算）；
2. `batch.results{i}.summary`：每任务的迭代数、筛选最好值、复核最好值；
3. `result.mat` → `history`：每次迭代的模式（`beam`/`mutate`）、状态、当前最好值；
   如果 `beam` 迭代大量返回 `dead_end`/`time_limit`，说明扇出或 `beam_seconds` 需要调；
4. 与历史对比：V1 `campaign16` 16 任务最好 **16.288 km/s**（独立复核），
   目前全项目最好 **13.9765 km/s**。

## 参数调优速查

| 现象 | 调整 |
|---|---|
| 束构造常常 `time_limit`，没产出完整候选 | 降低 `beam_targets_per_node` / `beam_times_per_target`，或提高 `beam_seconds` |
| 束常常 `dead_end`（某一层没有可行边） | **提高** `beam_times_per_target` 与 `beam_branches_per_edge`（可行性主要靠更多时间/分支候选） |
| 想要更多不同序列 | 提高 `beam_width`（同时提高 `restart_probability` 会让更多迭代用于束构造） |
| 想更省算力 | 降低 `refine_every`（改成 0 可完全关掉 SQP 精修） |

## 已知风险（先说清楚）

1. **已在 MATLAB R2024a 实测**：`selfTest()`、单任务生产参数、2 任务并行三条路径都跑通了
   （上表）。但**没有**做过 8 worker 的完整批次，第一次请先跑 `selfTest()` 再跑小规模。
2. **部分候选是"刀锋轨迹"**：实测里有出现"筛选通过但独立复核不过"的候选（如某任务
   `best_verified = Inf`）。这是本项目已知现象（见
   `tmp\design_studies\consistency\robustness_criteria.md`）。默认 `export_audit_count=32`
   会复核前 32 条，能筛掉大部分；要更保险可以叠加那份文档里的 G1–G4 闸门。
3. **质量还不能与历史批次直接比**：单任务 120 s 得到 42 km/s，历史 V1 `campaign16`
   最好 16.288 km/s 是 **16 任务 × 480 s + SQP 精修**的结果。本工具默认
   `refine_every=20`（含精修），放大到 8 任务 × 480 s 才有可比性。
4. **没有任何性能承诺**：这是"把搜索空间扩大"的实验，不是"一定能得到更好解"的承诺。
