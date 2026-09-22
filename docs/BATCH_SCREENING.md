# CTOC14 大批次筛选方案与编程设计

版本说明（2026-09-21）：本文保留 V1 批次设计。用户新要求的试射、多目标片段微调和自由机动时刻模型及仿真方案见 [MODEL_V2.md](MODEL_V2.md)；本文固定 35 段和 42 个连续变量的限制不用于 V2。两版均继续使用二体动力学，V2 设计不代表已实现。

日期：2026-09-21。状态：原始设计归档；当前已实现范围、配置与使用见 [BATCH_IMPLEMENTATION.md](BATCH_IMPLEMENTATION.md)，实测见 [BATCH_RESULTS_20260921.md](BATCH_RESULTS_20260921.md)。本文其余建议不自动代表已经实现。本阶段按用户指示不进入 ATK 仿真，忽略奖励系数 B，以总 ΔV 筛选候选方案。

正式物理模型见 [MODEL.md](MODEL.md)，强制约定见 [../AGENTS.md](../AGENTS.md)。本文规定其近似搜索阶段，不将二体结果等同于正式 J2 可行解。算法配置中的建议初值尚未经性能实验，不是赛题限制。

## 1. 本阶段产物与边界

产物是一批可复现的完整候选解：初始轨道、35 目标排列、初始滑行时间、各段飞行时长、转移分支、脉冲矢量、实际到达状态、二体复核结果和搜索记录。

- 主任务仍是 35 颗全访问、总时长不超过 864000 s、全过程高度不低于 200 km。
- 本阶段目标为 `total_dv_km_s = sum(norm(delta_v[k]))`；不计算 B，不以 F2 排序，不把耗时作为第二目标。
- 初始轨道与 42 个连续自由度保持 MODEL.md 的范围。
- 不启动、连接、运行或验证 ATK，不生成已通过官方验证的结论。
- 不建设完整 J2 精修器、不训练神经网络代理模型，也不直接优化 105 个独立脉冲分量。
- 不因当前排行榜参考值设置 ΔV 硬上限，也不将二体 ΔV 与排行榜 F2 直接作同口径胜负比较。
- 原始场景和输入数据保持不变，筛选结果写入独立运行目录。

## 2. 近似动力学：一致的二体模型

### 2.1 当前明确选择

目标卫星、初始滑行和所有转移弧统一采用地球中心引力二体模型：

```text
dr/dt = v
dv/dt = -mu * r / norm(r)^3
```

沿用 EGM96 的 mu 和 RE，以减少未来切换到 J2 时的无谓参数差异，但当前不施加 J2。参数从原始重力场文件读出或记录经核验的值与来源，不依赖库的默认地球常数。

对非圆目标及可能出现的高偏心率转移采用通用变量 Kepler 传播器；Lambert 求解器需要支持有效的零圈和多圈解族，明确区分不收敛、无解、共线退化和分支不存在。不得只支持圆轨道，不得自动禁止双曲转移；双曲/抛物转移无多圈分支，但可满足有限时段内的任务约束。

选择纯二体的目的是提高候选评价吞吐量和保证模型内部一致性。10 天内忽略 J2 可能改变访问几何、可行性与优劣排序，尤其初始低轨滑行对摄动敏感，不能给出未经验证的误差保证。保留多个不同序列、时间安排和初始轨道的精英，避免只留一个二体最优解。

未来若采用 J2 目标星历配合二体转移，应单独命名为混合近似模型，重新标记和评价，不能与本阶段同名缓存或混合排名。

### 2.2 输入加载

读取 `data/ctoc14b_targets.csv` 的以下字段：

```text
序号、名称、历元 (UTC)
X (km)、Y (km)、Z (km)
Vx (km/s)、Vy (km/s)、Vz (km/s)
```

检查唯一 ID 1..35、名称映射、统一历元、有限状态、非零位置和合理单位；核对 CSV 的惯性坐标系与原场景。不得因为只运行二体就省略坐标系确认。CSV 内根数与分类只作诊断和种子启发，传播以笛卡尔状态为准。

读取原场景和重力场文件属于离线资料核对，不需要启动 ATK。记录 CSV 文件散列、原场景来源、frame 标识、epoch、mu、RE。不得使用速度模长列替代速度分量。

## 3. 核心数据结构

用户已确定编程语言为 MATLAB。主程序、动力学、Lambert 适配、搜索、复核和测试均以 MATLAB `.m` 文件实现，不引入 Python 运行依赖。只有性能剖析证明必要且另行明确后，才考虑可选 MEX 加速；首版不以 MEX 为前提。

数据对象优先使用结构体，数值数组使用 MATLAB `double`。目标 ID 和数组下标采用 1..35，仍明确保存 ID 到行号的映射，避免输入重排改变含义。时间为相对 t0 的秒。MAT 文件保存完整数值状态；JSON 摘要中不可写 NaN/Infinity，失败数值写 null，状态另行保存。下方类型描述为协议伪代码，不是 MATLAB 语法。

```text
ProblemData
  target_ids[35], target_names[35], epoch_utc, frame
  states0[35,6], mu_km3_s2, re_km, horizon_s=864000
  input_hash, model_id="two_body_screen_v1"

Candidate
  initial_q[6]                  # MODEL.md 编码及 encoding_version
  order[35]                     # 目标 ID 的排列
  wait_s                        # >= 0
  tof_s[35]                     # > 0
  branch_ids[35]                # 分支求解后填充；不以数组顺序冒充稳定标识
  candidate_id, parent_id, seed, generation_method

TransferBranch
  branch_id                     # 圈数、几何方向、解族等规范化元数据
  v_depart[3], v_arrive[3]
  min_radius_km, end_residual_km
  status, solver_metadata

Evaluation
  status, model_id, validation_level
  planned_count, reached_planned_count, unique_visit_count
  total_dv_km_s, duration_s
  depart_times_s[35], arrive_times_s[35]
  departure_states[35,6], arrival_states[35,6], delta_v[35,3]
  min_altitude_km, max_endpoint_error_km
  constraint_diagnostics, failure_leg, failure_reason
  lambert_calls, propagation_calls, elapsed_s
```

约定 `unique_visit_count` 未做接近事件检查时为 null，不伪造为 35。`reached_planned_count` 必须来自实际连续传播的命中检查。全部 35 个不同计划目标均被命中即可证明二体模型下的全访问；全轨迹非计划接近扫描可补充首次访问时刻。

状态至少区分 `invalid_input`、`partial`、`no_branch`、`solver_failure`、`constraint_violation`、`nominal_complete`、`two_body_verified`。本阶段不得输出 `j2_verified` 或 `atk_verified`。程序异常与数值无解分别处理，不能吞掉异常当作一般差解。

## 4. 模块与接口

建议开发目录如下（仅规划，本文不创建代码文件）：

```text
src/+ctocscreen/                 # MATLAB 包命名空间
  loadProblem.m                 # 数据、单位、frame 与来源校验
  makeCandidate.m               # 候选结构构造与字段规范
  propagateTwoBody.m            # 二体传播
  initialState.m                # 初轨编码到笛卡尔状态
  targetStates.m                # 批量目标星历
  enumerateBranches.m           # Lambert 分支枚举与求解器适配
  checkArc.m                    # 弧内最低高度与末端残差
  findFlybys.m                  # 全轨迹接近事件
  chooseBranchesDP.m            # 名义固定边界上的分支动态规划
  evaluate.m                    # 连续真实状态重放、约束与 ΔV
  constructCandidates.m         # 随机化构造、束搜索、初值采样
  proposeNeighbor.m             # 重插、交换、片段搬移、破坏与修复
  repairAndRefine.m             # 时间修复、固定分支连续精修
  runSearch.m                   # 搜索预算、接受策略、多起点调度
  updateArchive.m               # 精英池更新
  saveCheckpoint.m              # 保存搜索和随机状态
  loadCheckpoint.m              # 检查与恢复
configs/screening/              # 返回配置结构体的 .m 函数
scripts/run_screening.m         # 批次入口，不承载数值核心逻辑
tests/                         # MATLAB 单元测试
runs/screening/<run_id>/
```

核心函数语义：

```matlab
problem = ctocscreen.loadProblem(csvPath, sourceMetadata);
state = ctocscreen.propagateTwoBody(state0, dtSeconds, mu);
states = ctocscreen.targetStates(problem, ids, timesSeconds);
branches = ctocscreen.enumerateBranches(r0, rf, dtSeconds, mu, branchPolicy);
check = ctocscreen.checkArc(stateDepart, dtSeconds, problem, tolerances);
path = ctocscreen.chooseBranchesDP(branchSets, initialVelocity);
evaluation = ctocscreen.evaluate(candidate, problem, config, verify);
candidates = ctocscreen.constructCandidates(problem, config, stream);
result = ctocscreen.repairAndRefine(candidate, problem, config, affectedWindow, budget);
summary = ctocscreen.runSearch(problem, config, resumeState);
```

函数不得隐式更改 ProblemData 或输入候选。随机数发生器以 `RandStream` 显式传入使用随机操作的函数。Evaluation 应包含足够诊断，不只返回一个罚函数标量。Lambert 适配器不得依赖 ATK 进程。将 `src` 加入路径，通过 `ctocscreen.*` 调用，不对包内部目录盲目使用 `genpath`。核心函数禁止使用 `clear all`、`close all` 或依赖 base workspace 的隐式变量。

实现开始时记录 MATLAB release 和已安装/许可的工具箱。基础数据、二体内核、评价器与串行搜索应能独立运行；连续精修优先计划使用 Optimization Toolbox 的 `fmincon`（SQP），并行依赖 Parallel Computing Toolbox，额外算法不能默认工具箱存在。若缺少对应依赖，明确报告：可运行串行管线、使用另行记录的自实现替代或暂时禁用该组件，但不得默默改变算法后仍使用原算法标签。测试采用 `matlab.unittest`；无需引入 Python 测试环境。

## 5. 一次候选评价

### 5.1 硬边界预检

验证初始轨道、排列、时间正值和总时限，转换初始状态，传播初始滑行。时间建议以 `wait_s/864000`、`tof_s/864000` 归一化后交给优化器，保留线性总和约束；不对所有时长自动缩放到恰好 10 天。

ec0/es0 使用半径 0.001 缩放；角变量周期处理。对严格 e0 < 0.001，以显式可配置的微小数值边界余量实现，记录它造成的受限域，不能悄悄改成允许 e0 = 0.001。

### 5.2 名义分支选择

固定候选的初轨、排列和时间，以目标精确位置为节点枚举每段有效分支。每段先排除弧内低于最低高度的分支，再构造链式动态规划：

```text
D[1,b] = norm(v_depart[1,b] - v_initial_after_wait)
D[k,b] = min_a {D[k-1,a] + norm(v_depart[k,b] - v_arrive[k-1,a])}
J_nominal = min_b D[35,b]
```

保存回溯指针，不加末端交会脉冲。连接计算复杂度为各相邻层分支数乘积之和；这不包括 Lambert 枚举成本。固定名义边界、时刻和已枚举分支集合时，该递推最优；不等于全问题最优。

不得仅因一个分支单段出发 ΔV 偏高就删除它。分支数量裁剪应保留不同到达速度的选择，并记录策略。多圈分支上限是搜索配置，需要分阶段扩大或抽查，而非物理约束。

### 5.3 连续真实状态重放

精英入正式筛选档案前，按选定分支逐段从真实位置重新求解并传播：

```text
state_before = propagate(initial_state, wait_s)
for k in 1..35:
    target_position = target_state(order[k], arrival_time[k]).position
    branch = solve_selected_branch(state_before.position, target_position, tof[k])
    delta_v[k] = branch.v_depart - state_before.velocity
    state_after = propagate([state_before.position, branch.v_depart], tof[k])
    check actual endpoint error and entire arc minimum altitude
    state_before = state_after
```

重新累计实际脉冲成本，不能沿用名义 DP 成本充当真实重放成本。分支消失或残差过大时，返回失败并允许重新分支选择/修复。全流程禁止将真实末位置或末速度重置为目标状态。

### 5.4 弧内高度与事件

二体弧的最低半径应由有限飞行区间内的候选极值决定：段端点，以及该区间内发生的近地点。椭圆多圈弧须正确处理绕行次数；双曲/近抛物情形须相应处理。由 `r·v = 0` 判定径向极值时，应可靠枚举或定位事件，不能仅用粗网格排除漏检。

不能因为完整延拓轨道的近地点低于阈值，就直接拒绝一个实际飞行区间未经过近地点的弧；这种做法只能作为明确标注的保守预筛，会损失候选。最终二体档案使用实际区间最小值。

飞越检查使用真实传播端点误差，阈值仍为 1 km。建议求解目标残差远小于该阈值；暂定 `endpoint_solver_tol_km=1e-5` 为起步配置，需通过内核验证与吞吐量实验调整。正式 1 km 阈值不得随数值配置改动。

对需要完整接近报告的精英，使用相对距离极值条件 `(r-r_target)·(v-v_target)=0` 与区间端点定位事件，合并重复目标，保存首次有效飞越。事件扫描不能替代分段状态连续性检查。

## 6. 搜索流程

### 6.1 生成不同起点

对初始轨道面、相位、初始滑行和时间分配作分层或低差异采样。均匀采样 cos(i) 可用于初始轨道面覆盖，但只是一种种子分布，不限制倾角搜索范围。

构造方式混合：随机排列、带随机性的低成本扩展、有限宽度束搜索。每个部分状态包含已访问集合、当前真实状态、绝对时刻、累计 ΔV 和剩余时间。用 64 位位掩码表示 35 个目标集合，不按累计 ΔV 一项合并终端速度不同的状态。

扩展时尝试多个下一目标、转移时长与分支；优先考虑低成本候选，同时保留随机探索。对剩余目标保留时间的启发式可以使用，但未经证明的剩余成本估计不能当作严格下界永久剪枝。对相同完整访问深度与验证层级的方案比较 ΔV；不能把较短部分序列的低成本当作更优完整解。

### 6.2 全访问序列改进

从精英或探索状态选择父解，使用：

1. 单点重插、两点交换。
2. 2—4 个连续目标片段搬移或短片段重排。
3. 随机移除、围绕高成本衔接移除、按连续时间窗口移除一批目标，再重新插入。

起步的破坏规模建议为 3—8 个目标，可自适应扩大，但不是固定模型限制。重插必须完成所有目标，不能把删除目标后的低 ΔV 记为成功改善。

序列修改后先做局部时间修复，再比较方案。局部窗口固定外部访问时刻，重分配内部相邻段时长；至少包含与前后片段连接的脉冲代价，尤其需要计入窗口末速度变化对下一段出发脉冲的影响。

使用累计 dt 时，一段时长变化可能平移后续所有访问时刻，缓存失效范围必须随之扩大。不能按静态 TSP 的两条边差值直接更新成本。可以在局部优化中使用等价的绝对访问时刻坐标，结束后转换回 dt 存储，以减少不必要的远端变化。

### 6.3 连续变量精修

固定排列与分支，先对局部时间窗口精修；有潜力的候选再对全部 42 个连续变量精修。MATLAB 实现优先评估 `fmincon` 的 SQP 方法，以线性不等式表示时间预算、非线性约束表示偏心率圆盘等条件。遇到非光滑或分支附近病态问题可在局部采用明确记录的无导数方法；使用 `patternsearch` 等工具箱方法前核验依赖。差分进化作为初值或困难片段的备选，可用 MATLAB 实现，不为每个序列重新启动大规模全维种群。

梯度计算期间固定分支标识；分支不可用则返回明确诊断，由可行性恢复机制处理。不能在差分扰动两侧任意切换最小成本分支而仍假设目标光滑。角度差采用周期差值，步长在归一化尺度上配置，数值传播误差应小于导数差分信号。

交替运行：固定分支连续精修 -> 重新枚举并 DP 选分支 -> 连续真实状态重放。设置每个候选的评价预算，避免一个难解候选耗尽全部资源。

### 6.4 接受与精英维护

始终保存最好二体全访问且通过复核的方案。探索链可在温度调度下接受较差但仍可行的完整方案，例如按 `exp(-(J_new-J_old)/temperature)` 接受，temperature 单位为 km/s，并从观测到的成本变化尺度初始化。

部分解与不可行解放在独立构造/修复池，不能通过低 ΔV 替换可行精英。需要时使用约束违反量指导恢复，但不把它混入最终 ΔV 或宣称改变正式目标。

精英池分为最低成本集合和多样性集合。相同序列的不同时间/初轨解可能进入不同局部最优区域，不能仅按排列去重；相同目标函数值也不代表重复轨迹。多样性可由有向相邻目标对、归一化访问时刻差、初始轨道差共同衡量。

### 6.5 主循环伪代码

```text
load and validate problem/config
construct or resume initial candidate pools
while global budget remains:
    parent = choose elite or exploration candidate
    trial = change sequence / perturb times / perturb initial orbit
    trial = repair permutation and time feasibility
    quick = nominal two-body evaluate(trial)
    if promising or selected for exploration audit:
        trial = local time refine with fixed branches
        trial = branch DP and sequential replay
        if promising:
            trial = full continuous refine within budget
            trial = branch DP and sequential replay again
        if two-body complete and verified:
            update best, diversity archive and exploration state
    update operator statistics
    checkpoint periodically
replay final archive using declared verification tolerances
export screening results, never official feasibility claims
```

## 7. 批处理、缓存与并行

- 先测量单次 Kepler 传播、单次 Lambert 枚举和完整评价的耗时、失败率，再确定总预算；不承诺未经测量的每秒候选数。
- 批量计算目标状态，先使用直接二体传播；仅在测出瓶颈后增加插值，并用随机时刻误差检查保证其筛选精度。
- 无并行工具箱时先运行串行 `for`。具备 Parallel Computing Toolbox 时，独立起点使用 `parfor`，或以 `parfeval` 调度有预算的独立任务；按实际环境建立进程池。worker 持有只读 ProblemData 和有界本地缓存，必要时通过 `parallel.pool.Constant` 复用。避免嵌套并行和线性代数线程过量。此处为未来 MATLAB 程序并行，不要求调用多代理开发。
- 使用由 master seed 和稳定任务 ID 确定的 `RandStream` 子流，而非根据 worker 调度顺序分配随机序列。记录流类型、任务 ID、seed/substream 与状态。异步精英交换可能改变搜索轨迹，不宣称改变 worker 数仍能逐位复现；严格复现模式按确定性批次和任务 ID 合并。
- 目标星历缓存键含 model_id、输入散列、目标 ID 和时间；Lambert 缓存键还含两端位置/来源、两端时刻、mu、求解容差及分支策略。
- 最终精英复核使用精确输入缓存或禁用近似缓存。网格量化时间/位置的缓存只允许作为明确标记的粗筛，不用于细粒度求导或最终二体复核。
- 缓存只按需填充并设置 LRU/内存上限，不预建所有目标对与两个时间轴的高维张量。
- 对局部修改可复用未变化的名义弧；实际状态重放从首个受影响点继续，除非确认边界状态完全一致，不能盲目复用后缀。

## 8. 起步运行配置

配置文件须完整保存并带 schema_version。建议先提供 smoke、pilot、batch 三档；下列数字仅用于起步实验，不是优化过的参数。

| 参数 | smoke | pilot | batch |
|---|---:|---:|---:|
| 独立种子数 | 2 | 8 | 32 |
| 墙钟时间预算/次运行 | 60 s | 1800 s | 21600 s |
| 束宽（启用束搜索时） | 8 | 32 | 64 |
| 总精英档案目标容量 | 8 | 32 | 64 |
| 局部精修最大函数评价/候选 | 30 | 150 | 300 |
| 全变量精修最大函数评价/候选 | 0 | 600 | 1500 |

所有档位仍使用 35 目标正式输入；smoke 不保证找到完整解，只检查管线。小规模数值测试另外标记目标子集。总预算、每候选预算和停止条件同时生效，长批次支持取消与恢复。

另外必须配置：

```text
model_id = two_body_screen_v1
objective = total_dv_km_s
master_seed, workers, deterministic_batches
max_lambert_calls, max_wall_time_s
branch_policy = staged                     # 先零圈，再扩大分支集合
max_revolutions = 2                        # 起步筛选上限，非题目限制
tof_numerical_min_s = 1.0                   # 可调算法下界，需做放宽敏感性实验
eccentricity_margin, lambert_tolerances, kepler_tolerances
endpoint_solver_tol_km = 1e-5               # 起步目标精度
flyby_limit_km = 1.0                        # 正式阈值，不随调参改变
min_altitude_km = 200.0                     # 正式阈值
cache_memory_limit_mb, checkpoint_interval_s
operator_weights, destroy_size_range, annealing_policy
```

多圈上限与正时间下界的敏感性实验是批次报告的一部分。扩充分支集合后必须重新计算候选成本，不能复用不带分支策略版本的旧排名。初始滑行允许严格为零，不强行加下界。初轨角度与轨道面不能只在某类卫星附近采样而排除全局探索。

## 9. 运行产物与断点恢复

```text
runs/screening/<run_id>/
  config.json           # 完整配置、软件版本、代码版本/散列
  problem_metadata.json # 来源、单位、frame、数据散列、model_id
  progress.jsonl        # 时刻、调用数、最优 ΔV、失败率、算子统计
  summary.json          # 总耗时、有效全访问数、档案最优与验证层级
  elites.csv            # 候选 ID、ΔV(km/s 与 m/s)、时长、验证级别等
  candidates/<id>.json  # 初轨、排列、时长、分支、seed 与父解
  trajectories/<id>.mat # MATLAB 各段状态、脉冲与复核诊断
  checkpoint.mat        # 搜索状态、随机状态、操作权重和档案
```

检查点以 MATLAB `save` 写入同目录临时 MAT 文件，完成后在文件系统支持时以原子重命名替换；不得先删除唯一有效检查点。恢复时核对模型版本、数据散列与配置兼容性。普通记录使用常规 MAT 格式，大型数组再考虑 `-v7.3`；使用 `jsonencode` 导出元数据、`writetable` 导出摘要。输出目录唯一，不覆盖此前运行。未存完整密集轨迹时，可由候选参数重算；复现依据不能只有总 ΔV。

报告至少包含：最好二体全访问 ΔV 随墙钟时间/调用数变化、不同种子的最好值分布、名义到连续重放的成本差、各类失败数、每个算子收益与花费、缓存命中率、多圈/时间下界敏感性。终止时没有完整二体解，应如实报告，不导出部分解冒充最佳可行解。

## 10. 开发顺序与验收

1. 数据与数值内核：加载校验、轨道转换、Kepler 传播、Lambert 分支。用已知圆轨道周期、传播往返、能量/角动量守恒及 Lambert 端点重构验证；覆盖椭圆、多圈、近共线与双曲情形。
2. 评价器与 DP：小分支集合用穷举核对 DP；构造“局部最省分支并非全链最省”的测试。确认脉冲使用真实前段速度，确认无末端匹配费用。
3. 约束检查：覆盖端点安全但中段低于 200 km、延拓近地点低但实际弧未经过、时间超限、重复目标、非有限输出与求解失败。用独立数值积分交叉核对精英的二体传播结果，无需 ATK。
4. 搜索最小闭环：多起点构造、单点重插、局部时间修复、精英存档、断点恢复；先跑 smoke 与 pilot，测吞吐量再扩展。
5. 强化搜索：加入交换/片段搬移、破坏—重建、全变量精修、分支扩大、多进程批次。相同墙钟预算和多个种子比较组件收益，不只看单次最优值。
6. 输出候选档案：连续重放所有导出精英，明确标记 two_body_verified。后续 J2 或 ATK 检查需要单独阶段、单独结果，不在本阶段执行。

本次只归档方案与编程细节；未创建求解器、安装依赖、执行优化批次或声称获得任何新成绩。
