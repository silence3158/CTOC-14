# MATLAB 大批次筛选实现与使用

日期：2026-09-21。实现版本：two_body_screen_v2。物理模型遵守 MODEL.md；本文件记录搜索近似和实际代码。所有单位为 km、s、rad、km/s。

## 已实现流程

1. `loadProblem` 按 ID 排序读取 35 个目标，检查名称和精确历元，计算真实 SHA-256；离线读取原始 ATK XML 核对笛卡尔状态及 EGM96 引力参数。没有启动 ATK。原始资料不改动。
2. `constructGreedy` 从可变初始轨道与初始滑行开始，在随机目标子集、扰动时间网格、短程/长程及配置圈数范围内生成可行访问路径。每次使用巡检器实际到达速度，随机受限候选选择用于增加多样性。这是构造启发式，不保证全局或序列最优。
3. `mutateCandidate` 对精英执行交换、移位、局部时长扰动、初始轨道六变量和初始滑行扰动；`prepareBranches` 对固定顺序和时刻进行名义分支 DP；`evaluate` 使用实际到达位置重新求解指定分支并执行真实脉冲，不重置位置或速度。
4. `refineCandidate` 固定排列与分支，用 Optimization Toolbox 的 `fmincon` SQP 优化全部 42 个连续变量。变量经过尺度归一化；时限为线性约束，偏心率圆盘为非线性约束。求解器内的不可行罚值仅用于搜索，不能进入可行档案；保存实际访问过的最好可行点和 exitflag/output，不把求解器退出等同于可行或最优。
5. `runSearch` 执行多起点与精英改进、历史记录、MAT 检查点、预算与 STOP 文件检查。`runBatch` 将独立种子任务分配给进程池，最后重新重放并合并候选。
6. `verifyIndependent` 使用 `ode113` 从历元连续传播飞行器，逐段施加已记录的原脉冲矢量，不再重新瞄准；目标也独立积分。逐个检查 35 个唯一目标的真实段末接近事件。ODE 相对容差为 3e-14、绝对容差为 1e-14；通用变量传播默认相对残差容差为 5e-15，并处理浮点括区间收敛。每段高度仍用有限弧解析径向极值检查。只有独立复核通过的路径进入 `export_archive` 和 `elite.mat`。

## 数值内核与工具箱选择

本机 MATLAB R2024a 已核验 Optimization Toolbox、Global Optimization Toolbox、Parallel Computing Toolbox 的许可。实际使用 Optimization Toolbox (`fmincon`) 和 Parallel Computing Toolbox (`parpool`/`parfor`)；没有借用不存在的轨道 Lambert 接口。安装路径中的 `lambert` 是 Mapping Toolbox 的地图投影函数。

Lambert 内核仍为项目 MATLAB 实现，但已去掉固定 z 网格扫描和根四舍五入：零圈自适应括根，多圈按相邻奇点区间用 `fminbnd` 求最短时间，再用 `fzero` 分别求两侧根。每个根都重构并传播检查端点，分支 ID 编码几何方向、圈数与左右根，缺失指定分支直接失败。共线退化明确返回状态，不假装求解成功。

二体传播使用带括区间的通用变量 Kepler 求解与稳定 Stumpff 函数，椭圆传播先化简周期；失败路径的诊断结构保持一致。高度检查改为椭圆、双曲线和近抛物线有限时间区间内的解析极值，不再用稀疏采样宣布安全。

参考工具文档：[fmincon](https://www.mathworks.com/help/optim/ug/fmincon.html)、[fmincon 的并行有限差分](https://www.mathworks.com/help/optim/ug/using-parallel-computing-with-fmincon-fgoalattain-and-fminimax.html)。本实现选择在独立搜索任务之间并行，内部 SQP 不嵌套并行。

## 运行方法

在 MATLAB 中：

```matlab
cd('D:\CTOC-14\simulation')
addpath('src','scripts')

% 完整小批次：两个任务串行，各最多 5 次外层迭代。
batch = run_batch_screening('smoke', 0, 'smoke_01', false);

% 较大批次：8 个独立任务，4 个进程，每任务最多 1000 次或 3600 秒。
batch = run_batch_screening('batch', 4, 'night_01', false);

% 相同配置继续已有任务；已有检查点的任务从其随机状态恢复。
batch = run_batch_screening('batch', 4, 'night_01', true);
```

精确控制预算：

```matlab
p = ctocscreen.loadProblem();
c = ctocscreen.defaultConfig('batch');
c.run_id = 'custom_01';
c.tasks = 16;
c.workers = 4;
c.max_candidates = 2000;
c.max_wall_s = 1800;
c.refine_every = 25;
c.refine_evaluations = 150;
b = ctocscreen.runBatch(p, c);
```

预算以每个任务计，不是整批总预算；包含各次候选构造与精修，但不包含 MATLAB/进程池启动及结束时独立导出复核。时间上限在迭代或求解器回调边界检查，因此是软预算。每任务独立种子为 master_seed + 104729*(任务序号-1)，模 2^32。相同迭代预算下结果可复现；时间预算可能因机器负载而停在不同迭代。

不提供并行许可时会警告并使用相同独立任务的串行路径；缺少 Optimization Toolbox 且启用 SQP 时明确报错，不替换算法。可显式设置 `refine_every=0` 只运行构造、DP 与局部变异。已有进程池大小与 workers 不符时明确报错，由调用者决定是否重建。

## 可配置搜索近似

- `branch_policy.max_revolutions=2`：每段最多枚举 2 圈，是搜索裁剪，不是赛题限制。
- `branch_policy.endpoint_tol_km=1e-4`：Lambert 数值端点容差；正式访问阈值仍为 1 km。
- `min_tof_s=60`：优化搜索下界，不是物理模型约束；直接评价器接受所有正时长。
- 构造时间网格 `[900 1800 3600 7200 14400]` 秒再加对数扰动，不是优化时长上界。后续时长可自由变化，只受总时限和搜索下界限制。
- `construct_target_count=5`：先探索最多 5 个剩余目标；无解时扩展到全部剩余目标。
- `construction_min_altitude_after_first_km=6000`：仅初始构造时，第二段起偏好较高的弧内最低高度，以减少过于敏感的掠地路径。后续变异和精修仍按正式 200 km 下限评价。可改为 200，属于构造启发式选择，不改变模型或导出标准。
- 初始构造偏心率采样小于 0.0009、初始滑行先采样 0–3600 s；变异与全变量 SQP 可扩展到模型允许域，不固定 600 km 圆轨道。
- `export_audit_count=16`：结束时独立积分检查前 16 个搜索精英。未检查者不能视为独立复核通过。设为 0 可做纯吞吐实验，此时没有新导出精英。

## 输出与恢复

每任务目录：`runs/screening/<run_id>_jobNNNN/`。

- `checkpoint.mat`：迭代、随机流状态、档案、历史、SQP 报告、数据与配置签名；保存前复制上一代到 `checkpoint.mat.bak`，不先删除旧文件。此为可恢复双代保存，不声称操作系统原子事务。
- `screened_elite.mat`：按逐段瞄准重放得到的最好搜索候选。它不代表独立固定脉冲重放通过。
- `elite.mat`：仅当独立 ODE 连续重放通过时存在，包含报告；失败不得冒充导出解。
- `result.mat`：全部搜索结果，`archive` 是搜索库，`export_archive` 是独立二体复核库，`export_audits` 保留通过或失败报告。
- 主目录 `batch.mat`：合并结果；选择 `batch.export_archive` 获取可用导出候选。

恢复时核验 CSV 哈希、MATLAB 包源代码哈希、模型版本、全部输入状态、引力常数与配置。仅迭代/时间预算和显示/保存周期可调整，其他不同明确拒绝。旧 v1 检查点不兼容。损坏时可显式读取 `.bak` 恢复，不静默跳过损坏。任务目录创建名为 `STOP` 的文件可在下一轮边界停止并保存；继续前删除该标记。

## 验证边界

连续固定脉冲轨迹可能对微小传播误差极敏感；逐段 Lambert 端点很准并不足以保证其独立积分仍能访问全部目标。开发时已实际发现该问题，所以新增独立导出门槛。`evaluate.status=two_body_verified` 的兼容字段仅表示逐段二体筛选通过；查看 `evaluation.validation_level` 和 `independent_report` 区分导出复核。

本模块从未执行 J2 或竞赛版 ATK 验证，也没有计算 B 或网站 F2。坐标保持原 ATK 笛卡尔轴系，CSV 与场景已逐项比对；其更具体的 ICRF/J2000 标识尚未独立认证，在 J2/ATK 阶段必须确认。二体中心引力筛选不另行变换这些轴。

仍未实现：J2 微分修正、ATK 场景生成、网站提交、全局最优证明、近共线专用 Lambert 处理、自动扩大圈数边界。当前结果不是竞赛成绩。

## 回归测试

```matlab
results = runtests('tests/TestScreening.m');
assertSuccess(results);
```

覆盖独立 ODE 多分支 Lambert 检查、双曲/抛物弧内近地点、长时间双曲传播、完整访问与速度连续传递、初始约束拒绝、指定分支缺失、部分路径拒绝入库、DP 真正失败段、恢复与不中断运行一致性、圆轨道全链独立积分。实际批次成绩和运行证据见 `BATCH_RESULTS_20260921.md`。