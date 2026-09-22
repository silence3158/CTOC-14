# CTOC14 运行归档：内核自检与 smoke 闭环

日期：2026-09-21  
环境：MATLAB R2024a（24.1.0.2537033），Windows  
用户配置目录：`simulation/.matlab_pref`

## 1. 本次运行目标

验证早期 MATLAB 二体筛选闭环是否能正常启动、加载目标数据、执行数值内核并跑通一次 smoke 流程。

## 2. MATLAB 启动问题

最初 MATLAB 启动失败：

```text
Fatal Startup Error: failed to load settings errors_warnings plugin
```

在用户提供 `simulation/.matlab_pref` 后，MATLAB R2024a 可以正常启动并执行 `-batch` 命令。

本次运行使用的调用形式：

```powershell
matlab -batch "cd('D:\CTOC-14\simulation'); addpath('src'); run('tests\run_kernel_checks.m')"
matlab -batch "cd('D:\CTOC-14\simulation'); run('scripts/run_screening.m')"
```

## 3. 内核自检结果

命令：

```matlab
cd D:\CTOC-14\simulation
addpath('src')
tests\run_kernel_checks
```

结果：**通过**，exit code = 0。

关键输出：

```text
Loading problem data...
  target count OK
  circular period round-trip OK (err 1.375e-11)
  energy conservation OK
  Lambert zero-rev reconstruction OK (err 4.944e-11)
  DP smoke OK (cost 3.236)
All kernel self-checks passed.
```

覆盖内容：

- `loadProblem` 读取 35 个目标；
- `initialState` 与二体传播圆轨道周期往返；
- 能量守恒检查；
- Lambert 零圈分支重构；
- `chooseBranchesDP` 小规模 DP smoke。

原始日志：

```text
simulation/docs/logs/20260921_kernel_test.log
```

## 4. smoke 闭环运行结果

命令：

```matlab
cd D:\CTOC-14\simulation
run scripts\run_screening.m
```

结果：**闭环正常执行**，exit code = 0；但随机初始候选均未生成可行档案。

关键输出：

```text
Loaded 35 targets from D:\CTOC-14\simulation\data\ctoc14b_targets.csv
CTOC14 two-body screening run 20260921_154412_smoke
candidates: 3, seed: 20260921
iter 1/3: status=no_branch leg=1 reason=No preselected safe branch reference on leg 1. dv=Inf
iter 2/3: status=no_branch leg=1 reason=No preselected safe branch reference on leg 1. dv=Inf
iter 3/3: status=no_branch leg=1 reason=No preselected safe branch reference on leg 1. dv=Inf

Screening summary:
            model_id: "two_body_screen_v1"
    validation_level: 'two_body_screen'
           n_archive: 0
        best_dv_km_s: Inf
         best_status: 'none'
```

运行产物：

```text
simulation/runs/screening/20260921_154412_smoke/checkpoint.mat
```

原始日志：

```text
simulation/docs/logs/20260921_smoke.log
```

## 5. 本轮发现并修复的问题

1. **`constructCandidates` 结构体数组预分配错误**  
   空 struct 数组不能接收带字段的候选结构体，已改为首个候选直接赋值，后续元素再扩展。

2. **`propagateTwoBody` 与 `targetStates` 的 info 结构体字段不一致**  
   成功返回包含 `f/g/fdot/gdot`，失败返回缺少这些字段，导致结构体数组赋值失败。已统一所有返回路径字段。

3. **`enumerateBranches` 标量结构体初始化错误**  
   `emptyBranchStruct()` 返回 0×1 数组，不能直接赋字段。已增加 `branchPrototype()` 标量结构体。

4. **测试脚本变量遮蔽内置 `inf`**  
   `[st, inf] = ...` 覆盖了内置 `inf`，导致 `bestErr = inf` 和比较失败。已改用变量名 `info` 和数值 `1e300`。

5. **名义 DP 未排除低高度分支**  
   原 DP 仅按 ΔV 选择分支，可能选中弧内低于 200 km 的长程分支。现已在 DP 前用 `checkArc` 逐段过滤安全分支。

6. **名义分支与真实重放的分支对应不稳定**  
   `chooseBranchesDP` 现返回选中的 branch 结构体，`evaluate` 按 `geometry`、`rev_estimate` 和最近 `z` 匹配实际重放分支，不再只依赖重新编号的 `branch_id`。

7. **缺少安全分支时错误回退为贪心选择**  
   若预选安全分支不存在，`evaluate` 现在显式返回 `no_branch`，不会退回选择低高度分支。

## 6. 当前结论

- 二体筛选闭环已经可以运行到结束，并能显式报告失败。
- 当前 smoke 随机候选失败的原因是：在 200 km 最低高度硬约束下，随机初始轨道与时间安排导致第一段没有满足条件的 Lambert 安全分支。
- 这属于候选构造与搜索质量问题，不是本次的启动或语法问题。
- 当前结果不能作为可行解或优化成绩。

## 7. 下一步

1. 改进 `constructCandidates`：按 Lambert 可达性筛选目标和时间，不要纯随机生成 35 段。
2. 在 `runSearch` 中加入简单的可行构造策略，至少先生成 `partial` 候选。
3. 实现非计划接近事件扫描，填充 `unique_visit_count`。
4. 对非椭圆弧实现近地点事件检查，替换采样回退。
5. 在 smoke 跑通可行候选后，再扩展到 pilot/batch 搜索算子与性能测试。

