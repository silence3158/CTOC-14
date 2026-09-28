# GitHub 推送修复记录

日期：2026-09-28。用户要求解决本地仓库无法推送到 GitHub 的问题，SSH 密钥口令由用户在本机窗口输入；随后明确要求那个大文件不上传云端。

## 原因与处理

1. Windows ssh-agent 未启动，现有 SSH 密钥未加载。GitHub 能识别公钥，但最初无法完成签名认证。启动 agent 后，用户自行解锁密钥，`ssh -T` 确认账号为 `silence3158`。未读取、记录或复制私钥/口令。
2. 仓库原 main 的历史包含 `runs/fragments/fragment120_20260921/batch.mat`，大小 125.64 MiB。虽然提交 `55726ec` 已删除当前版本中的该文件，初始提交仍包含它，普通推送依然会触发 GitHub 的 100 MiB 限制。

仓库本地 `core.sshCommand` 设置为 `C:/Windows/System32/OpenSSH/ssh.exe`，统一使用 Windows agent。agent 的启动类型保留为 Manual；如果重启后没有运行，可在 PowerShell 中运行 `Start-Service ssh-agent`，然后用 `ssh-add "$env:USERPROFILE\.ssh\id_ed25519"` 由用户解锁密钥。如启动服务提示权限不足，再使用管理员终端。

## 原始历史如何保留

- 修复前 main：`53fc9e137dbbf5730918a0fd44c57327ec18c783`，103 条提交。
- 完整备份：`D:/CTOC-14/git-backups/simulation-before-github-20260928-154828.bundle`。`git bundle verify` 确认完整历史可用。
- 原 main 改名为 `codex/archive-before-github-20260928`；提交内容不变，旧回滚标签不变。
- 旧标签清单：`D:/CTOC-14/git-backups/original-tags-20260928.txt`。
- 独立上传副本：`D:/CTOC-14/github-publish-20260928`。只在此副本中过滤指定文件，再导入处理后的 main。没有在原始归档上执行过滤或清理。
- 过滤后、添加本记录前的提交：`fb2eb026b242895a60f89599b2e3abe871d50c9e`，仍为 103 条提交；最新文件树与原 main 完全相同，树哈希均为 `610a396bfe842efcdfe2ac467f9c56124fd5157c`。
- 新旧提交对应表：`D:/CTOC-14/git-backups/github-commit-map-20260928.txt`。去掉历史文件会改变发布副本的提交编号，旧编号仍可通过原分支、标签和 bundle 恢复。

工具为官方 [git-filter-repo](https://github.com/newren/git-filter-repo) v2.47.0 独立脚本，未全局安装；SHA-256：`67447413e273fc76809289111748870b6f6072f08b17efe94863a92d810b7d94`。操作仅为 `--path runs/fragments/fragment120_20260921/batch.mat --invert-paths`。其他历史文件未排除，不使用 Git LFS。

## 以后怎么推送

工作目录仍为 `D:/CTOC-14/simulation`，日常使用 main，通过 origin 推送至 `git@github.com:silence3158/CTOC-14.git`。SSH 解锁后正常执行 `git push` 即可（最终连接状态在收尾检查确认）。

旧归档分支和旧标签会重新引入被排除的大文件，因此保留在本地，不执行 `git push --all`、`git push --mirror` 或批量推送旧标签，也不要把旧归档分支合并回 main。未来需要迁移某个旧改动时，应只迁移具体改动。现有 `.gitignore` 已忽略 runs 目录。

原远程别名 `orgin` 与 `origin` 指向同一仓库；本轮使用 origin，不擅自删除旧别名。此次文件排除与历史接入为用户明确授权的发布处理，不扩展为实验记录清理。

GitHub 对普通 Git 文件大小的规定见[官方说明](https://docs.github.com/en/repositories/working-with-files/managing-large-files/about-large-files-on-github)。本轮修复不涉及算法改动或实验。
