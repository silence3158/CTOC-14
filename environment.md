# 运行环境要求

文件目录必须按照下面的编排，因为代码会自检引力常数文件

'''text
<根目录>\
├── problem\ctoc14乙题\
│   ├── CTOC14B.atk                              728.9 KB
│   └── ATK-CTOC14\AstroData\Earth\EGM96.grv     137.2 KB
└── simulation\
    ├── src\+ctocscreen\*.m                        必需
    ├── beamsearch\*.m                             必需
    └── data\ctoc14b_targets.csv                   11.4 KB
'''

检查目录是否正确：

'''
cd('<解压目录>\simulation')
addpath('src'); addpath('beamsearch')
p = ctocscreen.loadProblem();
p.target_ids
p.source_audit
selfTest()
'''



需要matlab优化工具包Optimization Toolbox，并行运算工具包Parallel Computing Toolbox

运行方式：

在matlab终端中运行下面的代码，即可开始计算

'''
cd('D:\你的根目录\simulation')
addpath('beamsearch')
b = run_beam_v1('beam_v1_s1b', 32, 8, 480, 8, false, 'master_seed', 11451400);
'''

参数解释：

'''
cd('D:\你的根目录\simulation')
addpath('beamsearch')
b = run_beam_v1('beam_v1_s1b',  运行标签，结果写到 runs\screening\beam_v1_s1b\
                    32,  独立任务数，32 个互不相干的搜索，各自一个种子、一份检查点
                    8,  并行任务数，现在开8个并行
                    480, 每任务搜索软预算（秒），到点在迭代边界停；一次束构造/SQP 可能让它略超
                    8, 束宽 W，每层保留 8 条半成品解
                    false,  是否续跑，false=全新跑；true 需 run id 与参数完全一致
                    'master_seed',
                    11451400);  设置主种子数值，任务 i 的种子 = 11451400 + 104729×(i-1)
'''

这条命令里标签后面那四个数字分别是：

32 是这次要跑多少个互相独立的搜索任务，每个任务有自己的随机种子，跑得越多越可能撞到更好的解；

8 是同时并行跑几个任务；480 是每个任务最多搜多少秒，搜得越久每个任务做得越深、结果越好；

最后的 8 是束宽，也就是搜索过程中每一步同时保留几条候选路线。

总时间可以按"任务数 × 每任务秒数 ÷ 并行数"估，所以现在这条命令大约是 32 × 480 ÷ 8，也就是 32 分钟左右，实际运行会多几分钟用于启动和结束时的独立复核。

想加大算力，可以把任务数从 32 改成 64，时间翻倍到约 64 分钟。其次是把每任务 480 秒改成 900 秒，时间同样翻倍。再就是把束宽从 8 改成 16。

如果你机器的物理核心足够多，还可以把并行数从 8 改成 16，时间直接减半，但先运行 feature('numcores') 确认核数，并行数超过物理核心反而会更慢。

最后三条注意：每次改完都要把第一个参数（运行标签，比如 beam_v1_s1b）换成新名字，否则会报目录已存在；最好一次只改一个数字，不然分不清是哪个起了作用；跑完不要看屏幕上那个筛选值，成绩要看 load('runs\screening\你的标签\batch.mat') 之后的 batch.export_archive(1).total_dv_km_s，那才是独立复核出来的数。

运行完成后，使用下面的代码绘图并生成报告

'''
report_beam_v1('beam_v1_s1b')
'''

report的内容是你定下来的运行标签