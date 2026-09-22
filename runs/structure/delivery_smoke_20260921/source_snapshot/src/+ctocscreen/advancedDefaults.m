function cfg=advancedDefaults()
cfg=ctocscreen.fragmentDefaults();
cfg.multi_seconds=6; cfg.multi_evaluations=1500; cfg.multi_polish_evaluations=120;
cfg.multi_velocity_window=2; cfg.templates={[0 3],[2 2],[1 1 2],[0 2 2]};
cfg.beam_width=6; cfg.beam_roots=4; cfg.beam_times=[.65 1 1.35];
cfg.beam_keep_per_parent=6; cfg.beam_template_trials=1; cfg.beam_seconds=120;
cfg.beam_warm_prefixes=false;
cfg.joint_seconds=120; cfg.joint_evaluations=5000; cfg.joint_iterations=150;
cfg.multi_total_seconds=120; cfg.seed=20320921;
end
