function cfg=targetedDefaults()
cfg=ctocscreen.advancedDefaults();
cfg.max_wall_s=480;cfg.window_seconds=12;cfg.window_evaluations=450;
cfg.lookahead_arcs=2;cfg.time_radius=.42;cfg.feasibility_fraction=.65;
cfg.discovery_anchors=10;cfg.discovery_samples=45;cfg.discovery_refresh=10;
cfg.relocation_weight=.015;cfg.merge_probability=.60;cfg.candidate_pool=30;
cfg.seed=20330921;cfg.checkpoint_s=30;
cfg.escape_enabled=false;cfg.stagnation_seconds=75;cfg.restart_period_seconds=150;
cfg.meaningful_improvement_km_s=.001;cfg.kick_seconds=10;cfg.kick_attempts=24;
cfg.kick_time_probability=.65;cfg.kick_time_sigma=.35;cfg.kick_duration_range=[.82 1.03];
cfg.kick_cost_ratio=2;
end
