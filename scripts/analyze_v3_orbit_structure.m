function report=analyze_v3_orbit_structure(inputFile,outputDir)
%ANALYZE_V3_ORBIT_STRUCTURE Diagnose plane, shape, and phase-like changes.
% Uses a saved fixed-pulse replay only; it never searches or propagates.
if nargin<1 || isempty(inputFile)
    inputFile=fullfile(fileparts(fileparts(mfilename('fullpath'))), ...
        'runs/v3/diagnostics/tail27_20260925_method_search/comparison.mat');
end
if nargin<2 || isempty(outputDir)
    outputDir=fullfile(fileparts(fileparts(mfilename('fullpath'))), ...
        'runs/v3/diagnostics/tail27_20260925_orbit_structure');
end
assert(isfile(inputFile),'Input diagnostic MAT is missing: %s',inputFile);
assert(~isfile(fullfile(outputDir,'orbit_structure.mat')), ...
    'Refusing to overwrite existing diagnostic evidence.');
if ~isfolder(outputDir), mkdir(outputDir); end
sim=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(sim,'src'));
loaded=load(inputFile,'result'); result=loaded.result;
eph=ctocscreen.v3LoadTargetEphemeris(fullfile(sim, ...
    'runs/v3/preprocessing/targets_20260923_release/target_ephemeris.mat'));
mu=eph.model.mu; re=eph.model.re;
s=result.reconstructed_schedule; tr=result.reconstructed_trace;
n=numel(s.maneuver_times_s); assert(n==35);
% The saved trace stores the initial state followed by each post-burn state.
assert(numel(tr.times)>=n+1 && size(tr.states,1)>=n+1);
x0=ctocscreen.initialState(s.initial_q,mu).';
post=tr.states(2:n+1,:); pre=post;
pre(:,4:6)=post(:,4:6)-s.delta_v_km_s;

preEl=zeros(n,9); postEl=zeros(n,9); rt=zeros(n,3); dvnorm=zeros(n,1);
plane=zeros(n,1); incChange=zeros(n,1); fullPlane=zeros(n,1); planeFromInitial=zeros(n,1);
for k=1:n
    preEl(k,:)=elements(pre(k,:),mu);
    postEl(k,:)=elements(post(k,:),mu);
    rvec=pre(k,1:3); vvec=pre(k,4:6); h=cross(rvec,vvec); R=rvec/norm(rvec); N=h/norm(h); T=cross(N,R);
    rt(k,:)=[dot(s.delta_v_km_s(k,:),R),dot(s.delta_v_km_s(k,:),T),dot(s.delta_v_km_s(k,:),N)];
    dvnorm(k)=norm(s.delta_v_km_s(k,:));
    plane(k)=acosd(max(-1,min(1,dot(N,postEl(k,7:9)))));
    incChange(k)=angledeg(preEl(k,7:9),postEl(k,7:9));
    fullPlane(k)=angledeg(preEl(k,7:9),postEl(k,7:9));
end
% 9 element order: a, e, peri, apo, period, energy, hx, hy, hz-normal.
% Replace normal components with a unit vector for the plane columns.
for k=1:n
    preEl(k,7:9)=cross(pre(k,1:3),pre(k,4:6)); preEl(k,7:9)=preEl(k,7:9)/norm(preEl(k,7:9));
    postEl(k,7:9)=cross(post(k,1:3),post(k,4:6)); postEl(k,7:9)=postEl(k,7:9)/norm(postEl(k,7:9));
    if k==1, n0=preEl(k,7:9); end
    planeFromInitial(k)=angledeg(n0,postEl(k,7:9));
end
% Cumulative phase-like proxy: Kepler period and osculating true anomaly.
phasePre=zeros(n,1); phasePost=zeros(n,1);
for k=1:n
    phasePre(k)=trueAnomaly(pre(k,:),mu);
    phasePost(k)=trueAnomaly(post(k,:),mu);
end
targetId=s.witness_times_s; %#ok<NASGU>
T=table((1:n).',s.delta_v_km_s(:,1),s.delta_v_km_s(:,2),s.delta_v_km_s(:,3),dvnorm, ...
    preEl(:,1),postEl(:,1),preEl(:,2),postEl(:,2),preEl(:,3),postEl(:,3), ...
    preEl(:,4),postEl(:,4),preEl(:,5),postEl(:,5),phasePre,phasePost, ...
    incChange,fullPlane,planeFromInitial,rt(:,1),rt(:,2),rt(:,3), ...
    'VariableNames',{'burn','dv_x','dv_y','dv_z','dv_norm','a_pre','a_post','e_pre','e_post', ...
    'peri_pre','peri_post','apo_pre','apo_post','period_pre','period_post', ...
    'true_anomaly_pre_deg','true_anomaly_post_deg','inclination_change_deg', ...
    'plane_change_deg','plane_from_initial_deg','dv_radial','dv_tangential','dv_normal'});
writetable(T,fullfile(outputDir,'orbit_structure.csv'));

% Summaries after fixed visit cutoffs, including the full task.
cutoffs=[2 8 17 27 35]; S=zeros(numel(cutoffs),12);
for j=1:numel(cutoffs)
    q=1:cutoffs(j); S(j,:)=[cutoffs(j),sum(dvnorm(q)),sum(abs(rt(q,1))),sum(abs(rt(q,2))), ...
        sum(abs(rt(q,3))),sum(abs(rt(q,3)))/sum(dvnorm(q)), ...
        min(postEl(q,1)),max(postEl(q,1)),min(postEl(q,2)),max(postEl(q,2)), ...
        min(postEl(q,5)),max(postEl(q,5))];
end
summary=table(S(:,1),S(:,2),S(:,3),S(:,4),S(:,5),S(:,6),S(:,7),S(:,8),S(:,9),S(:,10),S(:,11),S(:,12), ...
 'VariableNames',{'through_burn','dv_total','abs_dv_radial','abs_dv_tangential','abs_dv_normal', ...
 'normal_fraction','a_min','a_max','e_min','e_max','period_min','period_max'});
writetable(summary,fullfile(outputDir,'orbit_structure_summary.csv'));
fig=figure('Visible','off','Color','w','Position',[100 100 1200 800]);
tiledlayout(2,2,'Padding','compact');
nexttile; plot(1:n,postEl(:,1),'o-'); xlabel('Burn'); ylabel('a (km)'); grid on;
nexttile; plot(1:n,postEl(:,3),'o-',1:n,postEl(:,4),'o-'); xlabel('Burn'); ylabel('Osculating apsis (km)'); legend('peri','apo','Location','best'); grid on;
nexttile; plot(1:n,postEl(:,2),'o-'); xlabel('Burn'); ylabel('e'); grid on;
nexttile; plot(1:n,cumsum(dvnorm),'k-',1:n,cumsum(abs(rt(:,3))),'b--'); xlabel('Burn'); ylabel('km/s'); legend('total DV','abs normal component','Location','best'); grid on;
exportgraphics(fig,fullfile(outputDir,'orbit_structure.png')); close(fig);
    report=struct('input_file',inputFile,'output_dir',outputDir,'n_burns',n, ...
    'total_dv_km_s',sum(dvnorm),'sum_abs_normal_km_s',sum(abs(rt(:,3))), ...
    'normal_fraction',sum(abs(rt(:,3)))/sum(dvnorm),'max_plane_change_deg',max(fullPlane), ...
    'max_plane_from_initial_deg',max(planeFromInitial), ...
    'max_inclination_change_deg',max(incChange),'summary',summary,'table',T);
save(fullfile(outputDir,'orbit_structure.mat'),'report','T','summary','pre','post','-v7');
fprintf('ORBIT STRUCTURE total %.9f km/s; abs normal %.9f (%.3f%%); max per-burn plane %.6f deg; max from initial %.6f deg\n', ...
    report.total_dv_km_s,report.sum_abs_normal_km_s,100*report.normal_fraction,report.max_plane_change_deg,report.max_plane_from_initial_deg);
disp(summary);
end

function y=elements(x,mu)
r=x(1:3); v=x(4:6); rr=norm(r); vv=dot(v,v); h=cross(r,v); hn=norm(h);
ev=cross(v,h)/mu-r/rr; e=norm(ev); a=1/(2/rr-vv/mu); peri=a*(1-e); apo=a*(1+e); period=2*pi*sqrt(abs(a)^3/mu); energy=vv/2-mu/rr;
y=[a,e,peri,apo,period,energy,h/hn];
end
function d=angledeg(a,b), d=acosd(max(-1,min(1,dot(a,b)/(norm(a)*norm(b))))); end
function nu=trueAnomaly(x,mu)
r=x(1:3); v=x(4:6); rr=norm(r); h=cross(r,v); evec=cross(v,h)/mu-r/rr; e=norm(evec);
if e<1e-10, nu=atan2d(r(3),r(1)); return; end
nu=atan2d(dot(cross(evec,r),h)/norm(h),dot(evec,r)); if nu<0,nu=nu+360;end
end
