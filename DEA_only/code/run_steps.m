function run_steps()
% Run Figure 2 with the same DEA-only physics and equilibrium as Figure 3.
root = fileparts(fileparts(mfilename('fullpath')));
out = fullfile(root,'results','step_damping');
if ~isfolder(out), mkdir(out); end
original_path = path; restore_path = onCleanup(@()path(original_path));
addpath(fullfile(root,'code'),'-begin');
[base,study] = params(); f = physics(); rows = cell(10,1); curves=cell(2,5); kk=0;
supports={'spring','gas'}; dampings=[100 300 1000 3000 10000];
for i=1:2
 name=supports{i}; p=rmfield(base,{'spring','gas'}); p.support=name;
 specific=base.(name); keys=fieldnames(specific);
 for j=1:numel(keys),p.(keys{j})=specific.(keys{j});end
 if strcmp(name,'gas'),p.c1_Pa_s_m=p.c1_prime_Pa_s_m;p.c2_Pa=p.c2_prime_Pa;end
 for j=1:numel(dampings)
  p.c1_Pa_s_m=dampings(j);p=f.initial(p);
  T=solver(p,[0 study.on_s .05],[0 base.pulse_V],p.dt_s);
  fine=p;fine.rel_tol=study.check_rel_tol;fine.abs_tol=study.check_abs_tol;
  D=solver(fine,[0 study.on_s .05],[0 base.pulse_V],p.dt_s/2);
  assert(all(isfinite(T{:,:}),'all'));
  err=max(abs(T.radius_m-D.radius_m));assert(err<study.radius_check_limit_m);
  area=(T.radius_m+T.thickness_m).^2-T.radius_m.^2;
  area_error=max(abs(area/area(1)-1));assert(area_error<1e-12);
  drift=max(abs(T.radius_m(T.time_s<study.on_s)-p.r0_m));assert(drift<study.equilibrium_limit_m);
  T.displacement_um=1e6*(T.radius_m-p.r0_m);curves{i,j}=T;kk=kk+1;
  [peak,ip]=max(T.displacement_um);
  rows{kk}=struct('support',name,'damping_N_s_m3',dampings(j),...
   'peak_um',peak,'peak_time_ms',1000*T.time_s(ip),'final_um',T.displacement_um(end),...
   'tight_radius_error_m',err,'area_rel_error',area_error,'pre_onset_drift_m',drift);
  writetable(T,fullfile(out,sprintf('step_%s_c1_%g.csv',name,dampings(j))));
 end
end
summary=struct2table(vertcat(rows{:}));writetable(summary,fullfile(out,'step_summary.csv'));
save(fullfile(out,'step.mat'),'base','study','curves','summary','dampings');
plot_steps(out,fullfile(root,'figure'));disp(summary);
end
