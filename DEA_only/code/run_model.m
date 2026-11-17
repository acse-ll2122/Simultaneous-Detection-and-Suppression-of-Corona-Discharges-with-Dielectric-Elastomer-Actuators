function run_model()
% Run the current Figure 3 pulse study and verify every saved trajectory.
root = fileparts(fileparts(mfilename('fullpath')));
out = fullfile(root,'results','pulse_response');
if ~isfolder(out)
    mkdir(out);
end
[base,study] = params();
supports = {'spring','gas'};
support_rows = cell(size(supports));
f = physics();
for i = 1:numel(supports)
    p = support_params(base,supports{i});
    edges_s = [0,study.on_s,study.pulse_off_s,study.end_s.(p.support)];
    assert(all(diff(edges_s)>0),'Pulse times must increase within each simulation.');
    levels_V = [0,base.pulse_V,0];
    radii_mm = unique([base.(p.support).a_ref_m*1e3,study.a_ref_mm]);
    rows = cell(1,numel(radii_mm));
    for j = 1:numel(radii_mm)
        p.a_ref_m = radii_mm(j)*1e-3;
        radius_label = strrep(sprintf('%g',radii_mm(j)),'.','p');
        stem = sprintf('pulse_%s_aref_%smm',p.support,radius_label);
        [T,row] = run_case(f,p,study,edges_s,levels_V,out,stem);
        row = add_recovery(row,T,edges_s(3),study.pulse_band_fraction);
        row.simulation_end_ms = 1e3*edges_s(4);
        row.displacement_at_end_um = T.displacement_um(end);
        row.mu_Pa = p.mu_Pa;
        row.eta_Pa_s = p.eta_Pa_s;
        row.c2_Pa = p.c2_Pa;
        row.post_off_ringing_ms = ringing_duration(T,edges_s(3), ...
            study.pulse_band_fraction*row.peak_um,p.abs_tol);
        rows{j} = row;
    end
    support_rows{i} = vertcat(rows{:});
end
P = struct2table(vertcat(support_rows{:}));
writetable(P,fullfile(out,'pulse_summary.csv'));
plot_results(out,fullfile(root,'figure'));
disp(P(:,{'support','a_ref_mm','peak_um','post_off_min_um', ...
    'post_off_ringing_ms','min_hoop_stress_Pa'}));
end

function duration_ms = ringing_duration(T,off_s,threshold_um,velocity_tol_m_s)
% Use the historical adjacent-extrema definition of significant ringing.
after = T.time_s>=off_s;
t = T.time_s(after);
u_um = T.displacement_um(after);
du_um = diff(u_um);
% Ignore slopes within the solver tolerance before looking for reversals.
active = find(abs(du_um)>1e6*velocity_tol_m_s*diff(t));
turns = find(du_um(active(1:end-1)).*du_um(active(2:end))<0);
extrema = active(turns)+1;
excursions_um = abs(diff(u_um(extrema)))/2;
last = find(excursions_um>threshold_um,1,'last');
duration_ms = 0;
if ~isempty(last)
    duration_ms = 1e3*(t(extrema(last+1))-off_s);
end
end

function p = support_params(base, support)
p = rmfield(base,{'spring','gas'});
p.support = support;
specific = base.(support);
keys = fieldnames(specific);
for j = 1:numel(keys)
    p.(keys{j}) = specific.(keys{j});
end
if strcmp(support,'gas')
    p.c1_Pa_s_m = p.c1_prime_Pa_s_m;
    p.c2_Pa = p.c2_prime_Pa;
    p = rmfield(p,{'c1_prime_Pa_s_m','c2_prime_Pa'});
end
end

function [T, row] = run_case(f,p,study,edges_s,levels_V,out,stem)
% Rebalance each a_ref at zero voltage, then verify each trajectory.
assert(p.a_ref_m>0 && p.a_ref_m<p.r0_m && p.rc_m<p.r0_m);
p = f.initial(p);
T = solver(p,edges_s,levels_V,p.dt_s);
assert(height(T)==round(edges_s(end)/(p.dt_s/p.output_substeps))+1);
assert(all(isfinite(T{:,:}),'all') && all(diff(T.time_s)>0));
assert(all(T.radius_m>p.rc_m) && all(T.thickness_m>0));
assert(all(T.radius_m<p.spring_free_m) || strcmp(p.support,'gas'));
on_s = edges_s(2);
expected_V = p.pulse_V*double(T.time_s>=on_s);
if numel(levels_V)==3
    expected_V(T.time_s>=edges_s(3)) = 0;
end
assert(isequal(T.voltage_V,expected_V));

tight = p;
tight.rel_tol = study.check_rel_tol;
tight.abs_tol = study.check_abs_tol;
D = solver(tight,edges_s,levels_V,p.dt_s/2);
radius_error_m = max(abs(T.radius_m-D.radius_m));
assert(radius_error_m<study.radius_check_limit_m);
pre_drift_m = max(abs(T.radius_m(T.time_s<on_s)-p.r0_m));
assert(pre_drift_m<study.equilibrium_limit_m);
area = (T.radius_m+T.thickness_m).^2-T.radius_m.^2;
area_error = max(abs(area/area(1)-1));
assert(area_error<1e-12);

T.displacement_um = 1e6*(T.radius_m-p.r0_m);
[peak_um, peak] = max(T.displacement_um);
row = struct('support',p.support,'a_ref_mm',1e3*p.a_ref_m, ...
    'damping_N_s_m3',p.c1_Pa_s_m,'ks_N_m3',p.ks_Pa_m,'p_eq_Pa',p.gas_eq_Pa, ...
    'peak_um',peak_um,'peak_time_ms',1e3*T.time_s(peak), ...
    'min_hoop_stress_Pa',min(T.hoop_stress_Pa), ...
    'compressive_stress_ms',1e3*sum(diff(T.time_s).* ...
    (T.hoop_stress_Pa(1:end-1)<0)), ...
    'tight_radius_error_m',radius_error_m,'pre_onset_drift_m',pre_drift_m, ...
    'area_rel_error',area_error);
writetable(T,fullfile(out,[stem '.csv']));
fprintf('Completed %s\n',stem);
end

function row = add_recovery(row,T,off_s,band_fraction)
% Measure pulse recovery relative to the zero-voltage initial radius.
u_um = T.displacement_um;
off = find(T.time_s==off_s,1,'first');
outside = find(abs(u_um(off:end))> ...
    band_fraction*max(abs(u_um)),1,'last');
row.displacement_at_off_um = u_um(off);
row.post_off_min_um = min(u_um(off:end));
row.post_off_settle_1pct_ms = NaN;
if isempty(outside)
    row.post_off_settle_1pct_ms = 0;
elseif off+outside<=height(T)
    row.post_off_settle_1pct_ms = 1e3*(T.time_s(off+outside)-off_s);
end
end
