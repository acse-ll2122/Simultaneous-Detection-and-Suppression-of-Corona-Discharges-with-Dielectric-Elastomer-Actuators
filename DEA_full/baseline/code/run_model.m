function run_model(base,out,mode,baseline_out,figure_out)
% Compute baseline supports or the spring-radius study with shared physics.
root = fileparts(fileparts(mfilename('fullpath')));
if nargin<1,base=params();end
if nargin<2,out=fullfile(root,'results');end
if nargin<3,mode='baseline';end
if nargin<4,baseline_out=fullfile(root,'results');end
if nargin<5,figure_out=fullfile(root,'figure');end
assert(ismember(mode,{'baseline','sensitive'}),'Unknown study mode.');
base.r0_m = base.rc_m+base.t1_m;
validateattributes(base.duration_s,{'numeric'},{'scalar','finite','positive'},mfilename,'duration_s');
validateattributes(base.spring.req_m,{'numeric'}, ...
    {'real','scalar','finite','positive'},mfilename,'spring.req_m');
validateattributes(base.spring_sweep_mm,{'numeric'},{'real','finite'},mfilename,'spring_sweep_mm');
validateattributes(base.max_step_s,{'numeric'},{'scalar','finite','positive'},mfilename,'max_step_s');
validateattributes(base.check_factor,{'numeric'}, ...
    {'scalar','finite','positive','<',1},mfilename,'check_factor');
validateattributes(base.integration_attempts,{'numeric'}, ...
    {'scalar','finite','integer','positive'},mfilename,'integration_attempts');
validateattributes(base.refinement_attempts,{'numeric'}, ...
    {'scalar','finite','integer','positive'},mfilename,'refinement_attempts');
if ~isfolder(out)
    mkdir(out);
end
sweep = strcmp(mode,'sensitive');
if sweep
    if ~isfile(fullfile(baseline_out,'run.mat'))
        run_model(base,baseline_out,'baseline',baseline_out,fullfile(root,'figure'));
    end
    reference = load_results(baseline_out);
    assert(isequaln(reference.base,base), ...
        'Baseline parameters changed. Run baseline/code/run_all.m before the sweep.');
    spring_mm = unique([base.spring.req_m*1e3,base.spring_sweep_mm(:)']);
else
    spring_mm = base.spring.req_m*1e3;
end
is_baseline = spring_mm==base.spring.req_m*1e3;
cases = cell(1,numel(spring_mm)+~sweep);
names = cell(size(cases));rows = cell(size(cases));skipped = cell(0,4);
for k = 1:numel(spring_mm)
    p = support_params(base,'spring');
    p.req_m = spring_mm(k)*1e-3;
    radius_label = strrep(sprintf('%.12g',spring_mm(k)),'.','p');
    names{k} = ['Spring_req_' radius_label 'mm'];
    if any(strcmp(names{k},names(1:k-1)))
        error('DEA:DuplicateCase','Two sweep radii produce the same label %s.',names{k});
    end
    if sweep && is_baseline(k)
        index = find(strcmp(reference.names,names{k}),1);
        assert(~isempty(index),'Baseline spring trajectory is missing.');
        cases{k} = reference.cases{index};
        row = table2struct(reference.summary(strcmp(reference.summary.case_name,names{k}),:));
        rows{k} = row;
        continue
    end
    try
        [cases{k},rows{k}] = run_case(p,names{k});
    catch failure
        case_failure = ismember(failure.identifier,{'DEA:InvalidSpringRadius', ...
            'DEA:SpringPreload','DEA:Geometry','DEA:IntegrationIncomplete'});
        if is_baseline(k) || ~case_failure,rethrow(failure);end
        skipped(end+1,:) = {names{k},spring_mm(k),failure.identifier,failure.message}; %#ok<AGROW>
        warning('DEA:SkippedSweepCase','Skipped %s. %s',names{k},failure.message);
    end
end
if ~sweep
    names{end} = 'Gas';
    [cases{end},rows{end}] = run_case(support_params(base,'gas'),names{end});
end
completed = ~cellfun(@isempty,cases);
spring_mm = spring_mm(completed(1:numel(spring_mm)));
baseline_spring = find(spring_mm==base.spring.req_m*1e3,1);
cases = cases(completed);names = names(completed);
computed_row = find(completed & ~is_baseline_padded(is_baseline,~sweep),1);
if isempty(computed_row),computed_row=find(completed,1);end
row_order = fieldnames(rows{computed_row});
for k=find(completed),rows{k}=orderfields(rows{k},row_order);end
summary = struct2table(vertcat(rows{completed}));
summary = movevars(summary,{'duration_s','samples','radius_change_end_um', ...
    'V_DEA_end_V','E4_end_MV_m','J_tot_end_A_m2'},'After','r_eq_mm');
shared_index = [];
shared_baseline = '';
if sweep
    shared_index = baseline_spring;
    shared_baseline = '../../../baseline/results/run.mat';
end
write_results(cases,names,summary,skipped,out,shared_index);
if sweep,cases{baseline_spring}=[];end % Resolve the shared trajectory on load.
save(fullfile(out,'run.mat'),'base','cases','names','summary','spring_mm','baseline_spring','shared_baseline');
parameter_file = fopen(fullfile(out,'parameters.json'),'w');
if parameter_file<0,error('DEA:ParameterWrite','Cannot write parameters in %s.',out);end
close_parameters = onCleanup(@()fclose(parameter_file));
fprintf(parameter_file,'%s\n',jsonencode(base));clear close_parameters
plot_reviewed(out,figure_out,mode);
disp(summary(:,{'case_name','duration_s','samples','radius_change_end_um', ...
    'V_DEA_end_V','E4_end_MV_m','J_tot_end_A_m2'}));
end

function flags = is_baseline_padded(flags,has_gas)
if has_gas,flags=[flags false];end
end

function p = support_params(base,support)
p = rmfield(base,{'spring','gas','spring_sweep_mm'});
p.support = support;
specific = base.(support);
keys = fieldnames(specific);
for k = 1:numel(keys)
    p.(keys{k}) = specific.(keys{k});
end
if strcmp(support,'gas')
    p.c1_Pa_s_m = p.c1_prime_Pa_s_m;
    p.c2_Pa = p.c2_prime_Pa;
end
if ~(p.a_ref_m>0 && p.a_ref_m<p.r0_m && p.r0_m>p.rc_m && p.d0_m>0)
    error('DEA:InitialGeometry','%s requires 0 < a_ref < r0, rc < r0 and d0 > 0.',support);
end
if ~(p.U2_V<0 && p.rg_m>p.r0_m+p.d0_m)
    error('DEA:ElectricalInputs','%s requires a negative supply and rg > r0 + d0.',support);
end
validateattributes(p.current_C,{'numeric'},{'real','scalar','finite','nonnegative'});
validateattributes(p.kappa_S_m,{'numeric'},{'real','scalar','finite','nonnegative'});
validateattributes(p.Eth_V_m,{'numeric'},{'real','scalar','finite','positive'});
validateattributes(p.current_n,{'numeric'},{'real','scalar','finite','>=',1});
end

function [T,row] = run_case(p,name)
[T,p,sol] = solver(p);
for comparison = 1:p.refinement_attempts
    tight = p;
    keys = {'rel_tol','abs_radius_m','abs_velocity_m_s','abs_strain','abs_charge_C_m'};
    for k = 1:numel(keys)
        tight.(keys{k}) = p.(keys{k})*p.check_factor;
    end
    tight.max_step_s = p.max_step_s/2;
    [fine_T,fine_p,refined] = solver(tight);
    [radius_error_um,voltage_error_V,field_error_V_m,current_error_A_m2] = ...
        refinement_errors(sol,refined,p);
    refinement_passed = radius_error_um<p.check_radius_um && voltage_error_V<p.check_voltage_V ...
        && field_error_V_m<p.check_field_V_m && current_error_A_m2<p.check_current_A_m2;
    if refinement_passed || comparison==p.refinement_attempts
        break
    end
    T = fine_T;
    p = fine_p;
    sol = refined;
end
if ~refinement_passed
    warning('DEA:Refinement','%s reached %.6g s but did not meet the refinement tolerances.', ...
        name,T.time_s(end));
end
% ====== Long-term diagnostics ======
T.delta_radius_um = 1e6*(T.radius_m-p.r0_m);
[peak_um,peak] = max(T.delta_radius_um);
tail_window_s = min(1,p.duration_s/5);
tail_start_s = p.duration_s-tail_window_s;
tail = deval(sol,tail_start_s);
tail_radius_m = [tail(1); T.radius_m(T.time_s>tail_start_s)];
f = physics();
tail_electrical = f.electrical(tail(1),tail(end),p);
tail_voltage_V = [tail_electrical.V_DEA_V; T.V_DEA_V(T.time_s>tail_start_s)];
threshold_time_s = NaN;
below = find(T.E4_V_m<=p.Eth_V_m,1,'first');
if ~isempty(below)
    threshold_time_s = T.time_s(below);
end
row = struct('case_name',name,'r_eq_mm',p.spring_free_m*1e3, ...
    'ks_N_m3',p.ks_Pa_m,'p_eq_Pa',p.gas_eq_Pa, ...
    'mu_Pa',p.mu_Pa,'a_ref_mm',1e3*p.a_ref_m,'damping_N_s_m3',p.c1_Pa_s_m, ...
    'supply_V',p.U2_V,'E_th_V_m',p.Eth_V_m,'current_n',p.current_n, ...
    'current_C',p.current_C,'kappa_S_m',p.kappa_S_m, ...
    'duration_s',T.time_s(end),'samples',height(T),'peak_delta_radius_um',peak_um, ...
    'peak_time_s',T.time_s(peak),'radius_change_end_um',T.delta_radius_um(end), ...
    'radius_span_tail_um',1e6*(max(tail_radius_m)-min(tail_radius_m)), ...
    'V_DEA_end_V',T.V_DEA_V(end), ...
    'V_DEA_span_tail_V',max(tail_voltage_V)-min(tail_voltage_V), ...
    'E4_min_MV_m',min(T.E4_V_m)*1e-6,'E4_end_MV_m',T.E4_V_m(end)*1e-6, ...
    'E4_end_minus_Eth_V_m',T.E4_V_m(end)-p.Eth_V_m, ...
    'first_below_threshold_sample_s',threshold_time_s, ...
    'jc_end_A_m2',T.jc_A_m2(end),'jo_end_A_m2',T.jo_A_m2(end), ...
    'J_tot_end_A_m2',T.J_tot_A_m2(end), ...
    'min_hoop_stress_Pa',min(T.sigma_DEA_Pa), ...
    'solver_max_step_s',p.max_step_s,'solver_rel_tol',p.rel_tol, ...
    'refinement_comparisons',comparison,'refinement_passed',refinement_passed, ...
    'refine_radius_um',radius_error_um,'refine_V',voltage_error_V, ...
    'refine_E4_V_m',field_error_V_m,'refine_J_tot_A_m2',current_error_A_m2);
fprintf('Completed %s, %.3f s, %d accepted samples.\n',name,T.time_s(end),height(T));
end

function [radius_um,voltage_V,field_V_m,current_A_m2] = refinement_errors(sol,refined,p)
f = physics();
time_s = union(sol.x,refined.x);
coarse = deval(sol,time_s);
fine = deval(refined,time_s);
radius_um = 1e6*max(abs(coarse(1,:)-fine(1,:)));
voltage_V = 0;
field_V_m = 0;
current_A_m2 = 0;
for k = 1:numel(time_s)
    q = f.electrical(coarse(1,k),coarse(end,k),p);
    q_fine = f.electrical(fine(1,k),fine(end,k),p);
    voltage_V = max(voltage_V,abs(q.V_DEA_V-q_fine.V_DEA_V));
    field_V_m = max(field_V_m,abs(q.E4_V_m-q_fine.E4_V_m));
    current_A_m2 = max(current_A_m2,abs(q.jt_A_m2-q_fine.jt_A_m2));
end
end

function write_results(cases,names,summary,skipped,out,shared_index)
% CSV files are primary trajectories; the workbook is a convenient export.
writetable(summary,fullfile(out,'summary.csv'));
staged = [tempname '.xlsx']; % Keep temporary workbook paths below Windows length limits.
cleanup = onCleanup(@()remove_staged(staged));
writetable(summary,staged,'Sheet','Summary');
for k = 1:numel(cases)
    if ismember(k,shared_index),continue;end
    writetable(cases{k},fullfile(out,[names{k} '.csv']));
    writetable(cases{k},staged,'Sheet',names{k});
end
if ~isempty(shared_index)
    shared = table(names(shared_index)',{['../../../baseline/results/' names{shared_index} '.csv']}, ...
        'VariableNames',{'case_name','result_file'});
    writetable(shared,fullfile(out,'shared_baseline.csv'));
    writetable(shared,staged,'Sheet','Shared_baseline');
end
if ~isempty(skipped)
    skipped_cases = cell2table(skipped,'VariableNames',{'case_name','r_eq_mm','error_id','reason'});
    writetable(skipped_cases,staged,'Sheet','Skipped');
end
[ok,message] = movefile(staged,fullfile(out,'coupled_results.xlsx'),'f');
if ~ok,error('DEA:WorkbookWrite','Cannot replace workbook. %s',message);end
clear cleanup
end

function remove_staged(path)
if isfile(path),delete(path);end
end

