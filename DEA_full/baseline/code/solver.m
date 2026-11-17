function [T,p,sol] = solver(p)
% Joint integration of radius, velocity, dashpot strain and signed line charge.
f = physics();
p = f.initial(p);
R0 = p.r0_m+p.d0_m;
ql0_C_m = 2*pi*p.eps0_F_m*p.U2_V/log(p.rg_m/R0);
strain0 = zeros(numel(p.ke_Pa),1);
y0 = [p.r0_m;0;strain0;ql0_C_m];
abs_tol = [p.abs_radius_m;p.abs_velocity_m_s; ...
    p.abs_strain+strain0;p.abs_charge_C_m];
options = odeset('RelTol',p.rel_tol,'AbsTol',abs_tol);
max_step_s = p.max_step_s;
for attempt = 1:p.integration_attempts
    try
        sol = ode15s(@(~,y) coupled_rhs(y,p,f),[0 p.duration_s],y0, ...
            odeset(options,'MaxStep',max_step_s));
        if abs(sol.x(end)-p.duration_s)>8*eps(p.duration_s)
            error('DEA:IntegrationIncomplete','Integration ended at %.9g of %.9g s.', ...
                sol.x(end),p.duration_s);
        end
        break
    catch failure
        recoverable = ismember(failure.identifier,{'DEA:Geometry','DEA:IntegrationIncomplete'});
        if ~recoverable || attempt==p.integration_attempts
            rethrow(failure);
        end
        max_step_s = max_step_s/10;
        fprintf('Retrying %s integration with MaxStep = %.3g s (%s).\n', ...
            p.support,max_step_s,failure.identifier);
    end
end
p.max_step_s = max_step_s;
% ====== Evaluate every accepted time step ======
diagnostics = {'U4_V','V_DEA_V','E4_V_m','jc_A_m2','jo_A_m2','jt_A_m2', ...
    'q_dot_C_m_s','surface_radius_m','thickness_m','lambda', ...
    'sigma_e0_Pa','sigma_relaxation_Pa','sigma_e_Pa','sigma_m_Pa','sigma_DEA_Pa', ...
    'rho_s_kg_m2','support_pressure_Pa','resistance_Pa','p1_Pa','p3_Pa','acceleration_m_s2', ...
    'E_signed_V_m','J_corona_A_m2','J_total_A_m2'};
Q = zeros(numel(sol.x),numel(diagnostics));
for i = 1:numel(sol.x)
    q = f.diagnostics(sol.y(1:end-1,i),sol.y(end,i),p);
    for j = 1:numel(diagnostics)
        Q(i,j) = q.(diagnostics{j});
    end
end
diagnostics{6} = 'J_tot_A_m2'; % Inward-positive manuscript convention.
strain_names = cellstr(compose('dashpot_strain_%d',1:numel(p.ke_Pa)));
names = [{'time_s','radius_m','velocity_m_s'},strain_names,{'ql_C_m'},diagnostics];
data = [sol.x',sol.y',Q];
if any(~isfinite(data(:))) || any(diff(sol.x)<=0)
    error('DEA:InvalidOutput','%s integration produced nonfinite data or unordered times.',p.support);
end
T = array2table(data,'VariableNames',names);
end

function dy = coupled_rhs(y,p,f)
% Use the same instantaneous geometry and charge for the two balances.
q = f.electrical(y(1),y(end),p);
dy = [f.mechanical_rhs(y(1:end-1),q.V_DEA_V,p);q.q_dot_C_m_s];
end
