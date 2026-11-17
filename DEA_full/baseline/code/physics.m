function f = physics()
% Figure 3 final mechanics plus the exterior-field closure from DEA_full.
% Governing balance: r_ddot = (p1-p3)/rho_s; sigma_DEA = sigma_m+sigma_e.
% Only this file owns geometry, material laws, support laws and static preload.
f.diagnostics = @diagnostics;
f.mechanical_rhs = @momentum_rhs;
f.electrical = @electrical;
f.initial = @initial_balance;
f.thickness = @thickness;
f.elastic = @elastic;
end

function [dy, q] = momentum_rhs(y, V_DEA, p)
% ====== State: r, r_dot, epsilon_e,i ======
r = y(1);
r_dot = y(2);
epsilon_e = y(3:end);
if r <= p.rc_m
    error('DEA:Geometry','DEA radius %.9g m reached the conductor.',r);
end

% ====== Eq. (3): incompressible thickness ======
d = thickness(r, p);

% ====== Eqs. (8) and (10): complete elastic stress sigma_e ======
lambda = r/p.a_ref_m;
sigma_e0 = elastic(lambda, p);
strain = r/p.r0_m - 1;          % Distinct from lambda=r/a_ref.
relaxation_sum = sum(p.ke_Pa(:).*(strain-epsilon_e));
sigma_e = sigma_e0 + relaxation_sum;

% ====== Eqs. (7) and (6): Maxwell stress and membrane pressure ======
sigma_m = -p.eps0_F_m*p.epsr*(V_DEA/d)^2;
sigma_DEA = sigma_e + sigma_m;
p3 = sigma_DEA*d/r;

% ====== Eqs. (5a)-(5b): p1 = support_pressure - resistance ======
% Gas uses c1_prime and c2_prime in the manuscript; common code fields below.
resistance = p.c1_Pa_s_m*r_dot + p.c2_Pa;
if strcmp(p.support, 'spring')
    support_pressure = p.ks_Pa_m*(p.spring_free_m-r);
else
    support_pressure = p.gas_eq_Pa*(p.r0_m^2-p.rc_m^2) ...
        /(r^2-p.rc_m^2) - p.pout_Pa;
end

% ====== Eqs. (4) and (1): surface density and radial momentum ======
rho_s = p.rho_kg_m3*d*(2*r+d)/(2*r);
r_ddot = (support_pressure-p3-resistance)/rho_s;

% ====== Eq. (9): dashpot strain evolution ======
epsilon_e_dot = (p.ke_Pa(:)./p.eta_Pa_s(:)).*(strain-epsilon_e);
dy = [r_dot; r_ddot; epsilon_e_dot];

% Existing output-column order is preserved; no additional stress is defined.
q = [d, lambda, sigma_e0, relaxation_sum, sigma_m, rho_s, ...
     support_pressure, p3, resistance, r_ddot];
end

function d = thickness(r, p)
% Eq. (3). Calculates d [m] from r [m] and constant annular film area.
area_term = p.d0_m*(2*p.r0_m+p.d0_m);
d = area_term./(sqrt(r.^2+area_term)+r);
end

function sigma_e0 = elastic(lambda, p)
% Eq. (8). Calculates Ogden stress [Pa]; lambda=r/a_ref, not r/r0.
sigma_e0 = sum(p.mu_Pa(:).*(lambda.^p.alpha(:)-lambda.^(-p.alpha(:))));
end

function p = initial_balance(p)
% ====== Initial-equilibrium paragraph after Eq. (10) ======
% V_DEA(0)=0, epsilon_e,i(0)=0, sigma_e(0)=sigma_e0(0).
% Spring: ks*(r_eq-r0)-c2 = sigma_e0(0)*d0/r0.
% Gas: p_eq-p_out-c2_prime = sigma_e0(0)*d0/r0.
sigma_e0_initial = elastic(p.r0_m/p.a_ref_m, p);
if strcmp(p.support, 'spring')
    p.gas_eq_Pa = NaN; % Not applicable to spring; output schema only.
    if p.req_min_m<=p.r0_m || p.req_max_m<p.req_min_m
        error('DEA:SpringBounds','The spring radius limits must satisfy r0 < min <= max.');
    end
    tolerance = 32*eps(p.req_max_m);
    if ~isfinite(p.req_m) || p.req_m<=p.r0_m ...
            || p.req_m<p.req_min_m-tolerance || p.req_m>p.req_max_m+tolerance
        error('DEA:InvalidSpringRadius','Input r_eq must exceed r0 and lie in [%.9g, %.9g] mm.', ...
            p.req_min_m*1e3,p.req_max_m*1e3);
    end
    p.spring_free_m = p.req_m; % Internal name used by the spring pressure law.
    preload = sigma_e0_initial*p.d0_m/p.r0_m+p.c2_Pa;
    p.ks_Pa_m = preload/(p.req_m-p.r0_m);
    if ~isfinite(p.ks_Pa_m) || p.ks_Pa_m<=0
        error('DEA:SpringPreload','Calculated spring stiffness must be positive.');
    end
else
    p.ks_Pa_m = NaN;
    p.spring_free_m = NaN;
    p.gas_eq_Pa = p.pout_Pa+sigma_e0_initial*p.d0_m/p.r0_m+p.c2_Pa;
end
end

function q = diagnostics(y,q_l,p)
% Output only: evaluate electrical terms and the current mechanical state.
e = electrical(y(1),q_l,p);
[~, mechanical] = momentum_rhs(y,e.V_DEA_V,p);
q = e;
q.thickness_m = mechanical(1);
q.lambda = mechanical(2);
q.sigma_e0_Pa = mechanical(3);
q.sigma_relaxation_Pa = mechanical(4);
q.sigma_e_Pa = mechanical(3)+mechanical(4);
q.sigma_m_Pa = mechanical(5);
q.sigma_DEA_Pa = q.sigma_e_Pa+q.sigma_m_Pa;
q.rho_s_kg_m2 = mechanical(6);
q.support_pressure_Pa = mechanical(7);
q.resistance_Pa = mechanical(9);
q.p1_Pa = mechanical(7)-mechanical(9);
q.p3_Pa = mechanical(8);
q.acceleration_m_s2 = mechanical(10);
end

function q = electrical(r,q_l,p)
% Exterior coaxial charge-to-voltage approximation from DEA_full.
d = thickness(r,p);
R = r+d;
if ~isfinite(R) || r<=p.rc_m || R>=p.rg_m
    error('DEA:Geometry','Coaxial geometry requires rc < r and R < rg (r = %.9g m, R = %.9g m).',r,R);
end
charge_scale = q_l/(2*pi*p.eps0_F_m);
U4 = charge_scale*log(p.rg_m/R);
V_DEA = p.U2_V-U4;                  % Eq. (18).
E_signed = charge_scale/R;          % Cancel the coaxial logarithm analytically.
E4 = abs(E_signed);
% C is the fixed input from params; geometry changes E4, never C.
J_corona = sign(E_signed)*p.current_C*max(E4-p.Eth_V_m,0)^p.current_n;
j_o = p.kappa_S_m*V_DEA/d;          % Eq. (14), signed thin-film conduction.
J_total = J_corona-j_o;            % Legacy netCurrent.
q_dot = -J_total*(2*pi*R);          % Legacy chargeIncrement divided by dt.
j_c = -J_corona;                  % Incoming-direction output only.
j_t = -J_total;
q = struct('U4_V',U4,'V_DEA_V',V_DEA,'E4_V_m',E4, ...
    'jc_A_m2',j_c,'jo_A_m2',j_o,'jt_A_m2',j_t, ...
    'q_dot_C_m_s',q_dot,'surface_radius_m',R, ...
    'E_signed_V_m',E_signed,'J_corona_A_m2',J_corona,'J_total_A_m2',J_total);
end

