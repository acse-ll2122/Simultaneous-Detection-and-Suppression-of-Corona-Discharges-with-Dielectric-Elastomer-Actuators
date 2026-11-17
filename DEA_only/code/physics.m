function f = physics()
% Corrected 101426 manuscript, Eqs. (1)-(10), in manuscript notation.
% Governing balance: r_ddot = (p1-p3)/rho_s; sigma_DEA = sigma_m+sigma_e.
% Only this file owns geometry, material laws, support laws and static preload.
f.rhs = @momentum_rhs;
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
    error('DEA radius reached the conductor; model invalid.');
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
rho_s = p.rho_kg_m3*((r+d)^2-r^2)/(2*r);
% Equivalent to (p1-p3)/rho_s; keep the validated arithmetic grouping.
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
d = sqrt(r.^2+2*p.r0_m*p.d0_m+p.d0_m^2) - r;
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
    assert(p.req_min_m>p.r0_m && p.req_max_m>=p.req_min_m);
    tolerance = 32*eps(p.req_max_m);
    assert(isfinite(p.req_m) && p.req_m>p.r0_m ...
        && p.req_m>=p.req_min_m-tolerance && p.req_m<=p.req_max_m+tolerance, ...
        'Input r_eq must exceed r0 and lie in [%.9g, %.9g] mm.', ...
        p.req_min_m*1e3,p.req_max_m*1e3);
    p.spring_free_m = p.req_m; % Internal name used by the spring pressure law.
    preload = sigma_e0_initial*p.d0_m/p.r0_m+p.c2_Pa;
    p.ks_Pa_m = preload/(p.req_m-p.r0_m);
    assert(isfinite(p.ks_Pa_m) && p.ks_Pa_m>0,'Calculated ks must be positive.');
else
    p.ks_Pa_m = NaN;
    p.spring_free_m = NaN;
    p.gas_eq_Pa = p.pout_Pa+sigma_e0_initial*p.d0_m/p.r0_m+p.c2_Pa;
end
end
