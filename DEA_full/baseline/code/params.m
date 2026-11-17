function p = params()
% Figure 3 final baseline mechanics and DEA_full corona inputs. SI units.

% ====== Spring DEA and corona ======
p.spring.mu_Pa = 1e4;
p.spring.a_ref_m = 1.2e-3;
p.spring.req_m = 3e-3;
p.spring.c1_Pa_s_m = 1000;
p.spring.c2_Pa = 0;
p.spring.req_min_m = 2e-3;
p.spring.req_max_m = 10e-3;
p.spring.current_C = 6e-15;  % Fixed in time and radius; n=2 gives A/V^2.
p.spring.kappa_S_m = 3.23e-14; % Foo et al. 2012, Sec. V; VHB 4905 low-field reference.

% ====== Gas DEA and corona ======
p.gas.mu_Pa = 1e4;
p.gas.a_ref_m = 1.2e-3;
p.gas.c1_prime_Pa_s_m = 1000;
p.gas.c2_prime_Pa = 0;
p.gas.current_C = 6e-15; % Same fixed C as Spring.
p.gas.kappa_S_m = 3.23e-14; % Same membrane conductivity as Spring.
p.gas.duration_s = 1;

% ====== Common geometry and material ======
p.rc_m = 1e-3;
p.t1_m = 0.5e-3;
p.r0_m = p.rc_m + p.t1_m;
p.d0_m = 0.1e-3;
p.alpha = 2;
p.ke_Pa = 0.48e6;
p.eta_Pa_s = 4400;
p.eps0_F_m = 8.8541878128e-12;
p.epsr = 3.92;
p.rho_kg_m3 = 1000;
p.pout_Pa = 101325;

% ====== Common corona and supply ======
p.U2_V = -24e3;
p.rg_m = 0.2;
p.Eth_V_m = 3e6;
p.current_n = 2;

% ====== Study and integration ======
p.spring_sweep_mm = [2 2.5 3.5];  % Also run the baseline req_m.
p.duration_s = 5; % Spring study horizon; Gas overrides to 1 s.
p.max_step_s = 0.01;
p.rel_tol = 1e-7;
p.abs_radius_m = 1e-12;
p.abs_velocity_m_s = 1e-8;
p.abs_strain = 1e-10;
p.abs_charge_C_m = 1e-16;
p.check_factor = 0.1;
p.integration_attempts = 3;    % Retry an incomplete solve at one-tenth the step.
p.refinement_attempts = 3;     % Maximum coarse versus fine comparisons.
p.check_radius_um = 0.01;
p.check_voltage_V = 0.1;
p.check_field_V_m = 30;
p.check_current_A_m2 = 1e-6;
end
