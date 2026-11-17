function [p, study] = params()
% Physical inputs transcribed from the supplied params.m, in SI units.

% ====== Spring ======
p.spring.mu_Pa = 1e4;
p.spring.a_ref_m = 1.2e-3;
p.spring.req_m = 3e-3;
p.spring.c1_Pa_s_m = 1000;
p.spring.c2_Pa = 0;
p.spring.req_min_m = 2e-3;
p.spring.req_max_m = 10e-3;

% ====== Gas ======
p.gas.mu_Pa = 1e4;
p.gas.a_ref_m = 1.2e-3;
p.gas.c1_prime_Pa_s_m = 1000;
p.gas.c2_prime_Pa = 0;

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

% ====== Drive, sweep and integration ======
p.pulse_V = -1500;
p.dt_s = 0.025e-3;
p.output_substeps = 5;
p.rel_tol = 1e-6;
p.abs_tol = 1e-11;
study.on_s = 10e-3;
study.pulse_off_s = 20e-3;
study.end_s.spring = 800e-3;
study.end_s.gas = 200e-3;
study.a_ref_mm = [1 1.2 1.4];
study.pulse_band_fraction = 0.01;

% ====== Verification ======
study.check_rel_tol = 1e-8;
study.check_abs_tol = 1e-13;
study.equilibrium_limit_m = 1e-11;
study.radius_check_limit_m = 1e-8;
end
