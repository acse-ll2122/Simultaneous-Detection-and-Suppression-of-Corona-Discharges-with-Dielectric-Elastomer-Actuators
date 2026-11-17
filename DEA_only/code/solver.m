function T = solver(p, edges_s, levels_V, max_step_s)
% Integrate successive constant-voltage segments on a common output grid.
f = physics();
dt_out_s = p.dt_s/p.output_substeps;
indices = round(edges_s/dt_out_s) + 1;
assert(edges_s(1)==0 && numel(edges_s)==numel(levels_V)+1);
assert(all(diff(indices)>=2) && ...
    all(abs((indices-1)*dt_out_s-edges_s)<1e-12));
t = (0:indices(end)-1)'*dt_out_s;
t(indices) = edges_s;
Y = zeros(numel(t), 2+numel(p.ke_Pa));
Y(1,:) = [p.r0_m, 0, zeros(1,numel(p.ke_Pa))];
V = zeros(size(t));
V(end) = levels_V(end);

% ====== Use exact switch times and carry the continuous state across ======
for j = 1:numel(levels_V)
    first = indices(j);
    last = indices(j+1);
    voltage_V = levels_V(j);
    V(first:last-1) = voltage_V;
    [~, q] = f.rhs(Y(first,:)', voltage_V, p);
    initial_step_s = min(max_step_s, 0.2*q(6)/p.c1_Pa_s_m);
    options = odeset('RelTol',p.rel_tol,'AbsTol',p.abs_tol, ...
        'MaxStep',max_step_s,'InitialStep',initial_step_s);
    [~, states] = ode45(@(~,y) f.rhs(y,voltage_V,p), ...
        t(first:last), Y(first,:)', options);
    Y(first+1:last,:) = states(2:end,:);
end

% ====== Diagnostics use the voltage active at each saved time ======
Q = zeros(numel(t),10);
for k = 1:numel(t)
    [~, Q(k,:)] = f.rhs(Y(k,:)', V(k), p);
end
columns = {'time_s','voltage_V','radius_m','velocity_m_s', ...
    'thickness_m','lambda','elastic_Pa','viscoelastic_Pa','maxwell_Pa', ...
    'areal_density_kg_m2','support_Pa','membrane_Pa','damping_Pa', ...
    'acceleration_m_s2','hoop_stress_Pa'};
dashpot_columns = cellstr(compose('dashpot_strain_%d',1:numel(p.ke_Pa)));
T = array2table([t,V,Y(:,1:2),Q,sum(Q(:,3:5),2),Y(:,3:end)], ...
    'VariableNames',[columns,dashpot_columns]);
end
