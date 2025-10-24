import casadi.*

z_ref = z_max*ones(1,H+1);%reference

if size(B_psi,2)-(n_r+n_u+n_w)~=0%just checing to be sure
   error('change in code, dimension inconsistent') 
end
A_MPC_c=[A_r,zeros(n_r,n_psi+1);...
        B_psi(:,n_w+1:n_w+n_r),A_psi,zeros(n_psi,1);...
        zeros(1,n_r+n_psi),-lambda_filtered];
B_MPC_c=[B_r,zeros(n_r,1);...
        B_psi(:,n_w+n_r+1:n_w+n_r+n_u,:),zeros(n_psi,1);...
        zeros(1,n_u),gamma_filtered];
[A_MPC_d,B_MPC_d]=c2d(A_MPC_c,B_MPC_c,T_s);

%% Decision variables
x_r = SX.sym('x_r', n_r, H+1);
psi = SX.sym('psi', n_psi, H+1);
delta = SX.sym('delta', 1, H+1);
delta_z = SX.sym('delta_z', 1, H+1);
delta_z_square = SX.sym('delta_z_square', 1, H+1);
cost_z = SX.sym('cost_z', 1, H+1);
cost_u = SX.sym('cost_u', 1, H);
u = SX.sym('u', n_u, H);
w = SX.sym('w', 1, H);
x_MPC = [x_r; psi; delta];
u_MPC = [u; w];
%% Dynamics loop
g_dyn = [];
X_next = A_MPC_d * x_MPC(:,1:H) + B_MPC_d * u_MPC;
g_dyn = x_MPC(:,2:H+1) - X_next;  
%compute bound on |r_psi|^2
stacked_input = [zeros(n_w,H); x_r(:,1:H); u]; 
r_psi = C_psi*psi(:,1:H) + D_psi*stacked_input;
g_r_psi = w - sum(r_psi.^2,1) - w_max^2;

z_r = C_r*x_r(:,:);
%ensure that delta_z 
g_delta_z = delta*lambda_filtered*gamma_filtered - delta_z.^2;
g_z_ub = z_r + delta_z;
g_cost_z = [z_r - z_ref - cost_z; -(z_r - z_ref) - cost_z];
g_cost_u = [u - cost_u; -u - cost_u];

%% Combine constraints
g = [x_MPC(:,1) - [x_0r; zeros(n_psi,1); delta_0];  % initial condition
    g_dyn(:);                                      % flatten dynamics
    g_r_psi(:);                                    % flatten r_psi
    g_delta_z(:);                                  % flatten delta_z
    g_z_ub(:);                                     % flatten z bounds
    g_cost_z(:);                                   % flatten cost_z
    g_cost_u(:)                                    % flatten cost_u
];

% Lower bounds
lb_g = [ ...
    zeros(n_r+n_psi+1,1);      % initial condition
    zeros(numel(g_dyn),1);     % dynamics equality
    zeros(numel(g_r_psi),1);   % r_psi >=0
    zeros(numel(g_delta_z),1); % delta_z equality
    -inf*ones(numel(g_z_ub),1);   % z upper bound <= z_max
    -inf*ones(numel(g_cost_z),1); % cost_z <=0
    -inf*ones(numel(g_cost_u),1)  % cost_u <=0
];

% Upper bounds
ub_g = [ ...
    zeros(n_r+n_psi+1,1);      % initial condition
    zeros(numel(g_dyn),1);     % dynamics equality
    inf*ones(numel(g_r_psi),1);   % r_psi >=0
    zeros(numel(g_delta_z),1);     % delta_z equality
    z_max*ones(numel(g_z_ub),1);  % z upper bound
    zeros(numel(g_cost_z),1);      % cost_z <=0
    zeros(numel(g_cost_u),1)       % cost_u <=0
];

% Objective
obj = R*sum(cost_u) + sum(cost_z) + sum(delta_z);

vars = vertcat(x_r(:), psi(:), delta(:), delta_z(:), delta_z_square(:), cost_z(:), cost_u(:), u(:), w(:));
nlp = struct('x', vars, 'f', obj, 'g', g);
opts = struct('ipopt', struct('print_level',0,'tol',1e-3,'max_iter',1000,'linear_solver','mumps'));
solver = nlpsol('solver','ipopt',nlp,opts);


x0 = zeros(length(vars),1);%zero warmstart (for fair comparison)

% Lower bounds (enforce non-negativity for some blocks)
lb = -inf(length(vars),1);
sizes = [n_r*(H+1), n_psi*(H+1), H+1, H+1, H+1, H+1, H, n_u*H, H];
start_idx = [1, cumsum(sizes(1:end-1))+1]; end_idx = cumsum(sizes);
pos_blocks = [3,4,6,7,9]; % delta, delta_z, cost_z, cost_u, w
for b = pos_blocks
    lb(start_idx(b):end_idx(b)) = 0;
end

ub = inf(length(vars),1);

%% Solve NLP
sol = solver('x0', x0, 'lbx', lb, 'ubx', ub, 'lbg', lb_g, 'ubg', ub_g);
stats=solver.stats;
stats.t_wall_total

x_opt = full(sol.x);
fprintf('Solve time: %.6f [s]', stats.t_wall_total);

%% Unpack variables
offset = 0;
x_r_opt = reshape(x_opt(offset + (1:n_r*(H+1))), n_r, H+1); offset = offset + n_r*(H+1);
psi_opt = reshape(x_opt(offset + (1:n_psi*(H+1))), n_psi, H+1); offset = offset + n_psi*(H+1);
delta_opt = reshape(x_opt(offset + (1:H+1)), 1, H+1); offset = offset + H+1;
delta_z_opt = reshape(x_opt(offset + (1:H+1)), 1, H+1); offset = offset + H+1;
delta_z_square_opt = reshape(x_opt(offset + (1:H+1)), 1, H+1); offset = offset + H+1;
cost_z_opt = reshape(x_opt(offset + (1:H+1)), 1, H+1); offset = offset + H+1;
cost_u_opt = reshape(x_opt(offset + (1:H)), 1, H); offset = offset + H;
u_opt = reshape(x_opt(offset + (1:n_u*H)), n_u, H); offset = offset + n_u*H;
w_opt = reshape(x_opt(offset + (1:H)), 1, H);

x_MPC_opt = [x_r_opt; psi_opt; delta_opt]; 
u_MPC_opt = [u_opt; w_opt];

%% Compute r_psi_opt and z_opt 
stacked_input = [zeros(n_w,H); x_r_opt(:,1:H); u_opt];  
r_psi_opt = C_psi*psi_opt(:,1:H) + D_psi*stacked_input; 
z_opt = C_r*x_r_opt(:,1:H);  



%% Simulation / plotting
x_MPC_sim = zeros(size(x_MPC_opt));
dyn_err = zeros(1,H);
for k=1:H
    dyn_err(k) = norm(A_MPC_d*x_MPC_opt(:,k) + B_MPC_d*u_MPC_opt(:,k) - x_MPC_opt(:,k+1));
    x_MPC_sim(:,k+1) = A_MPC_d*x_MPC_sim(:,k) + B_MPC_d*u_MPC_opt(:,k);
end

[A_d,B_d] = c2d(A,B,T_s);
x_full = zeros(n_x,H+1); z_exact = zeros(n_y,H);
for k=1:H
    x_full(:,k+1) = A_d*x_full(:,k) + B_d*u_opt(:,k);
    z_exact(:,k) = C*x_full(:,k);
end

% Plot
t = (0:H-1)*T_s;
figure(2); hold on
plot( t, z_opt-delta_z_opt(1:end-1), 'b',t, z_exact, 'r--', [0,H*T_s],[z_max,z_max],'k--',t, z_opt+delta_z_opt(1:end-1), 'b','linewidth',2);
fill([t, fliplr(t)], [z_opt+delta_z_opt(1:end-1), fliplr(z_opt-delta_z_opt(1:end-1))], [0.9 0.9 0.9], 'EdgeColor','none');
plot(t, z_exact, 'r--','linewidth',2)
%plot(,'LineWidth',2);
axis([0,H*T_s,min(z_opt-delta_z_opt(1:end-1))*1.2,z_max*1.05]);
legend({'Robust prediction','Full-order simulation','Constraint'},'location','southeast');
xlabel('Time $t$ [s]','interpreter','latex'); ylabel('Position $z [m]$','interpreter','latex');
set(gca,'fontname','Arial','fontsize',16,'YTick',[0 0.5 1]);
print('OCP_ROM','-depsc');
