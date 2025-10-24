%% Peak with filter
time_start_filter=tic;
freq=20;%cut of frequency of high-pass filter (user chosen)
n_psi=n_input_error-n_w;
%scaling chosen heuristically based on observered magnitudes
Input_Scaling=inv(diag([ones(4,1);10;10;0.1;0.1;1]));%eye(n_psi);% inv(diag([ones(4,1);10;1;1;1;50]));%d=[w;x_r;u], scale inptus in B,D s.t. smilar magnitude
%filter dynamics (with normalized input)
A_psi=-eye(n_psi)*freq;
B_psi=freq*[zeros(n_psi,n_w),eye(n_psi)*Input_Scaling];
C_psi=[zeros(n_w,n_psi);-1*eye(n_psi)];
D_psi=[eye(n_w),zeros(n_w,n_psi);zeros(n_psi,n_w),eye(n_psi)*Input_Scaling];
%D_psi=[eye(n_w),zeros(n_w,n_psi);zeros(n_psi,n_input_error)];
A_combined=[A,zeros(n_x,n_psi);zeros(n_psi,n_x),A_psi];
B_combined=[B_e;B_psi];
C_combined_output=[C,zeros(n_y,n_psi)];
C_combined_input=[zeros(n_input_error,n_x),C_psi];
D_combined_input=D_psi;
%2. LMI
n_lambda=10;
lambda_values = linspace(-2*max(real(eig(A)))/n_lambda, -2*max(real(eig(A_combined))), n_lambda);  %
% Initialize best gamma and lambda
best_gamma_filtered = inf;
best_lambda_filtered = NaN;
best_P_filtered = [];
feasible_found = false;

% Settings for solver
%'solver','sdpt3',...
ops = sdpsettings('verbose', 0);  % Set to 1 if you want solver output
ops.usex0 = 0;
% Loop over lambda values
for lambda = lambda_values
gamma_filtered = sdpvar(1);
P_filtered = sdpvar(n_x+n_psi);
% LMI constraints
con = [P_filtered >= 0];
% Define the middle and side matrices
middle_matrix = [lambda*P_filtered,P_filtered,zeros(n_x+n_psi, n_input_error);
                 P_filtered,zeros(n_x+n_psi),zeros(n_x+n_psi, n_input_error);
                     zeros(n_input_error, n_x+n_psi),zeros(n_input_error, n_x+n_psi), -gamma_filtered*eye(n_input_error)];
side_matrix = [eye(n_x+n_psi),           zeros(n_x+n_psi, n_input_error);
               A_combined,                  B_combined;
               C_combined_input,    D_combined_input];
% Main LMI constraint
con = [con;
       side_matrix' * middle_matrix * side_matrix <= 0];
% Additional LMI
con = [con;
      [sqrt(lambda)*gamma_filtered * eye(n_y), C_combined_output;
       C_combined_output', sqrt(lambda)*P_filtered] >= 0];
% Solve
%ops = sdpsettings('verbose', 0);%'solver', 'mosek', 
diagnostics = optimize(con, gamma_filtered, ops);

% Check if solution is feasible
   if diagnostics.problem ~= 1 && diagnostics.problem ~= 2
        gamma_filtered_val = value(gamma_filtered);
        if gamma_filtered_val < best_gamma_filtered
            best_gamma_filtered = gamma_filtered_val;
            best_lambda_filtered = lambda;
            best_P_filtered = value(P_filtered);
            feasible_found = true;
        end
   else
   end
end
% Report result
fprintf('Filter: gamma: %.6f at lambda = %.6f\n in %.6f seconds', best_gamma_filtered, best_lambda_filtered,toc(time_start_filter));
P_filtered=best_P_filtered;
lambda_filtered=best_lambda;
gamma_filtered=best_gamma_filtered;
time_filter=toc(time_start_filter);