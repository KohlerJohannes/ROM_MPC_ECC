%% Compute peak-to-peak
time_start_peak=tic;
% Define range of lambda
n_lambda=10;
lambda_values = linspace(-2*max(real(eig(A)))/n_lambda, -2*max(real(eig(A))), n_lambda);  %
best_gamma = inf;
best_lambda = NaN;
best_P = [];
ops = sdpsettings('verbose', 0); 
for lambda = lambda_values
    gamma = sdpvar(1);
    P = sdpvar(n_x);
    % LMI constraints
    con = [P >= 0];
    % Define the middle and side matrices
    middle_matrix = [lambda*P,         P,              zeros(n_x, n_input_error);
                     P,               zeros(n_x),      zeros(n_x, n_input_error);
                     zeros(n_input_error, n_x),  zeros(n_input_error, n_x), -gamma*eye(n_input_error)];

    side_matrix = [eye(n_x),           zeros(n_x, n_input_error);
                   A,                  B_e;
                   zeros(n_input_error, n_x),    eye(n_input_error)];
    % Main LMI constraint
    con = [con;
           side_matrix' * middle_matrix * side_matrix <= 0];
    % Additional LMI
    con = [con;
           [gamma * eye(n_y), C;
            C',              lambda * P] >= 0];
    % Solve
    diagnostics = optimize(con, gamma, ops);
   % Check if solution is feasible
   if diagnostics.problem ~= 1 && diagnostics.problem ~= 2
        gamma_val = value(gamma);
        if gamma_val < best_gamma
            best_gamma = gamma_val;
            best_lambda = lambda;
            best_P = value(P);
        end
    end
end
% Report result
fprintf('Peak: %.6f at lambda = %.6f in %.6f seconds\n', best_gamma, best_lambda, toc(time_start_peak));
P_peak=best_P;
lambda_peak=best_lambda;
gamma_peak=best_gamma;
time_peak=toc(time_start_peak);