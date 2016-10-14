clear all;
close all;

num_q = 3;
num_cstate = 4;
gamma = 0.999;

SimTime = 200;

x_sim = zeros(num_cstate, SimTime+1);
belief_sim = zeros(num_q, SimTime+1);
control = zeros(SimTime+1, 1);

MatlabSimfile_Discrete = fopen('OutFiles/MatlabSimfile_Discrete.txt', 'w');

x_sim(:, 1) = [-330; 30; -360; 30];
belief_sim(:, 1) = [1/3; 1/3; 1/3];
q_sim = 3;

horizon = 30;

deltaT = 0.2;
mReactionDist = 25;

alpha = 0; % for objective cost
beta = 1; % for u cost
lambda = 20; % for KL-div
eta = 1e8; % for constraints

epsilon = 1e-2;
old_horizon = horizon;

for k = 1:SimTime
     k
    if mod(k, 10) ==0
        k
    end
    horizon = max(5, min(horizon, floor(abs(min(x_sim(1, k), x_sim(3, k)))/30/deltaT)))
    if k == 1
        params = initialize_params_New(q_sim, x_sim(:, k), mReactionDist, deltaT, horizon);
    else
        params = initParamsFromOldParams(params_out, deltaT, num_q, num_cstate, ...
                                         old_horizon, horizon);
    end
    old_horizon = horizon;
  
    state_cov_inv = initialize_cov_inv(params, deltaT, num_cstate, num_q, horizon);

    % setting up the mpc problem
    funcs.objective = @(vars)costFn_Ipopt(vars, state_cov_inv, belief_sim(:, k), ...
                                          alpha, beta, lambda, eta, ...
                                          num_cstate, num_q, horizon);
    funcs.gradient = @(vars)costGradient_Ipopt(vars, state_cov_inv, belief_sim(:, k), ...
                                               alpha, beta, lambda, eta, ...
                                               num_cstate, num_q, horizon);

    options.lb = [-Inf * ones(num_cstate * horizon * num_q, 1); -ones(num_q * horizon, 1)];
    options.ub = [Inf * ones(num_cstate * horizon * num_q, 1); ones(num_q * horizon, 1)];

    funcs.constraints = @(vars)constrantsFn_Ipopt(vars, x_sim(:, k), epsilon, ...
                                      deltaT, num_cstate, num_q, horizon, horizon);
    funcs.jacobian = @(x) constrantsJacobian_Ipopt(deltaT, num_cstate, num_q, horizon);
    funcs.jacobianstructure = @() constrantsJacobian_Ipopt(deltaT, num_cstate, num_q, horizon);

    num_ceq = num_q-1 + horizon * num_q * num_cstate;
    options.cl = zeros(num_ceq, 1);   % Lower bounds on the constraint functions.
    options.cu = zeros(num_ceq, 1);   % Upper bounds on the constraint functions.

    options.ipopt.hessian_approximation = 'limited-memory';
    options.ipopt.mu_strategy = 'adaptive';
    options.ipopt.linear_solver = 'ma97';
    options.ipopt.print_level = 3;
    options.ipopt.tol         = 1e-5;
    options.ipopt.max_iter    = 100;

    % Run IPOPT.
    disp('Solving optimization problem ...')
    
    [params_out, info] = ipopt(params,funcs,options);
    [x_params, u_params] = parseParams(params_out, num_q, num_cstate, horizon);

    sigma = u_params(1);

    % apply output to the system
    control(k) = sigma;
    str = sprintf('q = %d, x = (%.2f, %.2f, %.2f, %.2f), belief = (%.2f, %.2f, %.2f), sigma = %f', ...
                  q_sim, x_sim(1, k), x_sim(2, k), x_sim(3, k), x_sim(4, k), ...
                  belief_sim(1, k), belief_sim(2, k), belief_sim(3, k), control(k));
    disp(str);
    fprintf(MatlabSimfile_Discrete, ...
        '%d  %f  %f  %f  %f  %f  %f  %f  %f\n', ...
        q_sim, x_sim(1, k), x_sim(2, k), x_sim(3, k), x_sim(4, k), ...
        belief_sim(1, k), belief_sim(2, k), belief_sim(3, k), control(k));
    
    if (IsSimEnded(q_sim, x_sim(:, k)))
        break;
    end
    
    [x_next, b_next] = NextStateAndBelief(q_sim, x_sim(:, k), belief_sim(:, k), control(k), deltaT);
    
    x_sim(:, k+1) = x_next;
    belief_sim(:, k+1) = b_next;
end

fclose(MatlabSimfile_Discrete);
