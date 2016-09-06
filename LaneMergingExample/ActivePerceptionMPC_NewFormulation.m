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

x_sim(:, 1) = [-270; 30; -300; 32];
belief_sim(:, 1) = [1/3; 1/3; 1/3];
q_sim = 3;

horizon = 30;

deltaT = 0.2;
mReactionDist = 25;

lambda = 50;
epsilon = 1e-2;


for k = 1:SimTime
    k
    horizon = max(5, min(horizon, floor(abs(min(x_sim(1, k), x_sim(3, k)))/30/deltaT)))
    
    params = initialize_params_New(q_sim, x_sim(:, k), mReactionDist, deltaT, horizon);
    for iter = 1:1
%         iter
        state_cov_inv = initialize_cov_inv(params, deltaT, num_cstate, num_q, horizon);

%         if iter == 1
%             disp('Try to find a feasible point.');
%             fObj = @(vars)emptyFn(vars, num_cstate, num_q, horizon);
%         else
%             disp('Do the original optimization problem.');
            fObj = @(vars)costFn_New(vars, state_cov_inv, belief_sim(:, k), lambda, epsilon, ...
                                 num_cstate, num_q, horizon);
%         end

        A = [];
        b = [];
        Aeq = [];
        beq = [];
        lb = [-Inf * ones(num_cstate * horizon * num_q * 2, 1); -ones(num_q * horizon + horizon, 1)];
        ub = [Inf * ones(num_cstate * horizon * num_q * 2, 1); ones(num_q * horizon + horizon, 1)];
        nonlcon = @(vars)constrantsFn_New(vars, x_sim(:, k), belief_sim(:, k), epsilon, ...
                                          deltaT, num_cstate, num_q, horizon, horizon);
        
        options = optimoptions('fmincon', 'Display','iter', ...
                               'Algorithm','sqp', 'GradObj','on');%, 'MaxFunEvals', 20000);

        disp('Solving optimization problem...')
        [params, min_value] = fmincon(fObj, params, A, b, Aeq, beq, lb, ub, nonlcon, options);
        min_value
        [x_params, x_tilde_params, u_params, u_tilde_params] = ...
                                                    parseParams(params, num_q, num_cstate, horizon);

%         sigma = u_params(1)
    end

    sigma = u_params(1);
%     params = [reshape(params(1:num_cstate * horizon), num_cstate, []); 
%                params((num_cstate*horizon + 1):(num_cstate*horizon + horizon))';
%                reshape(params((num_cstate*horizon + horizon+1):end), num_cstate * num_q, [])];

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
