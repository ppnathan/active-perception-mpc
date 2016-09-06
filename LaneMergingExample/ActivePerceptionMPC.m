close all; 
clear all;


num_q = 3;
num_cstate = 4;
gamma = 0.999;

SimTime = 200;

x_sim = zeros(num_cstate, SimTime+1);
belief_sim = zeros(num_q, SimTime+1);
control = zeros(SimTime+1, 1);

MatlabSimfile_Discrete = ...
    fopen('OutFiles/MatlabSimfile_Discrete.txt', 'w');

x_sim(:, 1) = [-274; 30; -300; 31];
belief_sim(:, 1) = [1/3; 1/3; 1/3];
q_sim = 3;

horizon = 60;

deltaT = 0.2;
mReactionDist = 25;

lambda1 = 100;
lambda2 = 1;


for k = 1:SimTime
    k
    %horizon = max(5, min(horizon, floor(abs(min(x_sim(1, k), x_sim(3, k)))/30/deltaT)))
%     log(belief_sim(:, k))
    fObj = @(vars)costFn(vars, belief_sim(:, k), lambda1, lambda2, num_cstate, num_q, horizon);
    
%    x0 = [repmat(x_sim(:, k), horizon, 1); ones(horizon, 1); ones(num_cstate * num_q * horizon, 1)];
    x0 = initialize_params(q_sim, x_sim(:, k), mReactionDist, deltaT, horizon);
    
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    lb = [-Inf*ones(num_cstate*horizon, 1); -ones(horizon, 1); -Inf*ones(num_cstate * num_q * horizon, 1)];
    ub = [Inf*ones(num_cstate*horizon, 1); ones(horizon, 1); Inf*ones(num_cstate * num_q * horizon, 1)];
    nonlcon = @(vars)constrantsFn(vars, x_sim(:, k), deltaT, num_cstate, num_q, ...
                                  horizon, horizon);
    options = optimoptions('fmincon', ... %'Display','iter', ...
                           'Algorithm','sqp', 'GradObj','on', 'MaxFunEvals', 20000);
    
    [opt_arg, min_value] = fmincon(fObj, x0, A, b, Aeq, beq, lb, ub, nonlcon, options);
    min_value
    sigma = opt_arg(num_cstate*horizon + 1);
    opt_arg = [reshape(opt_arg(1:num_cstate * horizon), num_cstate, []); 
               opt_arg((num_cstate*horizon + 1):(num_cstate*horizon + horizon))';
               reshape(opt_arg((num_cstate*horizon + horizon+1):end), num_cstate * num_q, [])];
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