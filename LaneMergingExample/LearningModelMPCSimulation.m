close all; 
clear all;


num_q = 3;
num_cstate = 4;
gamma = 0.999;

SimTime = 200;

x_sim = zeros(SimTime+1, num_cstate);
belief_sim = zeros(num_q, SimTime+1);
%q_sim = ones(SimTime+1, 1);
%zq_sim = ones(SimTime+1, 1);
%b_q = zeros(1, num_qq); b_q(1) = 1;
%b_q = repmat(b_q, SimTime+1, 1);
control = zeros(SimTime+1, 1);

MatlabSimfile_Discrete = ...
    fopen('OutFiles/MatlabSimfile_Discrete.txt', 'w');

x_sim(1, :) = [-300 30 -300 30];
belief_sim(:, 1) = [1/3; 1/3; 1/3];
%q_sim(1, 1) = 1;
q_sim = 1;

horizon = 60;

deltaT = 0.2;
mReactionDist = 25;

lambda1 = 1;
lambda2 = 2;

SysSqrtCov = [0.1; 0.1; 0.1; 0.1];

A = [1, deltaT, 0, 0;
     0, 1,      0, 0;
     0, 0,      1, deltaT;
     0, 0,      0, 1];
B = [0; 0;      0; deltaT];

g{1} = [0; 0;       0; 0];
g{2} = [0; deltaT;  0; 0];
g{3} = [0; -deltaT; 0; 0];

for k = 1:SimTime
    k
    
    x = sdpvar(num_cstate, horizon+1);
    u = sdpvar(1, horizon);
    b = sdpvar(num_q, horizon+1);
    s0 = sdpvar(num_cstate + num_q, 1);
    xi = sdpvar(num_q, horizon+1);
    react = binvar(1, horizon+1);


    objective = -sum(x(3, :)) - lambda1 * (b(1, horizon+1) * log(b(1, horizon+1)) + ...
                                      b(2, horizon+1) * log(b(2, horizon+1)) + ...
                                      b(3, horizon+1) * log(b(3, horizon+1))) ...
                         + lambda2 * sum(b(1, :) .* xi(1, :) .* xi(1, :) + ...
                                         b(2, :) .* xi(2, :) .* xi(2, :) + ...
                                         b(3, :) .* xi(3, :) .* xi(3, :));
   
    contraints_horizon = max(5, min(horizon, floor(abs(min(x_sim(k, 1), x_sim(k, 3)))/30/deltaT)));
    constraints = [];
    constraints = [constraints, x(:, 1) == s0(1:num_cstate), ...
                   b(:, 1) == s0(num_cstate+1:num_cstate + num_q)];

    %     s0 == [-450; 30; -500; -30]];
    for i = 1:contraints_horizon
        constraints = [constraints, iff(react(1, i), abs(x(1, i) - x(3, i)) <= mReactionDist)]
        for j = 1:num_q
            constraints = [constraints, ...
                implies(~react(1, i), x(:, i+1) == A * x(:, i) + B * u(1, i)), ...
                implies(react(1, i), x(:, i+1) == A * x(:, i) + B * u(1, i) + g{q_sim})];

%                 x(:, i+1) == DynamicalModels(j, x(:, i), u(1, i), react(i, :), deltaT) + SysSqrtCov * xi(j, i)];
        end
        constraints = [constraints, b(:, i+1) == BeliefUpdate(b(:, i), xi(:, i))];
    end
    
    constraints = [constraints, -1 <= u <= 1, b >0, ...
                   (x(1, :) * cos(pi/4) - x(3, :) * sin(pi/4)) >= (7*sin(pi/4)) | ...
                   (x(1, :) * cos(pi/4) - x(3, :) * sin(pi/4)) <= -(7*sin(pi/4)) | ...
                   ((x(1, :) + x(3, :) + 6) <= 0)];

    % Set some options for YALMIP and solver
    options = sdpsettings('verbose',1);
    % Solve the problem
    MPCController = optimizer(constraints, objective, options, s0, u);


    sigma = MPCController{[x_sim(k, :)'; belief_sim(:, k)]}
    control(k) = value(sigma(1));
    str = sprintf('q = %d, x = (%f, %f, %f, %f), belief = (%f, %f, %f), sigma = %f', ...
                  q_sim, x_sim(k, 1), x_sim(k, 2), x_sim(k, 3), x_sim(k, 4), ...
                  belief_sim(1, k), belief_sim(2, k), belief_sim(3, k), control(k));
    disp(str);
    fprintf(MatlabSimfile_Discrete, ...
        '%d  %f  %f  %f  %f  %f  %f  %f  %f\n', ...
        q_sim, x_sim(k, 1), x_sim(k, 2), x_sim(k, 3), x_sim(k, 4), ...
        belief_sim(1, k), belief_sim(2, k), belief_sim(3, k), control(k));
    

    if (IsSimEnded(q_sim(k), x_sim(k, :)))
        break;
    end
    
    [q_next, x_next] = SampleNext(q_sim(k), x_sim(k, :), control(k), deltaT, false);
    
    x_sim(k+1, :) = x_next;
end

fclose(MatlabSimfile_Discrete);