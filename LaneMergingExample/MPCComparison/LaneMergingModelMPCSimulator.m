close all;
clear all;

num_qq = 12;
num_zq = 1;
num_cstate = 4;
gamma = 0.999;

SimTime = 200;

x_sim = zeros(SimTime+1, num_cstate);
q_sim = ones(SimTime+1, 1);
%zq_sim = ones(SimTime+1, 1);
%b_q = zeros(1, num_qq); b_q(1) = 1;
%b_q = repmat(b_q, SimTime+1, 1);
control = zeros(SimTime+1, 1);

MatlabSimfile_Discrete = ...
    fopen('./OutFiles/MatlabSimfile_Discrete.txt', 'w');

x_sim(1, :) = [-270 30 -300 32];
q_sim(1, 1) = 3;

horizon = 60;

deltaT = 0.2;

% x1 = sdpvar(horizon+1, 1);
% v1 = sdpvar(horizon+1, 1);
% x2 = sdpvar(horizon+1, 1);
% v2 = sdpvar(horizon+1, 1);
% %u = sdpvar(repmat(1,1,horizon),repmat(1,1,horizon));
% u = sdpvar(horizon, 1);
% s0 = sdpvar(4, 1);
% 
% 
% objective = -sum(x2) + norm(u, 1);
% 
% constraints = [];
% constraints = [constraints, x1(1) == s0(1), v1(1) == s0(2), x2(1) == s0(3), v2(1) == s0(4)];
% %     s0 == [-450; 30; -500; -30]];
% for i = 1:max(5, contraints_horizon)
%     constraints = [constraints, x1(i+1) == x1(i) + deltaT*v1(i), v1(i+1) == v1(i), ...
%         x2(i+1) == x2(i) + deltaT*v2(i), v2(i+1) == v2(i) + deltaT * u(i)];
% %         ((x1(i+1) + deltaT*i*(i+1)/2)*cos(pi/4) - x2(i+1)*sin(pi/4)) >= (7*sin(pi/4)) | ...
% %         ((x1(i+1) + deltaT*i*(i+1)/2)*cos(pi/4) - x2(i+1)*sin(pi/4)) <= -(7*sin(pi/4)) | ...
% %         (((x1(i+1) + deltaT*i*(i+1)/2) + x2(i+1) + 6) <= 0), ...
% %         ((x1(i+1) - deltaT*i*(i+1)/2)*cos(pi/4) - x2(i+1)*sin(pi/4)) >= (7*sin(pi/4)) | ...
% %         ((x1(i+1) - deltaT*i*(i+1)/2)*cos(pi/4) - x2(i+1)*sin(pi/4)) <= -(7*sin(pi/4)) | ...
% %         (((x1(i+1) - deltaT*i*(i+1)/2) + x2(i+1) + 6) <= 0)];
% end
% constraints = [constraints, -1 <= u <= 1, ...
%    (x1*cos(pi/4) - x2*sin(pi/4)) >= (7*sin(pi/4)) | (x1*cos(pi/4) - x2*sin(pi/4)) <= -(7*sin(pi/4)) | ((x1 + x2 + 6) <= 0)];
%     % Set some options for YALMIP and solver
% options = sdpsettings('verbose',1);
% % Solve the problem
% MPCController = optimizer(constraints, objective, options, [s0; contraints_horizon], u);

% Analyze error flags
% if sol.problem == 0
% % Extract and display value
%     solution = value(u(1))
% else
%     display('Hmm, something went wrong!');
%     sol.info
%     yalmiperror(sol.problem)
% end

for k = 1:SimTime
    k
    
    x1 = sdpvar(horizon+1, 1);
    v1 = sdpvar(horizon+1, 1);
    x2 = sdpvar(horizon+1, 1);
    v2 = sdpvar(horizon+1, 1);
    %u = sdpvar(repmat(1,1,horizon),repmat(1,1,horizon));
    u = sdpvar(horizon, 1);
    s0 = sdpvar(4, 1);


    objective = -sum(x2);
    contraints_horizon = max(5, min(horizon, floor(abs(min(x_sim(k, 1), x_sim(k, 3)))/30/deltaT)));
    constraints = [];
    constraints = [constraints, x1(1) == s0(1), v1(1) == s0(2), x2(1) == s0(3), v2(1) == s0(4)];

    for i = 1:contraints_horizon
        constraints = [constraints, x1(i+1) == x1(i) + deltaT*v1(i), v1(i+1) == v1(i), ...
            x2(i+1) == x2(i) + deltaT*v2(i), v2(i+1) == v2(i) + deltaT * u(i)];
    %         ((x1(i+1) + deltaT*i*(i+1)/2)*cos(pi/4) - x2(i+1)*sin(pi/4)) >= (7*sin(pi/4)) | ...
    %         ((x1(i+1) + deltaT*i*(i+1)/2)*cos(pi/4) - x2(i+1)*sin(pi/4)) <= -(7*sin(pi/4)) | ...
    %         (((x1(i+1) + deltaT*i*(i+1)/2) + x2(i+1) + 6) <= 0), ...
    %         ((x1(i+1) - deltaT*i*(i+1)/2)*cos(pi/4) - x2(i+1)*sin(pi/4)) >= (7*sin(pi/4)) | ...
    %         ((x1(i+1) - deltaT*i*(i+1)/2)*cos(pi/4) - x2(i+1)*sin(pi/4)) <= -(7*sin(pi/4)) | ...
    %         (((x1(i+1) - deltaT*i*(i+1)/2) + x2(i+1) + 6) <= 0)];
    end
    constraints = [constraints, -1 <= u <= 1, ...
                                (x1 - x2) >= 7 | (x1 - x2) <= -7 | ((x1 + x2 + 6) <= 0)];
    % Set some options for YALMIP and solver
    options = sdpsettings('verbose',1);
    % Solve the problem
    MPCController = optimizer(constraints, objective, options, s0, u);


    sigma = MPCController{[x_sim(k, :)']}
    control(k) = value(sigma(1));
    str = sprintf('q = %d, x = (%f, %f, %f, %f), sigma = %f', ...
                  q_sim(k), x_sim(k, 1), x_sim(k, 2), x_sim(k, 3), x_sim(k, 4), control(k));
    disp(str);
    fprintf(MatlabSimfile_Discrete, ...
        '%d  %f  %f  %f  %f  %f\n', ...
        q_sim(k), x_sim(k, 1), x_sim(k, 2), x_sim(k, 3), x_sim(k, 4), control(k));
    

    if (IsSimEnded(q_sim(k), x_sim(k, :)))
        break;
    end
    
    [q_next, x_next] = SampleNext(q_sim(k), x_sim(k, :), control(k), deltaT, true);
    
    x_sim(k+1, :) = x_next;
    q_sim(k+1) = q_next;
end

fclose(MatlabSimfile_Discrete);







