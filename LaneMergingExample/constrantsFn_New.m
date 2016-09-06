function [c, ceq] = constrantsFn_New(vars, x0, b0, epsilon, deltaT, num_cstate, num_q, ...
                                     horizon, constraint_horizon)


c = [];

x = zeros(num_cstate, horizon, num_q);
starting_idx = 0;
for i = 1:num_q
    x(:, :, i) = reshape(vars(starting_idx + 1 : starting_idx + num_cstate * horizon), ...
                   [num_cstate horizon]);
    starting_idx = starting_idx + num_cstate * horizon;
end

x_tilde = zeros(num_cstate, horizon, num_q);
for i = 1:num_q
    x_tilde(:, :, i) = reshape(vars(starting_idx + 1 : starting_idx + num_cstate * horizon), ...
                         [num_cstate horizon]);
    starting_idx = starting_idx + num_cstate * horizon;
end

u = zeros(horizon, num_q);
for i = 1:num_q
    u(:, i) = vars(starting_idx + 1 : starting_idx + horizon);
    starting_idx = starting_idx + horizon;
end

u_tilde = vars(starting_idx + 1 : starting_idx + horizon);


A = [1, deltaT, 0, 0;
     0, 1,      0, 0;
     0, 0,      1, deltaT;
     0, 0,      0, 1];
B = [0; 0;      0; deltaT];

g{1} = [0; 0;       0; 0];
g{2} = [0; deltaT;  0; 0];
g{3} = [0; -deltaT; 0; 0];

mReactionDist = 25;
SysSqrtCov = [0.1; 0.1; 0.1; 0.1];


num_ceq = num_q + horizon * num_q * num_cstate * 2 ;
count_ceq = 1;
ceq = zeros(num_ceq, 1);
% ceq = [];
for i = 1 : num_q
    ceq(count_ceq) = u(1, i) - u_tilde(1);
    count_ceq = count_ceq + 1;
%     ceq = [ceq; u(1, i) - u_tilde(1)];
end

for t = 1 : horizon
    for i = 1 : num_q
        if t == 1
            ceq(count_ceq:count_ceq+num_cstate-1) = double(abs(x0(1) - x0(3)) <= mReactionDist) * ...
                             (x(:, t, i) - (A * x0 + B * u(t, i) + g{i})) + ...
                             double(abs(x0(1) - x0(3)) > mReactionDist) * ...
                             (x(:, t, i) - (A * x0 + B * u(t, i)));
            count_ceq = count_ceq + num_cstate;         
           
            ceq(count_ceq:count_ceq+num_cstate-1) = double(abs(x0(1) - x0(3)) <= mReactionDist) * ...
                        (x_tilde(:, t, i) - (A * x0 + B * u_tilde(t) + g{i})) + ...
                        double(abs(x0(1) - x0(3)) > mReactionDist) * ...
                        (x_tilde(:, t, i) - (A * x0 + B * u_tilde(t))); 
            count_ceq = count_ceq + num_cstate;
%             ceq = [ceq; (abs(x0(1) - x0(3)) <= mReactionDist) * ...
%                         (x(:, t, i) - (A * x0 + B * u(t, i) + g{i})) + ...
%                         (abs(x0(1) - x0(3)) > mReactionDist) * ...
%                         (x(:, t, i) - (A * x0 + B * u(t, i)));
%                         
%                         (abs(x0(1) - x0(3)) <= mReactionDist) * ...
%                         (x_tilde(:, t, i) - (A * x0 + B * u_tilde(t) + g{i})) + ...
%                         (abs(x0(1) - x0(3)) > mReactionDist) * ...
%                         (x_tilde(:, t, i) - (A * x0 + B * u_tilde(t))); 
%                   ];
        else
            ceq(count_ceq:count_ceq+num_cstate-1) = double(abs(x(1, t-1, i) - x(3, t-1, i)) <= mReactionDist) * ...
                        (x(:, t, i) - (A * x(:, t-1, i) + B * u(t, i) + g{i})) + ...
                        double(abs(x(1, t-1, i) - x(3, t-1, i)) > mReactionDist) * ...
                        (x(:, t, i) - (A * x(:, t-1, i) + B * u(t, i)));
            count_ceq = count_ceq + num_cstate;
            
            ceq(count_ceq:count_ceq+num_cstate-1) = double(abs(x_tilde(1, t-1, i) - x_tilde(3, t-1, i)) <= mReactionDist) * ...
                        (x_tilde(:, t, i) - (A * x_tilde(:, t-1, i) + B * u_tilde(t) + g{i})) + ...
                        double(abs(x_tilde(1, t-1, i) - x_tilde(3, t-1, i)) > mReactionDist) * ...
                        (x_tilde(:, t, i) - (A * x_tilde(:, t-1, i) + B * u_tilde(t)));
            count_ceq = count_ceq + num_cstate;
%             ceq = [ceq; (abs(x(1, t-1, i) - x(3, t-1, i)) <= mReactionDist) * ...
%                         (x(:, t, i) - (A * x(:, t-1, i) + B * u(t, i) + g{i})) + ...
%                         (abs(x(1, t-1, i) - x(3, t-1, i)) > mReactionDist) * ...
%                         (x(:, t, i) - (A * x(:, t-1, i) + B * u(t, i)));
%                         
%                         (abs(x_tilde(1, t-1, i) - x_tilde(3, t-1, i)) <= mReactionDist) * ...
%                         (x_tilde(:, t, i) - (A * x_tilde(:, t-1, i) + B * u_tilde(t) + g{i})) + ...
%                         (abs(x_tilde(1, t-1, i) - x_tilde(3, t-1, i)) > mReactionDist) * ...
%                         (x_tilde(:, t, i) - (A * x_tilde(:, t-1, i) + B * u_tilde(t)));
%                   ];
        end
    end
end

% for i = 1:num_q
%     if b0(i) > epsilon
% %         c = [c;  ...
% %              - (((x_tilde(1, 1:constraint_horizon, i) * cos(pi/4) - x_tilde(3, 1:constraint_horizon, i) * sin(pi/4))' >= (7*sin(pi/4))) | ...
% %                 ((x_tilde(1, 1:constraint_horizon, i) * cos(pi/4) - x_tilde(3, 1:constraint_horizon, i) * sin(pi/4))' <= -(7*sin(pi/4))) | ...
% %                 ((x_tilde(1, 1:constraint_horizon, i) + x(3, 1:constraint_horizon, i) + 6)' <= 0)) + 0.5];
%         c = [c;  ...
%              - (((x_tilde(1, 1:constraint_horizon, i) - x_tilde(3, 1:constraint_horizon, i))' >= 7) | ...
%                 ((x_tilde(1, 1:constraint_horizon, i) - x_tilde(3, 1:constraint_horizon, i))' <= -7) | ...
%                 ((x_tilde(1, 1:constraint_horizon, i) + x(3, 1:constraint_horizon, i) + 6)' <= 0)) + 0.5];
% 
%     end
% end


end