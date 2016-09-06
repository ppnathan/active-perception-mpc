function [c, ceq] = constrantsFn(vars, x0, deltaT, num_cstate, num_q, horizon, constraint_horizon)

ceq = [];
c = [];

starting_idx = 1;
x = reshape(vars(starting_idx : num_cstate*horizon), [num_cstate horizon]);
starting_idx = num_cstate * horizon;
u = vars(starting_idx + 1 : starting_idx + horizon);
starting_idx = starting_idx + horizon;

xi = reshape(vars(starting_idx + 1 : starting_idx + num_cstate * num_q * horizon), ...
             [num_cstate num_q horizon]);

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
for t = 1: horizon
    for i = 1: num_q
        if t == 1
            
            ceq = [ceq; (abs(x0(1) - x0(3)) <= mReactionDist) * ...
                        (x(:, t) - (A * x0 + B * u(t) + g{i}) - SysSqrtCov .* xi(:, i, t)) + ...
                        (abs(x0(1) - x0(3)) > mReactionDist) * ...
                        (x(:, t) - (A * x0 + B * u(t)) - SysSqrtCov .* xi(:, i, t))]; 
%             if (abs(x0(1) - x0(3)) <= mReactionDist)
%                 ceq = [ceq; x(:, t) - (A * x0 + B * u(t) + g{i}) - SysSqrtCov .* xi(:, i, t)];
%             else
%                 ceq = [ceq; x(:, t) - (A * x0 + B * u(t)) - SysSqrtCov .* xi(:, i, t)];
%             end
        else
            ceq = [ceq; (abs(x(1, t-1) - x(3, t-1)) <= mReactionDist) * ...
                        (x(:, t) - (A * x(:, t-1) + B * u(t) + g{i}) - SysSqrtCov .* xi(:, i, t)) + ...
                        (abs(x(1, t-1) - x(3, t-1)) > mReactionDist) * ...
                        (x(:, t) - (A * x(:, t-1) + B * u(t)) - SysSqrtCov .* xi(:, i, t))];
%             if (abs(x(1, t-1) - x(3, t-1)) <= mReactionDist)
%                 ceq = [ceq; x(:, t) - (A * x(:, t-1) + B * u(t) + g{i}) - SysSqrtCov .* xi(:, i, t)];
%             else
%                 ceq = [ceq; x(:, t) - (A * x(:, t-1) + B * u(t)) - SysSqrtCov .* xi(:, i, t)];
%             end
        end
    end
end

c = [c;  ...
     - (((x(1, 1:constraint_horizon) * cos(pi/4) - x(3, 1:constraint_horizon) * sin(pi/4))' >= (7*sin(pi/4))) | ...
        ((x(1, 1:constraint_horizon) * cos(pi/4) - x(3, 1:constraint_horizon) * sin(pi/4))' <= -(7*sin(pi/4))) | ...
        (((x(1, 1:constraint_horizon) + x(3, 1:constraint_horizon) + 6)' <= 0))) + 0.5];


end