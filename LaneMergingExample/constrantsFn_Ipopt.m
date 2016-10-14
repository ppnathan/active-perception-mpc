function ceq = constrantsFn_Ipopt(vars, x0, epsilon, deltaT, num_cstate, num_q, ...
                                     horizon, constraint_horizon)

x = zeros(num_cstate, horizon, num_q);
starting_idx = 0;
for i = 1:num_q
    x(:, :, i) = reshape(vars(starting_idx + 1 : starting_idx + num_cstate * horizon), ...
                   [num_cstate horizon]);
    starting_idx = starting_idx + num_cstate * horizon;
end

u = zeros(horizon, num_q);
for i = 1:num_q
    u(:, i) = vars(starting_idx + 1 : starting_idx + horizon);
    starting_idx = starting_idx + horizon;
end


A = [1, deltaT, 0, 0;
     0, 1,      0, 0;
     0, 0,      1, deltaT;
     0, 0,      0, 1];
B = [0; 0;      0; deltaT];

g{1} = [0; 0;       0; 0];
g{2} = [0; deltaT;  0; 0];
g{3} = [0; -deltaT; 0; 0];

mReactionDist = 25;

num_ceq = num_q-1 + horizon * num_q * num_cstate;
count_ceq = 1;
ceq = zeros(num_ceq, 1);
for i = 2 : num_q
    ceq(count_ceq) = u(1, 1) - u(1, i);
    count_ceq = count_ceq + 1;
end

for i = 1 : num_q
    for t = 1 : horizon    
        if t == 1
            ceq(count_ceq:count_ceq+num_cstate-1) = double(abs(x0(1) - x0(3)) <= mReactionDist) * ...
                             (x(:, t, i) - (A * x0 + B * u(t, i) + g{i})) + ...
                             double(abs(x0(1) - x0(3)) > mReactionDist) * ...
                             (x(:, t, i) - (A * x0 + B * u(t, i)));
            count_ceq = count_ceq + num_cstate;

        else
            ceq(count_ceq:count_ceq+num_cstate-1) = double(abs(x(1, t-1, i) - x(3, t-1, i)) <= mReactionDist) * ...
                        (x(:, t, i) - (A * x(:, t-1, i) + B * u(t, i) + g{i})) + ...
                        double(abs(x(1, t-1, i) - x(3, t-1, i)) > mReactionDist) * ...
                        (x(:, t, i) - (A * x(:, t-1, i) + B * u(t, i)));
            count_ceq = count_ceq + num_cstate;
            
        end
    end
end


end