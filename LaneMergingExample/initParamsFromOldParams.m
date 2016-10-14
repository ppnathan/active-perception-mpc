function init_params = initParamsFromOldParams(oldParams, deltaT, num_q, num_cstate, ...
                                                   old_horizon, new_horizon)

x = zeros(num_cstate, old_horizon, num_q);
starting_idx = 0;
for i = 1:num_q
    x(:, :, i) = reshape(oldParams(starting_idx + 1 : starting_idx + num_cstate * old_horizon), ...
                   [num_cstate old_horizon]);
    starting_idx = starting_idx + num_cstate * old_horizon;
end

u = zeros(old_horizon, num_q);
for i = 1:num_q
    u(:, i) = oldParams(starting_idx + 1 : starting_idx + old_horizon);
    starting_idx = starting_idx + old_horizon;
end

safeDist = 7;

A = [1, deltaT, 0, 0;
     0, 1,      0, 0;
     0, 0,      1, deltaT;
     0, 0,      0, 1];
B = [0; 0; 0; deltaT];

g{1} = [0; 0;       0; 0];
g{2} = [0; deltaT;  0; 0];
g{3} = [0; -deltaT; 0; 0];

x_init = zeros(num_cstate, new_horizon, num_q);
u_init = ones(new_horizon, num_q);

if (new_horizon <= old_horizon-1)
    x_init = x(:, 2:(new_horizon+1), :);
    u_init = u(2:(new_horizon+1), :);
else
    x_init(:, 1:(old_horizon-1), :) = x(:, 2:old_horizon, :);
    u_init(1:(old_horizon-1), :) = u(2:old_horizon, :);
    for t = old_horizon:new_horizon
        for i = 1:num_q
            if (max(x_init(1, t-1, i), x_init(3, t-1, i)) > 0 && ...
                    (x_init(1, t-1, i) - x_init(3, t-1, i)) < safeDist)
                x_init(:, t, i) = x_init(:, t-1, i);
            else
                x_init(:, t, i) = A * x_init(:, t-1, i) + B * u_init(t, i) + g{i};
            end
        end
    end
end
    
init_params = [x_init(:); u_init(:)];


end