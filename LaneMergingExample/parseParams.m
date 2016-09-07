function [x, u] = parseParams(vars, num_q, num_cstate, horizon)

%x = zeros(num_cstate, horizon, num_q);
starting_idx = 0;
for i = 1:num_q
    x{i} = reshape(vars(starting_idx + 1 : starting_idx + num_cstate * horizon), ...
                   [num_cstate horizon]);
    starting_idx = starting_idx + num_cstate * horizon;
end

% x_tilde = zeros(num_cstate, horizon, num_q);
% for i = 1:num_q
%     x_tilde{i} = reshape(vars(starting_idx + 1 : starting_idx + num_cstate * horizon), ...
%                          [num_cstate horizon]);
%     starting_idx = starting_idx + num_cstate * horizon;
% end

u = zeros(horizon, num_q);
for i = 1:num_q
    u(:, i) = vars(starting_idx + 1 : starting_idx + horizon);
    starting_idx = starting_idx + horizon;
end

% u_tilde = vars(starting_idx + 1 : starting_idx + horizon);

end