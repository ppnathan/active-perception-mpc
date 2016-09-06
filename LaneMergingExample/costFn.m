function [cost, cost_grad] = costFn(vars, b0, lambda1, lambda2, num_cstate, num_q, horizon)

% b0 is column vector

starting_idx = 1;
x = reshape(vars(starting_idx : num_cstate*horizon), [num_cstate horizon]);
starting_idx = num_cstate * horizon;
u = vars(starting_idx + 1 : starting_idx + horizon);
starting_idx = starting_idx + horizon;
%b = reshape(vars(starting_idx + 1 : starting_idx + num_q*horizon), [num_q horizon]);
%starting_idx = starting_idx + num_q*horizon;

xi = reshape(vars(starting_idx + 1 : starting_idx + num_cstate * num_q * horizon), ...
             [num_cstate num_q horizon]);

cost_task = 0; %-sum(x(3, :));
cost_disturbance = lambda2 * sum(sum(sum(repmat(b0', num_cstate, 1, horizon) .* xi.^2)));

b_last = sum(sum(log(normpdf(xi))), 3) + log(b0');
b_last = ones(1, num_q)./ b_last;
b_last = exp(b_last') / sum(exp(b_last)); % now b_last is a column vector.
cost_information = -lambda1 * sum(b_last .* log(b_last));

cost = cost_task + cost_disturbance + cost_information;

if nargout > 1 % gradient required
    grad_cost_task = 0; %[repmat([0; 0; -1; 0], horizon, 1); zeros(horizon + num_q*horizon , 1)];
    grad_temp = lambda2 * 2 * repmat(b0', num_cstate, 1, horizon) .* xi;
    grad_cost_disturbance = [zeros(num_cstate * horizon, 1);
                             zeros(horizon, 1);
                             grad_temp(:)];
    
    dH_db_N = log(b_last)+ ones(num_q, 1);
    temp = reshape(repmat(b_last*b_last' - diag(b_last), num_cstate, horizon), num_q, []) .* ...
           repmat(xi(:)', num_q, 1);
    grad_cost_information = [zeros(num_cstate * horizon, 1); zeros(horizon, 1); ...
                             - lambda1 * (dH_db_N' * temp)'];
    
    cost_grad = grad_cost_task + grad_cost_disturbance + grad_cost_information;
end

end