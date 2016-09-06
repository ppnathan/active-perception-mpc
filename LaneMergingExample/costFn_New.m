function [cost, cost_grad] = costFn_New(vars, state_cov_inv, b0, lambda, epsilon, ...
                                        num_cstate, num_q, horizon)

% b0 is column vector

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

% u = vars(starting_idx + 1 : starting_idx + horizon); % Change notations
% starting_idx = starting_idx + horizon;
% 
% u_tilde = zeros(horizon, num_q);
% for i = 1:num_q
%     u_tilde(:, i) = vars(starting_idx + 1 : starting_idx + horizon);
%     starting_idx = starting_idx + horizon;
% end

cost_task_tmp = sum(x(3, :, :), 2); % result: [sum(J(x)); sum(J(y))]
cost_task = -sum(cost_task_tmp(:) .* b0); % result: b0(1)*sum(J(x)) + b0(2)*sum(J(y))
                                          % Use minus sign because it is a minimization problem

% calculate the entropy of our belief

b = b0(b0 > 0);
entropy = -sum(b.*log2(b));

cost_KL = 0;

for i = 2:num_q
    tmp = num2cell(state_cov_inv(:, :, :, i), [1,2]);
    blkmat_tmp = blkdiag(tmp{:}); 
    for j = 1:(i-1)
        differences = x_tilde(:, :, i) - x_tilde(:, :, j); % .* repmat(sqrt(horizon - (1:horizon) + 1), num_cstate, 1);
        cost_KL = cost_KL + differences(:)' * blkmat_tmp * differences(:);
    end
end

eta = 1e6;
D_eff = 7;
points = [-7  -3.5  0;
          -7  -3.5  0];
num_points = size(points, 2);
dist_to_points = zeros(num_points, horizon, num_q);
cost_constraints_tmp = zeros(1, horizon, num_q);
for i = 1:num_points
    dist_to_points(i, :, :) = sqrt((x(1, :, :) - points(1, i)).^2 + (x(3, :, :) - points(2, i)).^2);
    cost_constraints_tmp = cost_constraints_tmp + ...
                           0.5 * eta * (dist_to_points(i, :, :) < D_eff) .* (1./dist_to_points(i, :, :) - 1/D_eff).^2;
end
cost_constraints_tmp = sum(cost_constraints_tmp, 2);
cost_constraints = sum(cost_constraints_tmp(:) .* b0);

beta1 = 25;
beta2 = 25;
cost_u = sum(u.^2, 1);
cost_u = sum(cost_u(:) .* b0);
cost_u_tilde = sum(u_tilde.^2);

% cost_task
% cost_KL
% cost = cost_task + lambda * 0.5 * cost_KL;
cost = cost_task - entropy * lambda * 0.5 * cost_KL + cost_constraints + ...
       0.5 * beta * (cost_u + entropy * cost_u_tilde);

str = sprintf('task = %.3f, KL = %.3f, weighted KL = %.3f, constraints = %.3f, u = %.3f, u_tilde = %.3f', ...
    cost_task, cost_KL, entropy * lambda * 0.5 * cost_KL, cost_constraints, ...
    cost_u, entropy * cost_u_tilde);
disp(str);

if nargout > 1 % gradient required
    b0_rep = repmat(b0, 1, num_cstate * horizon)';
    grad_cost_task = [repmat([0; 0; -1; 0], num_q * horizon, 1) .* b0_rep(:); ...
                      zeros(num_cstate * horizon * num_q, 1); ...
                      zeros(num_q * horizon + horizon, 1)];
%     grad_cost_KL  = 0;
    grad_x_tilde = zeros(num_cstate * horizon, num_q);
    
    for i = 2:num_q
        tmp = num2cell(state_cov_inv(:, :, :, i), [1,2]);
        blkmat_tmp = blkdiag(tmp{:}); 
        for j = 1:i-1
            differences = x_tilde(:, :, i) - x_tilde(:, :, j);
            grad_tmp = blkmat_tmp * differences(:);
            grad_x_tilde(:, i) = grad_x_tilde(:, i) + grad_tmp;
            grad_x_tilde(:, j) = grad_x_tilde(:, j) - grad_tmp;
        end
    end
    
    grad_cost_KL = [zeros(num_cstate * horizon * num_q, 1);
                    grad_x_tilde(:);
                    zeros(num_q * horizon + horizon, 1)];


    grad_points = zeros(num_points, horizon, num_q);
    grad_cost_constraints_tmp = zeros(num_cstate, horizon, num_q);
    for i = 1:num_points
        grad_points(i, :, :) = eta * (dist_to_points(i, :, :) < D_eff) .* ...
                               (1/D_eff - 1./dist_to_points(i, :, :)) .* ...
                               1./(dist_to_points(i, :, :).^3);
        grad_cost_constraints_tmp(1, :, :) = grad_cost_constraints_tmp(1, :, :) + ...
                                             grad_points(i, :, :) .* (x(1, :, :) - points(1, i));
        grad_cost_constraints_tmp(3, :, :) = grad_cost_constraints_tmp(3, :, :) + ...
                                             grad_points(i, :, :) .* (x(3, :, :) - points(2, i));
    end

    grad_cost_constraints = [grad_cost_constraints_tmp(:).* b0_rep(:); ...
                             zeros(num_cstate * horizon * num_q, 1); ...
                             zeros(num_q * horizon + horizon, 1)];
    b0_rep2 = repmat(b0, 1, horizon)';                   
    grad_u = [zeros(num_cstate * horizon * num_q, 1);...
              zeros(num_cstate * horizon * num_q, 1); ...
              beta * u(:) .* b0_rep2(:);...
              beta * entropy * u_tilde(:)];
    cost_grad = grad_cost_task - entropy * lambda * grad_cost_KL + grad_cost_constraints + grad_u;
end

end