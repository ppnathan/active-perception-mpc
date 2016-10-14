function cost = costFn_Ipopt(vars, state_cov_inv, b0, alpha, beta, lambda, eta, ...
                                        num_cstate, num_q, horizon)

% b0 is column vector

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


cost_task_tmp = sum(x(3, :, :), 2); % result: [sum(J(x)); sum(J(y))]
cost_task = -sum(cost_task_tmp(:) .* b0); % result: b0(1)*sum(J(x)) + b0(2)*sum(J(y))
                                          % Use minus sign because it is a minimization problem
cost_u = sum(u.^2, 1);
cost_u = sum(cost_u(:) .* b0);

% calculate the entropy of our belief

b = b0(b0 > 0);
entropy = -sum(b.*log(b));

cost_KL = 0;

for i = 2:num_q
    tmp = num2cell(state_cov_inv(:, :, :, i), [1,2]);
    blkmat_tmp = blkdiag(tmp{:}); 
    for j = 1:(i-1)
        differences = x(:, :, i) - x(:, :, j); 
        cost_KL = cost_KL + differences(:)' * blkmat_tmp * differences(:);
    end
end

gamma = 1000;
cost_similarity = 0;
for i = 2:num_q
    cost_similarity = cost_similarity + sum((u(:, i) - u(:, 1)).^2);
end


D_eff = 11;
points = [-7.64, -5.73, -3.82, -1.91, 0, 1.91, 3.82, 5.73, 7.64;
          -7.64, -5.73, -3.82, -1.91, 0, 1.91, 3.82, 5.73, 7.64];
num_points = size(points, 2);
dist_to_points = zeros(num_points, horizon, num_q);
cost_constraints_tmp = zeros(1, horizon, num_q);
for i = 1:num_points
    dist_to_points(i, :, :) = sqrt((x(1, :, :) - points(1, i)).^2 + (x(3, :, :) - points(2, i)).^2);
    cost_constraints_tmp = cost_constraints_tmp + ...
                           0.5 * eta * (dist_to_points(i, :, :) < D_eff) .* ...
                                        (1./dist_to_points(i, :, :) - 1/D_eff).^2;
end
% cost_constraints_tmp = min (cost_constraints_tmp, 1e5);
cost_constraints_tmp = sum(cost_constraints_tmp, 2);
cost_constraints = sum(cost_constraints_tmp(:) .* b0);

cost = alpha * cost_task + 0.5 * beta * cost_u - entropy * lambda * 0.5 * cost_KL + ...
    0.5 * entropy * lambda * gamma * cost_similarity + cost_constraints;

% str = sprintf('task = %.3f, KL = %.3f, weighted KL = %.3f, constraints = %.3f, u = %.3f', ...
%     cost_task, cost_KL, entropy * lambda * 0.5 * cost_KL, cost_constraints, cost_u);
% disp(str);


end