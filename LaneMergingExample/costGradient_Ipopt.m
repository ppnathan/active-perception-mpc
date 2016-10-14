function cost_grad = costGradient_Ipopt(vars, state_cov_inv, b0, alpha, beta, lambda, eta, ...
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


b = b0(b0 > 0);
entropy = -sum(b.*log(b));

D_eff = 11;
points = [-7.64, -5.73, -3.82, -1.91, 0, 1.91, 3.82, 5.73, 7.64;
          -7.64, -5.73, -3.82, -1.91, 0, 1.91, 3.82, 5.73, 7.64];
num_points = size(points, 2);

% str = sprintf('task = %.3f, KL = %.3f, weighted KL = %.3f, constraints = %.3f, u = %.3f', ...
%     cost_task, cost_KL, entropy * lambda * 0.5 * cost_KL, cost_constraints, cost_u);
% disp(str);


b0_rep = repmat(b0, 1, num_cstate * horizon)';
grad_cost_task = [repmat([0; 0; -1; 0], num_q * horizon, 1) .* b0_rep(:); ...
                  zeros(num_q * horizon, 1)];

%     grad_cost_KL  = 0;
grad_x = zeros(num_cstate * horizon, num_q);

b0_rep2 = repmat(b0, 1, horizon)';                   
grad_u = [ zeros(num_cstate * horizon * num_q, 1); ...
           beta * u(:) .* b0_rep2(:);];

for i = 2:num_q
    tmp = num2cell(state_cov_inv(:, :, :, i), [1,2]);
    blkmat_tmp = blkdiag(tmp{:}); 
    for j = 1:i-1
        differences = x(:, :, i) - x(:, :, j);
        grad_tmp = blkmat_tmp * differences(:);
        grad_x(:, i) = grad_x(:, i) + grad_tmp;
        grad_x(:, j) = grad_x(:, j) - grad_tmp;
    end
end

grad_cost_KL = [grad_x(:);
                zeros(num_q * horizon, 1)];
            
gamma = 1000;
grad_similarity = [ zeros(num_cstate * horizon * num_q, 1); ...
                    u(:, 1) - u(:, 2) + u(:, 1) - u(:, 3); ...
                    u(:, 2) - u(:, 1); ...
                    u(:, 3) - u(:, 1)];

dist_to_points = zeros(num_points, horizon, num_q);
cost_constraints_tmp = zeros(1, horizon, num_q);
for i = 1:num_points
    dist_to_points(i, :, :) = sqrt((x(1, :, :) - points(1, i)).^2 + (x(3, :, :) - points(2, i)).^2);
    cost_constraints_tmp = cost_constraints_tmp + ...
                           0.5 * eta * (dist_to_points(i, :, :) < D_eff) .* (1./dist_to_points(i, :, :) - 1/D_eff).^2;
end

% is_max_cost_constraints = cost_constraints_tmp >= 1e5;

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
% grad_cost_constraints_tmp(1, :, :) = grad_cost_constraints_tmp(1, :, :) .* (~is_max_cost_constraints);
% grad_cost_constraints_tmp(3, :, :) = grad_cost_constraints_tmp(3, :, :) .* (~is_max_cost_constraints);

grad_cost_constraints = [grad_cost_constraints_tmp(:).* b0_rep(:); ...
                         zeros(num_q * horizon, 1)];



cost_grad = alpha * grad_cost_task + beta * grad_u - entropy * lambda * grad_cost_KL + ...
            entropy * lambda * gamma * grad_similarity + grad_cost_constraints;
end