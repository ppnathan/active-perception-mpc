function state_cov_inv = initialize_cov_inv(param_init, deltaT, num_cstate, num_q, horizon)

x = zeros(num_cstate, horizon, num_q);

starting_idx = num_cstate * horizon * num_q;
for i = 1:num_q
    x(:, :, i) = reshape(param_init(starting_idx + 1 : starting_idx + num_cstate * horizon), ...
                   [num_cstate horizon]);
    starting_idx = starting_idx + num_cstate * horizon;
end

state_cov_inv(:, :, :, :) = ones([num_cstate num_cstate horizon num_q]);
state_cov = ones([num_cstate num_cstate horizon num_q]);

NoiseCov = diag([0.01; 0.01; 0.01; 0.01]);

for i = 1:num_q
    for t = 1:horizon
        if t ==1
            state_cov(:, :, t, i) = NoiseCov;
            state_cov_inv(:, :, t, i) = inv(state_cov(:, :, t, i));
        else
            Jacobian = calJacobian(x(:, t-1, i), i, deltaT);
            state_cov(:, :, t, i) = Jacobian * state_cov(:, :, t-1, i) * Jacobian' + NoiseCov;
            state_cov_inv(:, :, t, i) = inv(state_cov(:, :, t, i));
        end
    end
end

end