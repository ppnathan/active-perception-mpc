function ceq_Jaco = constrantsJacobian_Ipopt(vars, x0, epsilon, deltaT, num_cstate, num_q, ...
                                     horizon, constraint_horizon)


x = zeros(num_cstate, horizon, num_q);
starting_idx = 0;
for i = 1:num_q
    x(:, :, i) = reshape(vars(starting_idx + 1 : starting_idx + num_cstate * horizon), ...
                   [num_cstate horizon]);
    starting_idx = starting_idx + num_cstate * horizon;
end

% x_tilde = zeros(num_cstate, horizon, num_q);
% for i = 1:num_q
%     x_tilde(:, :, i) = reshape(vars(starting_idx + 1 : starting_idx + num_cstate * horizon), ...
%                          [num_cstate horizon]);
%     starting_idx = starting_idx + num_cstate * horizon;
% end

u = zeros(horizon, num_q);
for i = 1:num_q
    u(:, i) = vars(starting_idx + 1 : starting_idx + horizon);
    starting_idx = starting_idx + horizon;
end

% u_tilde = vars(starting_idx + 1 : starting_idx + horizon);


A = [1, deltaT, 0, 0;
     0, 1,      0, 0;
     0, 0,      1, deltaT;
     0, 0,      0, 1];
B = [0; 0;      0; deltaT];

%g{1} = [0; 0;       0; 0];
%g{2} = [0; deltaT;  0; 0];
%g{3} = [0; -deltaT; 0; 0];

% mReactionDist = 25;

% num_ceq = num_q-1 + horizon * num_q * num_cstate;
count_ceq = 1;
% ceq_Jaco = zeros(num_ceq, num_cstate * horizon * num_q + num_q * horizon);

u_eq_con_grad = zeros(num_q-1, num_cstate * horizon * num_q + num_q*horizon);
for i = 2 : num_q
%     u(1, 1) - u(1, i);
    u_eq_con_grad(count_ceq, :) = [ zeros(1, num_cstate * horizon * num_q), ...
                               1, zeros(1, horizon-1 + (i-2) * horizon), ...
                               1, zeros(1, horizon-1 + (num_q-i) * horizon)];

    count_ceq = count_ceq + 1;
end

A_blk = num2cell(repmat(A, 1, 1, horizon - 1), [1,2]);
x_gradient_block = sparse( eye(num_cstate * horizon) + ...
                           [ zeros(num_cstate, num_cstate * horizon);
                           [ -blkdiag(A_blk{:}), zeros(num_cstate * (horizon-1), num_cstate)]]);
B_blk = num2cell(repmat(B, 1, 1, horizon), [1, 2]);
u_gradient_block =  sparse(blkdiag(B_blk{:}));
x_grad_zeros_blk = sparse(num_cstate * horizon, num_cstate * horizon);
u_grad_zeros_blk = sparse(num_cstate * horizon, horizon);
ceq_Jaco = [ sparse(u_eq_con_grad);
             x_gradient_block, x_grad_zeros_blk, x_grad_zeros_blk, u_gradient_block, u_grad_zeros_blk, u_grad_zeros_blk;
             x_grad_zeros_blk, x_gradient_block, x_grad_zeros_blk, u_grad_zeros_blk, u_gradient_block, u_grad_zeros_blk;
             x_grad_zeros_blk, x_grad_zeros_blk, x_gradient_block, u_grad_zeros_blk, u_grad_zeros_blk, u_gradient_block];

% for i = 1 : num_q
%     for t = 1 : horizon
%         if t == 1
%             ceq(count_ceq:count_ceq+num_cstate-1) = double(abs(x0(1) - x0(3)) <= mReactionDist) * ...
%                              (x(:, t, i) - (A * x0 + B * u(t, i) + g{i})) + ...
%                              double(abs(x0(1) - x0(3)) > mReactionDist) * ...
%                              (x(:, t, i) - (A * x0 + B * u(t, i)));
%             ceq_Jaco(count_ceq:count_ceq+num_cstate-1, :) = ...
%                 zeros(num_cstate, );
%             
%             count_ceq = count_ceq + num_cstate;         
% 
%         else
%             ceq(count_ceq:count_ceq+num_cstate-1) = double(abs(x(1, t-1, i) - x(3, t-1, i)) <= mReactionDist) * ...
%                         (x(:, t, i) - (A * x(:, t-1, i) + B * u(t, i) + g{i})) + ...
%                         double(abs(x(1, t-1, i) - x(3, t-1, i)) > mReactionDist) * ...
%                         (x(:, t, i) - (A * x(:, t-1, i) + B * u(t, i)));
%             count_ceq = count_ceq + num_cstate;
%         end
%     end
% end


end