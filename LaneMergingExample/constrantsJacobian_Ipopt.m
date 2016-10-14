function ceq_Jaco = constrantsJacobian_Ipopt(deltaT, num_cstate, num_q, horizon)

A = [1, deltaT, 0, 0;
     0, 1,      0, 0;
     0, 0,      1, deltaT;
     0, 0,      0, 1];
B = [0; 0;      0; deltaT];

% mReactionDist = 25;

%  num_ceq = num_q-1 + horizon * num_q * num_cstate;
% ceq_Jaco = zeros(num_ceq, num_cstate * horizon * num_q + num_q * horizon);
count_ceq = 1;

u_eq_con_grad = zeros(num_q-1, num_cstate * horizon * num_q + num_q*horizon);
for i = 2 : num_q
    u_eq_con_grad(count_ceq, :) = [ zeros(1, num_cstate * horizon * num_q), ...
                               1, zeros(1, horizon-1 + (i-2) * horizon), ...
                               -1, zeros(1, horizon-1 + (num_q-i) * horizon)];

    count_ceq = count_ceq + 1;
end

A_blk = num2cell(repmat(A, 1, 1, horizon - 1), [1,2]);
x_gradient_block = sparse( eye(num_cstate * horizon) + ...
                           [ zeros(num_cstate, num_cstate * horizon);
                           [ -blkdiag(A_blk{:}), zeros(num_cstate * (horizon-1), num_cstate)]]);
B_blk = num2cell(repmat(-B, 1, 1, horizon), [1, 2]);
u_gradient_block =  sparse(blkdiag(B_blk{:}));
x_grad_zeros_blk = sparse(num_cstate * horizon, num_cstate * horizon);
u_grad_zeros_blk = sparse(num_cstate * horizon, horizon);
ceq_Jaco = [ sparse(u_eq_con_grad);
             x_gradient_block, x_grad_zeros_blk, x_grad_zeros_blk, u_gradient_block, u_grad_zeros_blk, u_grad_zeros_blk;
             x_grad_zeros_blk, x_gradient_block, x_grad_zeros_blk, u_grad_zeros_blk, u_gradient_block, u_grad_zeros_blk;
             x_grad_zeros_blk, x_grad_zeros_blk, x_gradient_block, u_grad_zeros_blk, u_grad_zeros_blk, u_gradient_block];

end