function init_params = initialize_params_New(q_k, x_k, mReactionDist, deltaT, horizon)

num_cstate = length(x_k);
num_q = 3;
x_init = zeros(num_cstate, horizon, num_q);
% x_tilde_init = zeros(num_cstate, horizon, num_q);
u_init = ones(horizon, num_q);
% u_init(:, 2) = -ones(horizon, 1);
% u_tilde_init = ones(horizon, 1);

    
A = [1, deltaT, 0, 0;
     0, 1,      0, 0;
     0, 0,      1, deltaT;
     0, 0,      0, 1];
B = [0; 0; 0; deltaT];

g{1} = [0; 0;       0; 0];
g{2} = [0; deltaT;  0; 0];
g{3} = [0; -deltaT; 0; 0];

safeDist = 7;

for t = 1:horizon
    if t == 1
        for i = 1:num_q
            if (max(x_k(1), x_k(3)) > 0 && (x_k(1) - x_k(3)) < safeDist)
                x_init(:, t, i) = x_k;
%                 x_tilde_init(:, t, i) = x_k;
            else
                x_init(:, t, i) = A * x_k + B * u_init(t, i);
%                 x_tilde_init(:, t, i) = A * x_k + B * u_tilde_init(t);
            end
        end
    else
        for i = 1:num_q
            if (max(x_init(1, t-1, i), x_init(3, t-1, i)) > 0 && ...
                    (x_init(1, t-1, i) - x_init(3, t-1, i)) < safeDist)
                x_init(:, t, i) = x_init(:, t-1, i);
            else
                x_init(:, t, i) = A * x_init(:, t-1, i) + B * u_init(t, i);
            end
        end
        
%         for i = 1:num_q
%             if (max(x_tilde_init(1, t-1, i), x_tilde_init(3, t-1, i)) > 0 && ...
%                     (x_tilde_init(1, t-1, i) - x_tilde_init(3, t-1, i)) < safeDist)
%                 x_tilde_init(:, t, i) = x_tilde_init(:, t-1, i);
%             else
%                 x_tilde_init(:, t, i) = A * x_tilde_init(:, t-1, i) + B * u_tilde_init(t);
%             end
%         end
    end
end

init_params = [x_init(:)', u_init(:)'];


end