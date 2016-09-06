function init_params = initialize_params(q_k, x_k, mReactionDist, deltaT, horizon)

num_cstate = length(x_k);
num_q = 3;
x_init = zeros(num_cstate, horizon);
u_init = ones(horizon, 1);
    
A = [1, deltaT, 0, 0;
     0, 1,      0, 0;
     0, 0,      1, deltaT;
     0, 0,      0, 1];
B = [0; 0; 0; deltaT];

g{1} = [0; 0;       0; 0];
g{2} = [0; deltaT;  0; 0];
g{3} = [0; -deltaT; 0; 0];

safeDist = 7;

for i = 1:horizon
    if i == 1
        if (max(x_k(1), x_k(3)) > 0 && (x_k(1) - x_k(3)) < safeDist)
            x_init(:, i) = x_k;
        else
            x_init(:, i) = A * x_k + B * u_init(i);
%             if (abs(x_k(1) - x_k(3)) > mReactionDist)
%                 x_init(:, i) = A * x_k + B * u_init(i);
%             else
%                 x_init(:, i) = A * x_k + B * u_init(i) + g{q_k};
%             end
        end
    else
        if (max(x_init(1, i-1), x_init(3, i-1)) > 0 && (x_init(1, i-1) - x_init(3, i-1)) < safeDist)
            x_init(:, i) = x_init(:, i-1);
        else
            x_init(:, i) = A * x_init(:, i-1) + B * u_init(i);
%             if (abs(x_init(1, i-1) - x_init(3, i-1)) > mReactionDist)
%                 x_init(:, i) = A * x_init(:, i-1) + B * u_init(i);
%             else
%                 x_init(:, i) = A * x_init(:, i-1) + B * u_init(i) + g{q_k};
%             end
        end
    end
end

init_params = [x_init(:); u_init; ones(num_cstate * num_q * horizon, 1)];


end