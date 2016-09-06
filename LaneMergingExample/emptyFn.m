function [x, g] = emptyFn(vars, num_cstate, num_q, horizon)
x = 0;
g = zeros(num_cstate * horizon * num_q * 2 + horizon + horizon * num_q, 1);
end