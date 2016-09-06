function x_next = DynamicalModels(q_k, x_k, u, isReactive, deltaT)

A = [1, deltaT, 0, 0;
     0, 1,      0, 0;
     0, 0,      1, deltaT;
     0, 0,      0, 1];
B = [0; 0;      0; deltaT];

g{1} = [0; 0;       0; 0];
g{2} = [0; deltaT;  0; 0];
g{3} = [0; -deltaT; 0; 0];

mReactionDist = 25;

implies(~isReactive, x_next == A * x_k + B * u);
implies(isReactive, x_next == A * x_k + B * u + g{q_k});

% if (abs(x_k(1) - x_k(3)) >= mReactionDist)
%     x_next = A * x_k + B * u;
% else
%     x_next = A * x_k + B * u + g{q_k};
% end

end