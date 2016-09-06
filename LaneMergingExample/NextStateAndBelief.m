function [x_next, b_next] = NextStateAndBelief(q_k, x_k, b_k, u_k, deltaT)

safeDist = 7;
mReactionDist = 25;

num_q = 3;
num_cstate = 4;

if (max(x_k(1), x_k(3)) > 0 && abs(x_k(1) - x_k(3)) < safeDist)
    x_next(1) = x_k(1);
    x_next(2) = x_k(2);
    x_next(3) = x_k(3);
    x_next(4) = x_k(4);
    b_next = b_k;
else
    
    A = [1, deltaT, 0, 0;
         0, 1,      0, 0;
         0, 0,      1, deltaT;
         0, 0,      0, 1];
    B = [0; 0; 0; deltaT];

    g{1} = [0; 0;       0; 0];
    g{2} = [0; deltaT;  0; 0];
    g{3} = [0; -deltaT; 0; 0];

    mDistNoiseStd = 0.1;
    mVelNoiseStd = 0.1;

    noise = [randn(1) * mDistNoiseStd; randn(1) * mVelNoiseStd; ...
             randn(1) * mDistNoiseStd; randn(1) * mVelNoiseStd];

    if (abs(x_k(1) - x_k(3)) > mReactionDist)
        x_next = A * x_k + B * u_k + noise;
        b_next = b_k;
    else
        x_next = A * x_k + B * u_k + g{q_k} + noise;
        
        differences = zeros(num_cstate, num_q);
        sum_b = 0;
        b_temp = zeros(num_q, 1);
        for i = 1:num_q
            differences(:, i) = x_next - (A * x_k + B * u_k + g{i});
            b_temp(i) = b_k(i) * normpdf(differences(1, i), 0, 0.1) * ...
                                 normpdf(differences(2, i), 0, 0.1) * ...
                                 normpdf(differences(3, i), 0, 0.1) * ...
                                 normpdf(differences(4, i), 0, 0.1);
            sum_b = sum_b + b_temp(i);
        end

        if (sum_b ~= 0)
            b_next = b_temp/sum_b;
        else
            b_next = b_k;
        end
    end

end
end