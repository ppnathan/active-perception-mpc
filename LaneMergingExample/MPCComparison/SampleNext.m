function [q_next, x_next] = SampleNext(q_k, x_k, sigma, deltaT, hasNoise)

q_next = q_k;

safeDist = 7;
reactionDist = 25;

if (max(x_k(1), x_k(3)) > 0 && abs(x_k(1) - x_k(3)) < safeDist)
    x_next(1) = x_k(1);
    x_next(2) = x_k(2);
    x_next(3) = x_k(3);
    x_next(4) = x_k(4);
else
    humanDriver = floor(q_next / 3);
    if (humanDriver == 0)
        humanInput = 0;
    elseif (humanDriver == 1)
        if (abs(x_k(1) - x_k(3)) < reactionDist)      
            humanInput = 1;
        else
            humanInput = 0;
        end
    elseif (humanDriver == 2)
        if (abs(x_k(1) - x_k(3)) < reactionDist)
            humanInput = -1;
        else
            humanInput = 0;
        end
    else %if (humanDriver == REASONABLE) 
        if (abs(x_k(1) - x_k(3)) < reactionDist)
            if (x_k(2) > x_k(4))
                humanInput = 1;
            else
                humanInput = -1;
            end
        else
            humanInput = 0;
        end
    end
    
    if (hasNoise)
        x_next(1) = x_k(1) + deltaT * x_k(2) + randn(1)*0.1;
        x_next(2) = x_k(2) + deltaT * humanInput + randn(1)*0.1;
        x_next(3) = x_k(3) + deltaT * x_k(4) + randn(1)*0.1;
        x_next(4) = x_k(4) + deltaT * sigma + randn(1)*0.1;
    else
        x_next(1) = x_k(1) + deltaT * x_k(2);
        x_next(2) = x_k(2) + deltaT * humanInput;
        x_next(3) = x_k(3) + deltaT * x_k(4);
        x_next(4) = x_k(4) + deltaT * sigma;
    end
    
end

end