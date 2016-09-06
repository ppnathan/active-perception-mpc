function Jacobian = calJacobian(x_k, mode, deltaT)

num_cstate = length(x_k);
safeDist = 7;
mReactionDist = 25;

if (max(x_k(1), x_k(3)) > 0 && abs(x_k(1) - x_k(3)) < safeDist)
    Jacobian = eye(num_cstate);
else
    
    A = [1, deltaT, 0, 0;
         0, 1,      0, 0;
         0, 0,      1, deltaT;
         0, 0,      0, 1];
     
    Jacobian = A;

%     if (abs(x_k(1) - x_k(3)) > mReactionDist)
%         Jacobian = A;
%     else
%         Jacobian = A;
%     end

end

end