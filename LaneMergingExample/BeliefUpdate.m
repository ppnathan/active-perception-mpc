function b_next = BeliefUpdate(b_q, diffs)

sum_b = b_q(1) * normpdf(diffs(1)) + b_q(2) * normpdf(diffs(2)) + b_q(3) * normpdf(diffs(3));

% if (sum_b ~= 0)
    b_next(1, 1) = b_q(1) * normpdf(diffs(1)) / sum_b;
    b_next(2, 1) = b_q(2) * normpdf(diffs(2)) / sum_b;
    b_next(3, 1) = b_q(3) * normpdf(diffs(3)) / sum_b;
% else
%     b_next = b_q;
% end

end