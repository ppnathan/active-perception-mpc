function Jaco = dummyJaco(x)
num_q = 3;
horizon = 5;
num_cstate = 4;
num_ceq = num_q-1 + horizon * num_q * num_cstate;
Jaco = sparse(zeros(num_ceq, num_cstate * horizon * num_q + num_q*horizon));
end