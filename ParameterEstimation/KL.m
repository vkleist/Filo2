function D = KL(p1,p2)

% p1 and p2 are probabilities in some discrete state space (pdf)
%H = -sum([p1(1) p1(2:end)-p1(1:end-1)].*([p2(1) p2(2:end)-p2(1:end-1)]).*step);

p1(p1==0) = 1e-13;
p2(p2==0) = 1e-13;

D = sum(p1.*log(p1./p2));

