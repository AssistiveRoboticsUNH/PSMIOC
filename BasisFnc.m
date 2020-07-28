function m = BasisFnc(cp, p, functionsM, cut)
m1 = cell(length(p),1);
m = cell(length(p),1);
% m1
J = @(qPrime1) BasisObjFnc.getJ(cp, qPrime1);
for i = 1:length(p)
        m1{i} = 0;
        m{i} = 0;
   if p(i) ~= 0
       m1{i} = functionsM{p(i)};
       m{i} = @(qPrime1) m1{i}(qPrime1).*(J(qPrime1) == cut);
   end
end
end