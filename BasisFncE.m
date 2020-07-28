function g = BasisFncE(cp, g, cut)
m1 = {};
m = {};
J = @(qPrime1) BasisObjFnc.getJ(cp, qPrime1);
g = @(qPrime1) g(qPrime1).*(J(qPrime1) == cut);
end