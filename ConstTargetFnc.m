function v = ConstTargetFnc(cp, target, cut)
vIn = @(qPrime1) qPrime1*0+target;
J = @(qPrime1) BasisObjFnc.getJ(cp, qPrime1);
v = @(qPrime1) vIn(qPrime1).*( J(qPrime1) == cut);
end

