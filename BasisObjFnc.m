classdef BasisObjFnc < handle
    % Basis objective function for the PSM model, e.g.  
    properties
        J % j(x) is a function from q'_1 to varible index
        v % v(x) is the target function on q'_1
        m % vector of functions on q'_1
        H % quadraic cost
        f % linear cost
        c % constant term
        dH 
        df 
        dc 
        Aeq % equality constraint matrix 
        Ain % inequality constraint matrix 
        beq % equality constraint
        bin % inequality constraint
        qPrime1i % initial postion in qPrime1
        qPrime1f % final postion in qPrime1
        type % bookkeeping varible
        Wind % bookkeeping varible
    end
    
    methods(Static)
        function ind = getJ(cp, qPrime1)
            ind = zeros(length(qPrime1),1);
            for k = 1:length(qPrime1)
                if qPrime1(k) < cp(1)
                    ind(k) = 1;
                elseif qPrime1(k) > cp(end)
                    ind(k) = length(cp);
                else
                    for i = 1:length(cp)-1
                        if qPrime1(k) >= cp(i) && qPrime1(k) < cp(i+1)
                            ind(k) = i;
                        end
                    end
                end
            end
        end
    end
    
    methods
        function obj = BasisObjFnc(cp, m, v)
            obj.J = @(qPrime1) BasisObjFnc.getJ(cp, qPrime1);
            obj.v = v;
            obj.m = m;
        end
        
        function addConstraint2(obj, m1, m2, qPrime1)
            r = zeros(1,length(m1));
            for i = 1:length(m1)
                if ~isnumeric(m1{i})
                    f1 = m1{i};
                    r(i) = f1(qPrime1-0.00000001);
                end
                if ~isnumeric(m2{i})
                    f2 = m2{i};
                    r(i) = -f2(qPrime1);
                end
            end
            obj.Aeq = [obj.Aeq;r];
            obj.beq = [obj.beq;0];
        end
        
        function addConstraintE(obj, m, g, qPrime1)
            r = zeros(1,length(m));
            for i = 1:length(m)
                if ~isnumeric(m{i})
                    f1 = m{i};
                    r(i) = f1(qPrime1);
                end
            end
            b = g(qPrime1);
            obj.Aeq = [obj.Aeq;r];
            obj.beq = [obj.beq;b];
        end
        
        function addConstraintIn(obj, m, g, qPrime1)
            r = zeros(1,length(m));
            for i = 1:length(m)
                if ~isnumeric(m{i})
                    f1 = m{i};
                    r(i) = f1(qPrime1);
                end
            end
            b = g(qPrime1);
            obj.Ain = [obj.Ain;-r];
            obj.bin = [obj.bin;b];
        end
        
        function addConstraintIn2(obj, m, g, qPrime1)
            r = zeros(1,length(m));
            for i = 1:length(m)
                if ~isnumeric(m{i})
                    f1 = m{i};
                    r(i) = f1(qPrime1);
                end
            end
            b = g(qPrime1);
            obj.Ain = [obj.Ain;r];
            obj.bin = [obj.bin;b];
        end
        
        function build(obj)
            obj.H = zeros(length(obj.m));
            obj.f = zeros(length(obj.m),1);
            qPrime1 = linspace(obj.qPrime1i,obj.qPrime1f, 2000)';
            qPrime1Plot = linspace(obj.qPrime1i,obj.qPrime1f, 100)';
            f2 = obj.v;
            f2Val = f2(qPrime1);
            obj.c = trapz(qPrime1, f2Val.^2);
            obj.dH = zeros(length(obj.m),length(obj.m),length(qPrime1Plot));
            obj.df = zeros(length(obj.m),1,length(qPrime1Plot));
            obj.dc = f2(qPrime1Plot).^2;
            for i = 1:length(obj.m)
                if isnumeric(obj.m{i}) && obj.m{i} == 0 % this is the zero function   p'*m'*m*p-2*m*v*p+v^2
                    continue
                end
                f1 = obj.m{i};
                f1Val = f1(qPrime1);
                obj.f(i) = trapz(qPrime1,-2*f1Val.*f2Val);
                obj.df(i,1,:) = -2*f1(qPrime1Plot).*f2(qPrime1Plot);
                %                 for j = i:length(obj.m)
                for j = i:length(obj.m)
                    if isnumeric(obj.m{j}) && obj.m{j} == 0 % this is the zero function
                        continue
                    end
                    f3 = obj.m{j};
                    obj.H(i,j) = trapz(qPrime1,f1Val.*f3(qPrime1));
                    obj.H(j,i) = obj.H(i,j);
                    obj.dH(i,j,:) = f1(qPrime1Plot).*f3(qPrime1Plot);
                    obj.dH(j,i,:) = obj.dH(i,j,:);
                end
            end
        end
        
    end
end

