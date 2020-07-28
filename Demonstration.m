classdef Demonstration < handle
    % container for demonstration
    
    properties
        t,D,Dd,Ddd,DPrime,DdPrime,DddPrime,dt,V, dPrimeOffset, dPrimeScale
    end
    
    methods
        function obj = Demonstration(t,D,Dd,Ddd,DPrime,DdPrime,DddPrime,dt,V)
            obj.t = t;
            obj.D = D;
            obj.Dd = Dd;
            obj.Ddd = Ddd;
            obj.DPrime = DPrime;
            obj.DdPrime = DdPrime;
            obj.DddPrime = DddPrime;
            obj.dt = dt;
            obj.V = V;
        end
        
        function obj = crop(obj,ind1,ind2)
            obj.t = obj.t(ind1:ind2);%-obj.t(ind1);
            obj.D = obj.D(:,ind1:ind2);
            obj.Dd = obj.Dd(:,ind1:ind2);
            obj.Ddd = obj.Ddd(:,ind1:ind2);
            obj.DPrime = obj.DPrime(:,ind1:ind2);
            obj.DdPrime = obj.DdPrime(:,ind1:ind2);
            obj.DddPrime = obj.DddPrime(:,ind1:ind2);
        end
        
    end
end

