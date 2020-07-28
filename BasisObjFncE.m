classdef BasisObjFncE < handle
    
    properties
        g % g(x) is the expert attribute function on q'_1
        v % v(x) is the target function on q'_1
        c % expert basis fucntion value
        dc
        qPrime1i % starting point
        qPrime1f % ending point
        type
    end
    
    methods
        function obj = BasisObjFncE(g, v)
            obj.g = g;
            obj.v = v;
        end
        
        function build(obj)
            qPrime1 = linspace(obj.qPrime1i,obj.qPrime1f, 2000)';
            gVal = obj.g(qPrime1);
            vVal = obj.v(qPrime1);
            term1 =  gVal.^2;
            term2 = -2*gVal.*vVal;
            term3 = vVal.^2;
            obj.c = trapz(qPrime1,term1 +term2+term3); 
            qPrime1Plot = linspace(obj.qPrime1i,obj.qPrime1f, 100)';
            gVal = obj.g(qPrime1Plot);
            vVal = obj.v(qPrime1Plot);
            term1 =  gVal.^2;
            term2 = -2*gVal.*vVal;
            term3 = vVal.^2;
            obj.dc = term1 +term2+term3; 
        end
        
    end
end

