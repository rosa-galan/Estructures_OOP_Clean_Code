classdef ForceZComputer < handle

    properties (Access = private)
        rho
        V
        Cl
        lambda
        g
        b
        c
        Me
        liftForceDist 
        liftIntegral
        liftMomentDist
        weightForceDist
        weightMomentDist
        y
    end

    properties (Access = public)
         Sz
    end
    
    methods (Access = public)
        
        function obj = ForceZComputer(cParams)
            obj.init(cParams);  
        end

        function compute(obj)
            obj.computeForce();
        end

    end
    
    methods (Access = private)

        function init(obj,cParams)
            obj.rho = cParams.rho;
            obj.V = cParams.V;
            obj.Cl = cParams.Cl;
            obj.lambda = cParams.lambda;
            obj.g = cParams.g;
            obj.b = cParams.b;
            obj.c = cParams.c;
            obj.Me = cParams.Me;
            obj.liftForceDist = cParams.liftForceDist;
            obj.liftIntegral = cParams.liftIntegral;  
            obj.liftMomentDist = cParams.liftMomentDist;
            obj.weightMomentDist = cParams.weightMomentDist;
            obj.weightForceDist = cParams.weightForceDist;
        end

        function computeForce(obj)
            obj.Sz = double(obj.liftIntegral - obj.weightForceDist * obj.b/1000 - obj.Me * obj.g);
        end
        
    end
    
end
