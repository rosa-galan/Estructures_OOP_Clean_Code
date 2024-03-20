classdef ForceZComputer < handle

    properties (Access = public)
         Sz
    end
    
    properties (Access = private)
        gravity
        span
        chord
        engineMass
        density
        liftCoefficient
        volume
        weight
        lift
        lambda
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
            obj.gravity         = cParams.g;
            obj.span            = cParams.b;
            obj.chord           = cParams.c;
            obj.engineMass      = cParams.We;
            obj.density         = cParams.rho;
            obj.volume          = cParams.V;
            obj.lambda          = cParams.lambda;
            obj.liftCoefficient = cParams.Cl;
            
            forces = computeLiftAndWeight(obj);
            obj.lift = forces.liftIntegral;
            obj.weight = forces.weightForceDist;
        end

        function computeForce(obj)
            L  = obj.lift;
            W  = obj.weight;
            b  = obj.span;
            We = obj.engineMass;
            g  = obj.gravity;

            obj.Sz = double(L-W*b/1000 - We*g);
        end

        function forces = computeLiftAndWeight(obj)
            s.g = obj.gravity; 
            s.b = obj.span; 
            s.c = obj.chord;
            s.rho = obj.density;
            s.V = obj.volume;
            s.lambda = obj.lambda; 
            s.Cl = obj.liftCoefficient;

            forces = LiftWeightComputer(s);
            forces.compute();

        end
        
    end
    
end
