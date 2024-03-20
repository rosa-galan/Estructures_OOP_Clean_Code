classdef LiftWeightComputer < handle

    properties (Access = public)
        liftForceDist 
        liftIntegral
        liftMomentDist
        weightForceDist
        weightMomentDist
    end

    properties (Access = private)
        density
        volume
        liftCoefficient
        span
        chord
    end

    methods (Access = public)
        
        function obj = LiftWeightComputer(cParams)
            obj.init(cParams);  
        end

        function compute(obj)
            obj.computeLiftDistribution();
            obj.computeLiftForce();
            obj.computeLiftMoment();
            obj.computeWeightMoment();
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.density         = cParams.rho;
            obj.volume          = cParams.V;
            obj.liftCoefficient = cParams.Cl;
            obj.span            = cParams.b;
            obj.chord           = cParams.c;
            obj.weightForceDist = cParams.lambda * cParams.g;
        end
         
        function computeLiftDistribution(obj)
            syms y;
            rho = obj.density;
            V   = obj.volume;
            c   = obj.chord;
            Cl  = obj.liftCoefficient;
            b   = obj.span;

            obj.liftForceDist = 0.5 * rho * V^2 * (c/1000) * Cl * sqrt(1 - (y / (b/1000))^2);
        end

        function computeLiftForce(obj)
            b = obj.span;
            obj.liftIntegral = int(obj.liftForceDist, [0 b/1000]);           
        end

        function computeLiftMoment(obj)
            syms y;
            obj.liftMomentDist = y * obj.liftForceDist;
        end

        function computeWeightMoment(obj)
            syms y;
            obj.weightMomentDist = y * obj.weightForceDist;
        end   

    end

end