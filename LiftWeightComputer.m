classdef LiftWeightComputer < handle

    properties (Access = public)
        liftForceDist 
        liftIntegral
        liftMomentDist
        weightForceDist
        weightMomentDist
    end

    properties (Access = private)
        rho
        V
        Cl
        lambda
        g
        b
        c
        Me
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
            obj.rho = cParams.rho;
            obj.V = cParams.V;
            obj.Cl = cParams.Cl;
            obj.lambda = cParams.lambda;
            obj.g = cParams.g;
            obj.b = cParams.b;
            obj.c = cParams.c;
            obj.Me = cParams.Me;
            obj.weightForceDist = cParams.lambda * cParams.g;
        end
         
        function computeLiftDistribution(obj)
            syms y;
            obj.liftForceDist = 0.5 * obj.rho * obj.V^2 * (obj.c/1000) * obj.Cl * sqrt(1 - (y / (obj.b/1000))^2);
        end

        function computeLiftForce(obj)
            obj.liftIntegral = int(obj.liftForceDist, [0 obj.b/1000]);           
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