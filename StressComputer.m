classdef StressComputer < handle

    properties (Access = public)
        tanStress
        tanMaxStess
        normStress
    end

    properties (Access = private)
        q
        t1 
        t2 
        t3 
        bendingMoment
        liftMomentDist
        weightMomentDist
        b 
        Me 
        g 
        Ixx 
        Izz
        be 
        Ma
        J
    end

    methods (Access = public)
        
        function obj = StressComputer(cParams)
            obj.init(cParams);
        end

        function compute(obj)
            obj.computeStress();
        end
        
    end

    methods (Access = private)

        function init(obj,cParams)
            obj.q = cParams.q;
            obj.t1 = cParams.t1; 
            obj.t2 = cParams.t2; 
            obj.t3 = cParams.t3; 
            obj.liftMomentDist = cParams.liftMomentDist;
            obj.weightMomentDist = cParams.weightMomentDist;
            obj.b = cParams.b;
            obj.Me = cParams.Me;
            obj.g = cParams.g;
            obj.Ixx = cParams.Ixx; 
            obj.Izz = cParams.Izz;
            obj.be = cParams.be; 
            obj.Ma = cParams.Ma;
            obj.J = cParams.J;
            obj.bendingMoment = computeBendingMoment(obj);
        end

        function computeStress(obj)
            obj.tanStress(1) = obj.q(1) / obj.t1;
            obj.tanStress(2) = obj.q(2) / obj.t3;
            obj.tanStress(3) = obj.q(3) / obj.t2;
            obj.tanStress(4) = obj.q(4) / obj.t3;
            obj.tanStress(5) = obj.q(5) / obj.t1;
            obj.tanMaxStess = double(obj.Ma) * obj.t3/obj.J;
            obj.normStress = @(x,z) obj.bendingMoment * (obj.Izz*z - Ixz*x)/(obj.Ixx*obj.Izz);
        end

        function s = computeBendingMoment(obj)
            s = int(obj.liftMomentDist,[0, obj.b/1000])*1000 - int(obj.weightMomentDist, [0, obj.b/1000])*1000 - obj.Me*obj.g*obj.be;
        end

    end

end