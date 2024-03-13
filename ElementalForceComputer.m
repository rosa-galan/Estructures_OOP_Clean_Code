classdef ElementalForceComputer < handle

    properties (Access = public)
        fel
    end

    properties (Access = private)
        dim
        x
        Tn
        liftForceDist  
        liftMomentDist
        weightForceDist
        weightMomentDist
        gProp
        xa
        xc
        xm
    end

    methods (Access = public)

        function obj = ElementalForceComputer(cParams)
            obj.init(cParams)
        end

        function compute(obj)
            obj.computeElementalForce();
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.dim = cParams.dim;
            obj.x = cParams.x;
            obj.Tn = cParams.Tn;
            obj.xa = cParams.xa;
            obj.xc = cParams.xc;
            obj.xm = cParams.xm;
            obj.liftForceDist = cParams.liftForceDist;
            obj.liftMomentDist = cParams.liftMomentDist;
            obj.weightForceDist = cParams.weightForceDist;
            obj.weightMomentDist = cParams.weightMomentDist;
            obj.gProp = computeGeneralProperties(obj);
        end

        function computeElementalForce(obj)
            obj.fel = zeros(obj.dim.ni*obj.dim.nne,obj.dim.nel);
            for e=1:obj.dim.nel
                obj.fel(:,e) = obj.gProp.bendingFel(:,e) + obj.gProp.torsionFel(:,e);
            end
        end

        function s = computeGeneralProperties(obj)
            syms y 
            s.liftInt = int(obj.liftForceDist, y);
            s.weightInt = int(obj.weightForceDist, y);
            for e=1:obj.dim.nel
                s.le = abs(obj.x(obj.Tn(e,2),1)- obj.x(obj.Tn(e,1),1))/1000;
                s.intStart   = obj.x(obj.Tn(e,1),1)/1000;
                s.intEnd     = obj.x(obj.Tn(e,2),1)/1000;
                s.qe = double(subs(s.liftInt - s.weightInt, y, s.intEnd) - subs(s.liftInt - s.weightInt, y, s.intStart))/s.le;
                s.mbe = 0;
                s.msce = double(subs(s.liftInt, y, s.intEnd) - subs(s.liftInt, y, s.intStart))/s.le * (obj.xa-obj.xc)/1000 - ...
                double(subs(s.weightInt, y, s.intEnd) - subs(s.weightInt, y, s.intStart))/s.le * (obj.xm-obj.xc)/1000;
                s.bendingFel(:,e) = s.qe * s.le * [0.5; s.le/12; 0 ; 0.5; -s.le/12; 0] + ...
                s.mbe * s.le * [0  ;    0.5; 0 ;   0;     0.5; 0];
                s.torsionFel(:,e) = s.msce * s.le * [ 0; 0; 0.5; 0; 0; 0.5];
            end
        end

    end

end