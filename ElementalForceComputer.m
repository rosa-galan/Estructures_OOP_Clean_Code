classdef ElementalForceComputer < handle

    properties (Access = public)
        fel
    end

    properties (Access = private)
        dimensions
        xDist
        nodalConnectivities
        lift 
        weight
        xaDist
        xcDist
        xmDist
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
            obj.dimensions          = cParams.dim;
            obj.xDist               = cParams.x;
            obj.nodalConnectivities = cParams.Tn;
            obj.xaDist              = cParams.xa;
            obj.xcDist              = cParams.xc;
            obj.xmDist              = cParams.xm;
            obj.lift                = cParams.liftForceDist;
            obj.weight              = cParams.weightForceDist;
        end

        function computeElementalForce(obj)
            
            elementalForces = computeForces(obj);
            bendingFel = elementalForces.bending;
            torsionFel = elementalForces.torsion;

            dim = obj.dimensions;

            obj.fel = zeros(dim.ni*dim.nne, dim.nel);
            for e=1:dim.nel
                obj.fel(:,e) = bendingFel(:,e) + torsionFel(:,e);
            end
        end

        function elemForces = computeForces(obj)
            syms y 
            liftInt = int(obj.lift, y);
            weightInt = int(obj.weight, y);
            
            dim = obj.dimensions;
            x = obj.xDist;
            Tn = obj.nodalConnectivities;
            xa = obj.xaDist;
            xc = obj.xcDist;
            xm = obj.xmDist;
            
            for e=1:dim.nel
                le = abs(x(Tn(e,2),1)- x(Tn(e,1),1))/1000;
                intStart   = x(Tn(e,1),1)/1000;
                intEnd     = x(Tn(e,2),1)/1000;
                qe = double(subs(liftInt -weightInt, y, intEnd) - subs(liftInt - weightInt, y, intStart))/le;
                mbe = 0;
                msce = double(subs(liftInt, y, intEnd) - subs(liftInt, y, intStart))/le * (xa-xc)/1000 - ...
                double(subs(weightInt, y, intEnd) - subs(weightInt, y, intStart))/le * (xm-xc)/1000;
                
                elemForces.bending(:,e) = qe*le*[0.5; le/12; 0; 0.5; -le/12; 0] + mbe*le*[0; 0.5; 0; 0; 0.5; 0];
                elemForces.torsion(:,e) = msce*le*[0; 0; 0.5; 0; 0; 0.5];
            end
        end

    end

end