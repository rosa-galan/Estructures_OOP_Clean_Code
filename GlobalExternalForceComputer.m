classdef GlobalExternalForceComputer < handle

    properties (Access = public)
        Fext
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
        fdata
    end

    methods (Access = public)

        function obj = GlobalExternalForceComputer(cParams)
            obj.init(cParams);
        end

        function compute(obj)
            obj.computeExternalForce();
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.dimensions          = cParams.dim;
            obj.fdata               = cParams.fdata;
            obj.nodalConnectivities = cParams.Tn;
            obj.xaDist              = cParams.xa;
            obj.xcDist              = cParams.xc;
            obj.xmDist              = cParams.xm;
            obj.lift                = cParams.liftForceDist;
            obj.weight              = cParams.weightForceDist;
            obj.xDist               = cParams.x;
        end

        function computeExternalForce(obj)
            
            connectMatrix = computeConnectMatrix(obj);
            Td = connectMatrix.Td;

            elementalForce = computeElementalForce(obj);
            fel = elementalForce.fel;

            dim = obj.dimensions;
            fd  = obj.fdata;

            obj.Fext = zeros(dim.ndof,1);

            for e = 1:dim.nel
                for i = 1:dim.nne*dim.ni
                    I = Td(e,i);
                    obj.Fext(I) = obj.Fext(I) + fel(i,e);
                end
            end
            for i = 1:size(fd,1)
                I = computeNod2Dof(obj,i);
                obj.Fext(I) = obj.Fext(I) + fd(i,3);
            end
        end

        function s = computeNod2Dof(obj,i)
            dim = obj.dimensions;
            fd  = obj.fdata;

            s = dim.ni*fd(i,1) - dim.ni + fd(i,2);
        end

        function connectMatrix = computeConnectMatrix(obj)
              s.dim = obj.dimensions;
              s.Tn  = obj.nodalConnectivities;

              connectMatrix = ConnectDOFComputer(s);
              connectMatrix.compute();
        end

        function elementalForce = computeElementalForce(obj)
            s.dim             = obj.dimensions;
            s.x               = obj.xDist;
            s.Tn              = obj.nodalConnectivities;
            s.xa              = obj.xaDist;
            s.xc              = obj.xcDist;
            s.xm              = obj.xmDist;
            s.liftForceDist   = obj.lift;
            s.weightForceDist = obj.weight;

            elementalForce = ElementalForceComputer(s);
            elementalForce.compute();
        end

    end

end