classdef InternalForcesComputer < handle

    properties (Access = public)
        Q
        Mb
        Mt
    end

    properties (Access = private)
        dimensions
        nodalConnectivities
        displacement
        xDist
        matConnectivities
        matProp
    end

    methods (Access = public)

        function obj = InternalForcesComputer(cParams)
            obj.init(cParams);
        end

        function compute(obj)
            obj.computeInternalForces();
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.xDist               = cParams.x;
            obj.matConnectivities   = cParams.Tm;
            obj.matProp             = cParams.matProp;
            obj.nodalConnectivities = cParams.Tn;
            obj.dimensions          = cParams.dim;
            obj.displacement        = cParams.u;
        end

        function computeInternalForces(obj)

            connectDOF = computeDOFConnectivities(obj);
            Td = connectDOF.Td;

            elementalMatrix = computeElementalMatrix(obj);
            Kel = elementalMatrix.Kel;

            dim = obj.dimensions;
            u = obj.displacement;

            obj.Q = zeros(dim.nel,2);
            obj.Mb = zeros(dim.nel,2);
            obj.Mt = zeros(dim.nel,2);

            for e=1:dim.nel
                feInt = Kel(:,:,e)*u(Td(e,:));
                obj.Q(e,1)  = -feInt(1);
                obj.Mb(e,1) = -feInt(2);
                obj.Mt(e,1) = -feInt(3);
                obj.Q(e,2)  =  feInt(4);
                obj.Mb(e,2) =  feInt(5);
                obj.Mt(e,2) =  feInt(6);
            end
        end

        function elementalMatrix = computeElementalMatrix(obj)

            s.x = obj.xDist;
            s.dim = obj.dimensions;
            s.Tn = obj.nodalConnectivities;
            s.matProp = obj.matProp;
            s.Tm = obj.matConnectivities;

            elementalMatrix = ElementalStiffnessMatrixComputer(s);
            elementalMatrix.compute();
        end

        function connectDOF = computeDOFConnectivities(obj)

            s.dim = obj.dimensions;
            s.Tn = obj.nodalConnectivities;

            connectDOF = ConnectDOFComputer(s);
            connectDOF.compute();
        end

    end  
    
end