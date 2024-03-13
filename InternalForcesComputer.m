classdef InternalForcesComputer < handle

    properties (Access = public)
        Q
        Mb
        Mt
    end

    properties (Access = private)
        dim
        Td
        Kel
        u
        feInt
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
            obj.dim = cParams.dim;
            obj.Td = cParams.Td;
            obj.Kel = cParams.Kel;
            obj.u = cParams.u;
        end

        function computeInternalForces(obj)
            obj.Q = zeros(obj.dim.nel,2);
            obj.Mb = zeros(obj.dim.nel,2);
            obj.Mt = zeros(obj.dim.nel,2);
            for e=1:obj.dim.nel
                obj.feInt = obj.Kel(:,:,e) * obj.u(obj.Td(e,:));
                obj.Q(e,1)  = -obj.feInt(1);
                obj.Mb(e,1) = -obj.feInt(2);
                obj.Mt(e,1) = -obj.feInt(3);
                obj.Q(e,2)  =  obj.feInt(4);
                obj.Mb(e,2) =  obj.feInt(5);
                obj.Mt(e,2) =  obj.feInt(6);
            end
        end

    end  
    
end