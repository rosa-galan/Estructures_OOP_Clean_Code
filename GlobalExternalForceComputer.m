classdef GlobalExternalForceComputer < handle

    properties (Access = public)
        Fext
    end

    properties (Access = private)
        dim
        Td
        fel 
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
            obj.dim = cParams.dim;
            obj.Td = cParams.Td;
            obj.fel = cParams.fel;
            obj.fdata = cParams.fdata;
        end

        function computeExternalForce(obj)
            obj.Fext = zeros(obj.dim.ndof,1);
            for e = 1:obj.dim.nel
                for i = 1:obj.dim.nne*obj.dim.ni
                    I = obj.Td(e,i);
                    obj.Fext(I) = obj.Fext(I) + obj.fel(i,e);
                end
            end
            for i = 1:size(obj.fdata,1)
                I = computeNod2Dof(obj,i);
                obj.Fext(I) = obj.Fext(I) + obj.fdata(i,3);
            end
        end

        function s = computeNod2Dof(obj,i)
            s = obj.dim.ni*obj.fdata(i,1) - obj.dim.ni + obj.fdata(i,2);
        end

    end

end