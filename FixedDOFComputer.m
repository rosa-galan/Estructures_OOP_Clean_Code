classdef FixedDOFComputer < handle

    properties (Access = public)
        ur
        vr
        vf
    end

    properties (Access = private)
        dimensions
        fixNodes
    end

    methods (Access = public)

        function obj = FixedDOFComputer(cParams)
            obj.init(cParams);
        end

        function compute(obj)
            obj.computeFixedDOF();
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.dimensions = cParams.dim;
            obj.fixNodes = cParams.fixNodes;
        end

        function computeFixedDOF(obj)
            
            dim = obj.dimensions;
            DOF = 1:dim.ndof;
            fN = obj.fixNodes;
            
            obj.vr = zeros(1,size(fN, 1));
            obj.ur = zeros(1,size(fN, 1)); 

            for i=1:(size(fN, 1))
                fnd = fN(i,:);
                I = computeNod2Dof(obj, fnd);
                obj.vr(i) = I;
                obj.ur(i) = fnd(3);
            end
            obj.vf = setdiff(DOF, obj.vr);
        end

        function s = computeNod2Dof(obj,fnd)
            dim = obj.dimensions;

            s = dim.ni*fnd(1) - dim.ni + fnd(2);
        end 

    end

end
