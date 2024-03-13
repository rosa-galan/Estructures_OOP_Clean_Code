classdef FixedDOFComputer < handle

    properties (Access = public)
        ur
        vr
        vf
    end

    properties (Access = private)
        dim
        fixNodes
        allDOF
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
            obj.dim = cParams.dim;
            obj.fixNodes = cParams.fixNodes;
            obj.allDOF = 1:obj.dim.ndof;
        end

        function computeFixedDOF(obj)
            obj.vr = zeros(1,size(obj.fixNodes,1));
            obj.ur = zeros(1,size(obj.fixNodes,1)); 
            for i=1:(size(obj.fixNodes, 1))
                fnd = obj.fixNodes(i,:);
                I = computeNod2Dof(obj, fnd);
                obj.vr(i) = I;
                obj.ur(i) = fnd(3);
            end
            obj.vf = setdiff(obj.allDOF, obj.vr);
        end

        function s = computeNod2Dof(obj,fnd)
            s = obj.dim.ni*fnd(1) - obj.dim.ni + fnd(2);
        end 

    end

end
