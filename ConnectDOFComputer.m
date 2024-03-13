classdef ConnectDOFComputer < handle
    
    properties (Access = public)
        Td
    end

    properties (Access = private)
        dim
        Tn
    end

    methods (Access = public)
        
        function obj = ConnectDOFComputer(cParams)
            obj.init(cParams);
        end

        function compute(obj)
            obj.computeConnectMatrix();
        end

    end

    methods  (Access = private)

        function init(obj,cParams)
            obj.dim = cParams.dim;
            obj.Tn = cParams.Tn;
        end

        function computeConnectMatrix(obj)
            obj.Td=zeros(obj.dim.nel,obj.dim.ni*obj.dim.nne);
            for nElement=1:obj.dim.nel
                for nNode=1:obj.dim.nne
                    for nDegree=1:obj.dim.ni
                        obj.Td(nElement, obj.dim.ni*(nNode-1)+nDegree) = ...
                        computeNod2Dof(obj, nElement, nNode, nDegree);
                    end
                end
            end
        end

        function s = computeNod2Dof(obj, nElement, nNode, nDegree)
            s = obj.dim.ni*obj.Tn(nElement,nNode) - obj.dim.ni + nDegree;
        end
        
    end

end
