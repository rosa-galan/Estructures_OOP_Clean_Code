classdef ConnectDOFComputer < handle
    
    properties (Access = public)
        Td
    end

    properties (Access = private)
        dimensions
        nodalConnectivities
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
            obj.dimensions          = cParams.dim;
            obj.nodalConnectivities = cParams.Tn;
        end

        function computeConnectMatrix(obj)
            dim = obj.dimensions;
           
            obj.Td=zeros(dim.nel, dim.ni*dim.nne);

            for nElement=1:dim.nel
                for nNode=1:dim.nne
                    for nDegree=1:dim.ni
                        
                        obj.Td(nElement, dim.ni*(nNode-1)+nDegree) = ...
                        computeNod2Dof(obj, nElement, nNode, nDegree);
                    end
                end
            end
        end

        function s = computeNod2Dof(obj, nElement, nNode, nDegree)
            dim = obj.dimensions;
            Tn = obj.nodalConnectivities;

            s = dim.ni*Tn(nElement,nNode) - dim.ni + nDegree;
        end
        
    end

end
