classdef DeflexionComputer < handle

    properties (Access = public)
        defX
        defY
        defTheta
    end

    properties (Access = private)
        Nel
    end

    methods (Access = public)

        function obj = DeflexionComputer(cParams)
            obj.init(cParams);
        end

        function compute(obj)
            obj.computeDeflexion();
        end
        
    end

    methods (Access = private)

        function init(obj,cParams)
            obj.Nel = cParams.Nel;
        end

        function computeDeflexion(obj)
            obj.defX = 1:3:3*obj.Nel;
            obj.defY = obj.defX + 1;
            obj.defTheta = obj.defY + 1;
        end

    end
    
end