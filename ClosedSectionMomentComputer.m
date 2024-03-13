classdef ClosedSectionMomentComputer < handle
    
    properties (Access = public)
        closedMsc 
    end

    properties (Access = private)
        Mengine
        Mlift
        Mweight
        xm
        xc
        xa
        b
        Me
        g
        xe
        liftForceDist
        weightForceDist
    end

    methods (Access = public)
        
        function obj = ClosedSectionMomentComputer(cParams)
            obj.init(cParams);
        end

        function compute(obj)
            obj.computeMoment();
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.xe = cParams.xe;
            obj.xm = cParams.xm;
            obj.xc = cParams.xc;
            obj.xa = cParams.xa;
            obj.b = cParams.b;
            obj.g = cParams.g;
            obj.Me = cParams.Me;
            obj.liftForceDist = cParams.liftForceDist;
            obj.weightForceDist = cParams.weightForceDist;
            obj.Mengine = obj.Me*obj.g*(obj.xc/1000-(obj.xe/1000-obj.xc/1000));
            obj.Mlift = int(obj.liftForceDist, [0 obj.b/1000])*(obj.xc/1000-(obj.xc/1000-obj.xa/1000));
            obj.Mweight = obj.weightForceDist*16*(obj.xc/1000-(obj.xm/1000-obj.xc/1000));
        end

        function computeMoment(obj)
            obj.closedMsc = obj.Mengine + obj.Mlift + obj.Mweight;
        end
        
    end

end
