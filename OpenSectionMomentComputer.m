classdef OpenSectionMomentComputer < handle
    
    properties (Access = public)
        openMsc 
    end

    properties (Access = private)
        Mengine
        Mlift
        Mweight
        xm
        xsc
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
        
        function obj = OpenSectionMomentComputer(cParams)
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
            obj.xsc = cParams.xsc;
            obj.xa = cParams.xa;
            obj.b = cParams.b;
            obj.g = cParams.g;
            obj.Me = cParams.Me;
            obj.liftForceDist = cParams.liftForceDist;
            obj.weightForceDist = cParams.weightForceDist;
            obj.Mengine = obj.Me*obj.g*(obj.xsc/1000-(obj.xe/1000-obj.xc/1000));
            obj.Mlift = int(obj.liftForceDist, [0 obj.b/1000])*(obj.xsc/1000-(obj.xc/1000-obj.xa/1000));
            obj.Mweight = obj.weightForceDist*16*(obj.xsc/1000-(obj.xm/1000-obj.xc/1000));
        end

        function computeMoment(obj)
            obj.openMsc = obj.Mengine + obj.Mlift + obj.Mweight;
        end
        
    end

end
