classdef OpenSectionMomentComputer < handle
    
    properties (Access = public)
        openMsc 
    end

    properties (Access = private)
        xmDist
        shearCenter
        xcDist
        xaDist
        span
        engineMass
        gravity
        xeDist
        lift
        weight
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
            obj.xeDist      = cParams.xe;
            obj.xmDist      = cParams.xm;
            obj.xcDist      = cParams.xc;
            obj.xaDist      = cParams.xa; 
            obj.span        = cParams.b;
            obj.shearCenter = cParams.xsc; % ho deixo aixi si no he de definir masses coses
            obj.gravity     = cParams.g;
            obj.engineMass  = cParams.Me;
            obj.lift        = cParams.liftForceDist; % ho deixo aixi perquè si no s'ha de definir molt
            obj.weight      = cParams.weightForceDist;
        end

        function computeMoment(obj)
            engineMoment = computeEngineMoment(obj);
            liftMoment   = computeLiftMoment(obj);
            weightMoment = computeWeightMoment(obj);

            obj.openMsc = engineMoment + liftMoment + weightMoment;
        end

        function s = computeEngineMoment(obj)
            Me  = obj.engineMass;
            g   = obj.gravity;
            xsc = obj.shearCenter;
            xe  = obj.xeDist;
            xc  = obj.xcDist;

            s = Me*g*(xsc/1000-(xe/1000-xc/1000));
        end

        function s = computeLiftMoment(obj)
            xsc = obj.shearCenter;
            xa  = obj.xaDist;
            xc  = obj.xcDist;
            b   = obj.span;
            L   = obj.lift;

            s = int(L, [0 b/1000])*(xsc/1000-(xc/1000-xa/1000));
        end

        function s = computeWeightMoment(obj)
            xsc = obj.shearCenter;
            xm  = obj.xmDist;
            xc  = obj.xcDist;
            W  = obj.weight;

            s = W*16*(xsc/1000-(xm/1000-xc/1000));
        end
        
    end

end
