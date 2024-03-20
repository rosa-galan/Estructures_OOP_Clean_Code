classdef StressComputer < handle

    properties (Access = public)
        tanStress
        tanMaxStess
        normStress
    end

    properties (Access = private)
        flux
        thickness
        bendingMoment
        liftMoment
        weightMoment
        span 
        engineMass 
        gravity
        inertiaX
        inertiaZ
        beDist 
        Ma
        torsionalInertia
        height
        airfolDist
        xcDist
        xsDist
        width
        forceZ
    end

    methods (Access = public)
        
        function obj = StressComputer(cParams)
            obj.init(cParams);
        end

        function compute(obj)
            obj.computeStress();
        end
        
    end

    methods (Access = private)

        function init(obj,cParams)
            obj.thickness(1) = cParams.t1; 
            obj.thickness(2) = cParams.t2; 
            obj.thickness(3) = cParams.t3; 
            obj.liftMoment   = cParams.liftMomentDist;
            obj.weightMoment = cParams.weightMomentDist;
            obj.span         = cParams.b;
            obj.engineMass   = cParams.We;
            obj.gravity      = cParams.g;
            obj.beDist       = cParams.be; 
            obj.Ma           = cParams.Ma;
            obj.height(1)    = cParams.h1;
            obj.height(2)    = cParams.h2; 
            obj.airfolDist   = cParams.d;
            obj.xsDist       = cParams.xs;
            obj.xcDist       = cParams.xc;
            obj.width        = cParams.a;
            obj.forceZ       = cParams.Sz;

            inertias = computeInertias(obj);
            obj.inertiaX = inertias.IxxT;
            obj.inertiaZ = inertias.IzzT;
            obj.torsionalInertia = inertias.J;

            q = computeFlux(obj);
            obj.flux = q.q;

            obj.bendingMoment = computeBendingMoment(obj);
        end

        function computeStress(obj)
            t1  = obj.thickness(1);
            t2  = obj.thickness(2);
            t3  = obj.thickness(3);
            q   = obj.flux;
            Izz = obj.inertiaZ;
            Ixx = obj.inertiaX;
            Ixz = 0;
            bM  = obj.bendingMoment;
            M   = obj.Ma;
            J   = obj.torsionalInertia;

            obj.tanStress(1) = q(1) / t1;
            obj.tanStress(2) = q(2) / t3;
            obj.tanStress(3) = q(3) / t2;
            obj.tanStress(4) = q(4) / t3;
            obj.tanStress(5) = q(5) / t1;
            obj.tanMaxStess  = double(M) *t3/J;
            obj.normStress   = @(x,z) bM * (Izz*z - Ixz*x)/(Ixx*Izz);
        end

        function s = computeBendingMoment(obj)
            L  = obj.liftMoment;
            b  = obj.span;
            W  = obj.weightMoment;
            Me = obj.engineMass;
            be = obj.beDist;
            g  = obj.gravity;

            s = int(L,[0, b/1000])*1000 - int(W, [0, b/1000])*1000 - Me*g*be;
        end

         function inertias = computeInertias(obj)
            s.t1  = obj.thickness(1);
            s.t2  = obj.thickness(2);
            s.t3  = obj.thickness(3);
            s.d   = obj.airfolDist;
            s.h1  = obj.height(1);
            s.h2  = obj.height(2);
            s.xc  = obj.xcDist;
            s.xs  = obj.xsDist; 
            s.a   = obj.width;

            inertias = InertiaComputer(s);
            inertias.compute();  
         end

         function flux = computeFlux(obj)
            s.t1 = obj.thickness(1);
            s.t2 = obj.thickness(2);
            s.t3 = obj.thickness(3);
            s.h1 = obj.height(1);
            s.h2 = obj.height(2);
            s.d  = obj.airfolDist;
            s.a  = obj.width;
            s.xc = obj.xcDist;
            s.xs = obj.xsDist;
            s.Sz = obj.forceZ;

            flux = FluxComputer(s);
            flux.compute();
        end

    end

end