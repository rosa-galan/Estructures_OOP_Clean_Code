classdef FluxComputer < handle

    properties (Access = public)
        q
        qSection
    end

    properties (Access = private)
        forceZ 
        inertiaX 
        inertiaZ 
        thickness
        height 
        airfolDist
        xsDist 
        xcDist
        startFlux
        width
    end

    methods (Access = public)
        
        function obj = FluxComputer(cParams)
            obj.init(cParams); 
        end

        function compute(obj)
            obj.computeSectionFlux();
            obj.computeFlux();
        end

    end

    methods (Access = private)
        
        function init(obj,cParams)
            obj.thickness(1) = cParams.t1;
            obj.thickness(2) = cParams.t2;
            obj.thickness(3) = cParams.t3;
            obj.height(1)    = cParams.h1;
            obj.height(2)    = cParams.h2;
            obj.airfolDist   = cParams.d;
            obj.xcDist       = cParams.xc;
            obj.xsDist       = cParams.xs;
            obj.width        = cParams.a;
            obj.forceZ       = cParams.Sz;
            
            inertias     = computeInertias(obj);
            obj.inertiaX = inertias.IxxT;
            obj.inertiaZ = inertias.IzzT;

            % force      = computeForce(cParams);
            % force.compute();
            % obj.forceZ = force.Sz;

            obj.startFlux     = computeFirstFlux(obj);
            
        end

        function computeFlux(obj)
            h1 = obj.height(1);
            h2 = obj.height(2);
            d  = obj.airfolDist;
            qS = obj.qSection;
            
            obj.q(1) = double(int(qS(1), [0 h1/2]));
            obj.q(2) = double(int(qS(2), [0 d]));
            obj.q(3) = double(int(qS(3), [0 h2]));
            obj.q(4) = double(int(qS(4), [0 d]));
            obj.q(5) = double(int(qS(5), [0 h1/2]));
        end

         function computeSectionFlux(obj)
            syms y;
            t1  = obj.thickness(1);
            t2  = obj.thickness(2);
            t3  = obj.thickness(3);
            Sz  = obj.forceZ;
            Ixx = obj.inertiaX;
            d   = obj.airfolDist;
            h1  = obj.height(1);
            h2  = obj.height(2);
            q0  = obj.startFlux;
            
            obj.qSection = sym('qSection', [1,5]);
            obj.qSection(1) = -(Sz/Ixx)*int(t1*y, y) + q0(1);
            obj.qSection(2) = -(Sz/Ixx)*int(t3* (h1/2 - y*(h1-h2)/(2*d)), y) + q0(2);
            obj.qSection(3) = -(Sz/Ixx)*int(t2* (h2/2 - y), y) + q0(3);
            obj.qSection(4) = -(Sz/Ixx)*int(t3* (-h2/2 - y*(h1-h2)/(2*d)), y) + q0(4) ;
            obj.qSection(5) = -(Sz/Ixx)*int(t1*(-h1/2 + y), y) + q0(5);
        end

        function q0 = computeFirstFlux(obj)
            t1  = obj.thickness(1);
            t2  = obj.thickness(2);
            t3  = obj.thickness(3);
            Sz  = obj.forceZ;
            Ixx = obj.inertiaX;
            d   = obj.airfolDist;
            h1  = obj.height(1);
            h2  = obj.height(2);
            
            q0(1) = 0;
            q0(2) = -(Sz/Ixx) * integral(@(y) t1*y, 0, h1/2) + q0(1);
            q0(3) = -(Sz/Ixx) * integral(@(y) t3 * (h1/2 - y*(h1-h2)/(2*d)), 0, d) + q0(2);
            q0(4) = -(Sz/Ixx) * integral(@(y) t2 * (h2/2 - y), 0, h2) + q0(3);
            q0(5) = -(Sz/Ixx) * integral(@(y) t3 * (-h2/2 - y*(h1-h2)/(2*d)), 0, d) + q0(4);
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

        % function force = computeForce(cParams)
        %     s.ForceZComputer(cParams);
        %     s.compute();
        %     force = s.Sz;
        % end

    end
end