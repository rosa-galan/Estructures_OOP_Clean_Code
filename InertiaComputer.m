classdef InertiaComputer < handle

    properties (Access = public)
        IxxT 
        IzzT 
        J 
    end

    properties (Access = private)
        thickness
        height
        width
        airfolDist
        xcDist
        xsDist
    end

    methods (Access = public)

        function obj = InertiaComputer(cParams)
            obj.init(cParams);
        end

        function compute(obj)
            obj.computeInertia();
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.thickness(1) = cParams.t1;
            obj.thickness(2) = cParams.t2;
            obj.thickness(3) = cParams.t3;
            obj.height(1)    = cParams.h1;
            obj.height(2)    = cParams.h2;
            obj.width        = cParams.a;
            obj.airfolDist   = cParams.d;
            obj.xcDist       = cParams.xc;
            obj.xsDist       = cParams.xs;
        end

        function computeInertia(obj)
            Ixx = obj.computeInertiaX();
            obj.IxxT = Ixx(1)+Ixx(2)+2*Ixx(3);

            Izz = obj.computeInertiaZ();
            obj.IzzT = Izz(1)+Izz(2)+2*Izz(3);

            obj.J = obj.computeTorsionalInertia();
            
        end

        function inertiaX = computeInertiaX(obj)
            t1 = obj.thickness(1);
            t2 = obj.thickness(2);
            t3 = obj.thickness(3);
            h1 = obj.height(1);
            h2 = obj.height(2);
            a  = obj.width;
            d  = obj.airfolDist;
            
            inertiaX(1) = (t1*h1^3)/12;
            inertiaX(2) = (t2*h2^3)/12;
            inertiaX(3) = a^3*sind(4.76)^2*t3/12 + (h2/2+d/2*tand(4.76))^2*t3*a;
        end

        function inertiaZ = computeInertiaZ(obj)
            t1 = obj.thickness(1);
            t2 = obj.thickness(2);
            t3 = obj.thickness(3);
            h1 = obj.height(1);
            h2 = obj.height(2);
            a  = obj.width;
            d  = obj.airfolDist;
            xs = obj.xsDist;
            xc = obj.xcDist;

            inertiaZ(1) = (xs-xc)^2*t1*h1;
            inertiaZ(2) = (xs-xc+d)^2*t2*h2;
            inertiaZ(3) = a*d^2*t3/12 + (xs+d/2-xc)^2*t3*a;
        end

        function torsionalInertia = computeTorsionalInertia(obj)
            t1 = obj.thickness(1);
            t2 = obj.thickness(2);
            t3 = obj.thickness(3);
            h1 = obj.height(1);
            h2 = obj.height(2);
            a  = obj.width;

            torsionalInertia = 1/3*(t2^3*h2+t1^3*h1+2*a*t3^3);
        end

    end
    
end