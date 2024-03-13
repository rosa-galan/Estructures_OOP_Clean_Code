classdef InertiaComputer < handle

    properties (Access = public)
        Ixx 
        IxxT 
        Izz
        IzzT 
        J 
    end

    properties (Access = private)
        t1
        t2
        t3
        h1
        h2
        a
        d
        xc
        xs
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
            obj.t1 = cParams.t1;
            obj.t2 = cParams.t2;
            obj.t3 = cParams.t3;
            obj.h1 = cParams.h1;
            obj.h2 = cParams.h2;
            obj.a = cParams.a;
            obj.d = cParams.d;
            obj.xc = cParams.xc;
            obj.xs = cParams.xs;
        end

        function computeInertia(obj)
            obj.Ixx(1) = (obj.t1*obj.h1^3)/12;
            obj.Ixx(2) = (obj.t2*obj.h2^3)/12;
            obj.Ixx(3) = obj.a^3*sind(4.76)^2*obj.t3/12 + (obj.h2/2+obj.d/2*tand(4.76))^2*obj.t3*obj.a;
            obj.IxxT = obj.Ixx(1)+obj.Ixx(2)+2*obj.Ixx(3);
            obj.Izz(1) = (obj.xs-obj.xc)^2 * obj.t1*obj.h1;
            obj.Izz(2) = (obj.xs-obj.xc+obj.d)^2 * obj.t2*obj.h2;
            obj.Izz(3) = obj.a*obj.d^2*obj.t3/12 + (obj.xs + obj.d/2 - obj.xc)^2 * obj.t3*obj.a;
            obj.IzzT = obj.Izz(1)+obj.Izz(2)+2*obj.Izz(3);
            obj.J = 1/3*(obj.t2^3*obj.h2+obj.t1^3*obj.h1+2*obj.a*obj.t3^3);
        end

    end
    
end