classdef FluxComputer < handle

    properties (Access = public)
        q
        qSection
    end

    properties (Access = private)
        Sz 
        inertiaX
        inertiaZ
        thickness
        height 
        airfolDist
        xs 
        xc
        q0
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
            obj.height(1) = cParams.h1;
            obj.height(2) = cParams.h2;
            obj.airfolDist = cParams.d;
            obj.xc = cParams.xc;
            obj.xs = cParams.xs;
            obj.inertiaX = cParams.Ixx;
            obj.inertiaZ = cParams.Izz;
            obj.Sz = cParams.Sz;
            s = computeFirstFlux(obj);
            obj.q0 = s.q0;
        end

        function computeFlux(obj)
            h1 = obj.height(1);
            h2 = obj.height(2);
            d  = obj.airfolDist;
            obj.q(1) = double(int(obj.qSection(1),[0 h1/2]));
            obj.q(2) = double(int(obj.qSection(2),[0 d]));
            obj.q(3) = double(int(obj.qSection(3),[0 h2]));
            obj.q(4) = double(int(obj.qSection(4),[0 d]));
            obj.q(5) = double(int(obj.qSection(5),[0 h1/2]));
        end

         function computeSectionFlux(obj)
            syms y;
            t1 = obj.thickness(1);
            t2 = obj.thickness(2);
            t3 = obj.thickness(3);
            Fz = obj.Sz;
            Ixx = obj.inertiaX;
            d = obj.airfolDist;
            h1 = obj.height(1);
            h2 = obj.height(2);
            obj.qSection = sym('qSection', [1,5]);
            obj.qSection(1) = -(Fz/Ixx)*int(t1*y, y) + obj.q0(1);
            obj.qSection(2) = -(Fz/Ixx)*int(t3* (h1/2 - y*(h1-h2)/(2*d)), y) + obj.q0(2);
            obj.qSection(3) = -(Fz/Ixx)*int(t2* (h2/2 - y), y) + obj.q0(3);
            obj.qSection(4) = -(Fz/Ixx)*int(t3* (-h2/2 - y*(h1-h2)/(2*d)), y) + obj.q0(4) ;
            obj.qSection(5) = -(Fz/Ixx)*int(t1*(-h1/2 + y), y) + obj.q0(5);
        end

        function s = computeFirstFlux(obj)
            t1 = obj.thickness(1);
            t2 = obj.thickness(2);
            t3 = obj.thickness(3);
            Fz = obj.Sz;
            Ixx = obj.inertiaX;
            d = obj.airfolDist;
            h1 = obj.height(1);
            h2 = obj.height(2);
            s.q0(1) = 0;
            s.q0(2) = -(Fz/Ixx) * integral(@(y) t1*y, 0, h1/2) + s.q0(1);
            s.q0(3) = -(Fz/Ixx) * integral(@(y) t3 * (h1/2 - y*(h1-h2)/(2*d)), 0, d) + s.q0(2);
            s.q0(4) = -(Fz/Ixx) * integral(@(y) t2 * (h2/2 - y), 0, h2) + s.q0(3);
            s.q0(5) = -(Fz/Ixx) * integral(@(y) t3 * (-h2/2 - y*(h1-h2)/(2*d)), 0, d) + s.q0(4);
        end

    end
end