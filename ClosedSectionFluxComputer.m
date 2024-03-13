classdef ClosedSectionFluxComputer < handle

    properties (Access = public)
        qClosed
        Ma
    end

    properties (Access = private)
        Sz 
        thickness
        height 
        airfolDist 
        xsDist 
        xcDist
        qSection
        closedMoment
        distShearC
        interiorS
    end

    methods (Access = public)
        
        function obj = ClosedSectionFluxComputer(cParams)
            obj.init(cParams);
        end

        function compute(obj)
            obj.computeClosedSectionFlux();
            obj.computeMomentSum();
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.thickness(1) = cParams.t1;
            obj.thickness(2) = cParams.t2;
            obj.thickness(3) = cParams.t3;
            obj.height(1)    = cParams.h1;
            obj.height(2) = cParams.h2;
            obj.airfolDist = cParams.d;
            obj.xcDist = cParams.xc;
            obj.xsDist = cParams.xs;
            obj.qSection = cParams.qSection;
            obj.Sz = cParams.Sz;
            obj.closedMoment = computeClosedMoment(obj);
            obj.Ma = computeMomentSum(obj);
            obj.distShearC = double(obj.Ma)/obj.Sz;
            obj.interiorS = obj.height(2)*obj.airfolDist + obj.airfolDist*(obj.height(1)-obj.height(2))/2;
        end

        function computeClosedSectionFlux(obj)
            xsc = obj.distShearC;
            xc = obj.xcDist;
            Fz = obj.Sz;
            q = obj.qSection;
            h1 = obj.height(1);
            h2 = obj.height(2);
            d = obj.airfolDist;
            obj.qClosed = ((xsc-xc)*Fz - (q(1)+q(2)+q(3)+q(4)+q(5)))/2*(d*(h1+h2)/2);
        end

        function s = computeClosedMoment(obj)
            syms y
            q = obj.qSection;
            h1 = obj.height(1);
            h2 = obj.height(2);
            d = obj.airfolDist;
            xs = obj.xsDist;
            xc = obj.xcDist;
            s = int(q(1).*(xs-xc),[0.0, h1/2]) + int(q(2).*(h1/2 - y*(h1-h2)/(2*d)),[0, d]) + ...
            int(q(3).*(xc-d),[0 , h2]) + int(q(4).*(-h2/2 - y*(h1-h2)/(2*d)),[0 , d]) + int(q(5).*(xs-xc),[0 , h1/2]);
        end

        function s = computeMomentSum(obj)
            syms y
            q = obj.qSection;
            h1 = obj.height(1);
            h2 = obj.height(2);
            d = obj.airfolDist;
            t1 = obj.thickness(1);
            t3 = obj.thickness(3);
            s = int(q(1)*d*t1, [0 h1/2]) + int(q(2)*t3*(h1/2 - y*(h1-h2)/(2*d)),[0 d]) + int(q(3)*0, [h2/2  -h2/2]) + ... 
            int(q(4)*t3*(-h2/2 - y*(h1-h2)/(2*d)), [0 d]) + int(q(5)*d*t1, [0 h1/2]);
        end

    end

end