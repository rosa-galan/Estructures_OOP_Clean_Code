classdef ClosedSectionFluxComputer < handle

    properties (Access = public)
        qClosed
        Ma
    end

    properties (Access = private)
        forceZ
        width
        thickness
        height 
        airfolDist 
        xsDist 
        xcDist
        closedMoment
        distShearC
        interiorS
        fluxSection
    end

    methods (Access = public)
        
        function obj = ClosedSectionFluxComputer(cParams)
            obj.init(cParams);
        end

        function compute(obj)
            obj.computeMoments();
            obj.computeClosedSectionFlux();
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
            obj.forceZ       = cParams.Sz;
            obj.width        = cParams.a;
            
            flux  = computeFlux(obj);
            obj.fluxSection = flux.qSection;
            
            obj.interiorS    = computeInteriorSection(obj);
        end

        function computeClosedSectionFlux(obj)
            obj.qClosed = -obj.closedMoment/(2*obj.interiorS);
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


        function computeMoments(obj)
            s1 = computeFirstClosedMoment(obj);
            s2 = computeSecondClosedMoment(obj);
            s3 = computeThirdClosedMoment(obj);
            s4 = computeFourthClosedMoment(obj);
            s5 = computeFifthClosedMoment(obj);

            obj.closedMoment = s1(1) + s2(1) + s3(1) + s4(1) + s5(1);
            obj.Ma = s1(2) + s2(2) + s3(2) + s4(2) + s5(2);

        end

        function s1 = computeFirstClosedMoment(obj)
            q = obj.fluxSection(1);
            xs = obj.xsDist;
            xc = obj.xcDist;
            h1 = obj.height(1);
            d = obj.airfolDist;
            t1 = obj.thickness(1);

            s1(1) = int(q*(xs-xc), [0, h1/2]);
            s1(2) = int(q*d*t1, [0 h1/2]);
        end

        function s2 = computeSecondClosedMoment(obj)
            syms y
            q = obj.fluxSection(2);
            h1 = obj.height(1);
            h2 = obj.height(2);
            d = obj.airfolDist;
            t3 = obj.thickness(3);

            s2(1) = int(q*(h1/2 - y*(h1-h2)/(2*d)), [0, d]);
            s2(2) = int(q*t3*(h1/2 - y*(h1-h2)/(2*d)),[0 d]);
        end

        function s3 = computeThirdClosedMoment(obj)
            q = obj.fluxSection(3);
            h2 = obj.height(2);
            xc = obj.xcDist;
            d = obj.airfolDist;

            s3(1) = int(q*(xc-d), [0, h2]);
            s3(2) = int(q*0, [h2/2  -h2/2]);
        end

        function s4 = computeFourthClosedMoment(obj)
            syms y
            q = obj.fluxSection(4);
            h1 = obj.height(1);
            h2 = obj.height(2);
            d = obj.airfolDist;
            t3 = obj.thickness(3);

            s4(1) = int(q*(-h2/2 - y*(h1-h2)/(2*d)), [0, d]);
            s4(2) = int(q*t3*(-h2/2 - y*(h1-h2)/(2*d)), [0 d]);
        end

        function s5 = computeFifthClosedMoment(obj)
            q = obj.fluxSection(5);
            h1 = obj.height(1);
            xs = obj.xsDist;
            xc = obj.xcDist;
            d = obj.airfolDist;
            t1 = obj.thickness(1);

            s5(1) = int(q*(xs-xc), [0, h1/2]);
            s5(2) = int(q*d*t1, [0 h1/2]);
        end

        function s = computeInteriorSection(obj)
            h1 = obj.height(1);
            h2 = obj.height(2);
            d = obj.airfolDist; 
            
            s = h2*d + d*(h1-h2)/2;
        end


    end

end