classdef ShearCenterComputer < handle

    properties (Access = public)
        shearCenter
    end

    properties (Access = private)
        qSection
        forceZ
        Mc
        thickness 
        height 
        airfolDist 
        xsDist 
        xcDist
        centroideDist
        width
        fluxSection
        inertiaX
    end

    methods (Access = public)
        
        function obj = ShearCenterComputer(cParams)
            obj.init(cParams);
        end

        function compute(obj)
            obj.computeShearCenter();
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
            obj.width = cParams.a;
            obj.xcDist = cParams.xc;
            obj.xsDist = cParams.xs;
            obj.forceZ = cParams.Sz;
            
            flux = computeFlux(obj);
            obj.fluxSection = flux.qSection;
            
            inertias     = computeInertias(obj);
            obj.inertiaX = inertias.IxxT;
            
            obj.Mc = computeCentroideMoments(obj);
            obj.centroideDist = double(obj.Mc) / obj.forceZ;
        end

        function computeShearCenter(obj)
            obj.shearCenter = double(obj.xcDist + obj.centroideDist);
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

        function s = computeCentroideMoments(obj)
            syms y
            s1 = computeFirstSection(obj);
            s2 = computeSecondSection(obj);
            s3 = computeThirdSection(obj);
            s4 = computeFourthSection(obj);
            s5 = computeFifthSection(obj);

            s = s1 + s2 + s3 + s4 + s5;
        end

        function s1 = computeFirstSection(obj)
            t1 = obj.thickness(1);
            xs = obj.xcDist;
            xc = obj.xcDist;
            h1 = obj.height(1);
            q  = obj.fluxSection(1);

            s1 = int(q.*t1.*(xs-xc), [0, h1/2]);
        end

        function s2 = computeSecondSection(obj)
            syms y
            t3 = obj.thickness(3);
            d  = obj.airfolDist;
            h1 = obj.height(1);
            h2 = obj.height(2);
            q  = obj.fluxSection(2);

            s2 = int(q.*t3.*(h1/2 - y*(h1-h2)/(2*d)), [0, d]);
        end

        function s3 = computeThirdSection(obj)
            t2 = obj.thickness(2);
            xc = obj.xcDist;
            d  = obj.airfolDist;
            h2 = obj.height(2);
            q  = obj.fluxSection(3);

            s3 = int(q.*t2.*(xc-d), [0 , h2]);
        end

         function s4 = computeFourthSection(obj)
            syms y 
            t3 = obj.thickness(3);
            d  = obj.airfolDist;
            h1 = obj.height(1);
            h2 = obj.height(2);
            q  = obj.fluxSection(4);

            s4 = int(q.*t3.* (-h2/2 - y*(h1-h2)/(2*d)), [0 , d]);
         end

         function s5 = computeFifthSection(obj)
            t1 = obj.thickness(1);
            xs = obj.xcDist;
            xc = obj.xcDist;
            h1 = obj.height(1);
            q  = obj.fluxSection(5);

            s5 = int(q.*t1.*(xs-xc), [0, h1/2]);
         end


    end
    
end