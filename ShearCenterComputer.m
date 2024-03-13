classdef ShearCenterComputer < handle

    properties (Access = public)
        shearCenter
    end

    properties (Access = private)
        qSection
        Sz
        Mc
        t1 
        t2 
        t3 
        h1 
        h2 
        d 
        xs 
        xc
       distCentroide
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
            obj.t1 = cParams.t1;
            obj.t2 = cParams.t2;
            obj.t3 = cParams.t3;
            obj.h1 = cParams.h1;
            obj.h2 = cParams.h2;
            obj.d = cParams.d;
            obj.xc = cParams.xc;
            obj.xs = cParams.xs;
            obj.qSection = cParams.qSection;
            obj.Sz = cParams.Sz;
            obj.Mc = computeCentroideMoments(obj);
            obj.distCentroide = double(obj.Mc) / obj.Sz;
        end

        function computeShearCenter(obj)
            obj.shearCenter = double(obj.xc + obj.distCentroide);
        end

        function s = computeCentroideMoments(obj)
            syms y
            s = int(obj.qSection(1) .*obj.t1.* (obj.xs-obj.xc), [0.0, obj.h1/2]) + ...
            int(obj.qSection(2) .*obj.t3.* (obj.h1/2 - y*(obj.h1-obj.h2)/(2*obj.d)) ,[0, obj.d]) + ...
            int(obj.qSection(3) .*obj.t2.* (obj.xc-obj.d), [0 , obj.h2]) + ... 
            int(obj.qSection(4) .*obj.t3.* (-obj.h2/2 - y*(obj.h1-obj.h2)/(2*obj.d)), [0 , obj.d]) + ...
            int(obj.qSection(5) .*obj.t1.* (obj.xs-obj.xc), [0 , obj.h1/2]);
        end

    end
    
end