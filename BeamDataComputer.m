classdef BeamDataComputer < handle
  
    properties (Access = public)
        nElements
        materialProperties 
        fixNodes
        x 
        Tm 
        Tn 
        fData 
        dim
       
    end

    properties (Access = private)
        thickness
        height
        width
        airfolDist
        xcDist
        xsDist
        xeDist
        xmDist
        youngModulus
        closedTorsInertia
        shearModulus
        inertiaX
        engineMass
        gravity
    end

    methods (Access = public)

        function obj = BeamDataComputer(cParams)
            obj.init(cParams);
        end

        function compute(obj)
            obj.computeMaterialProperties();
            obj.computeConnectivites();
            obj.computePointLoads();
            obj.computeDimensions();
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.nElements = 32;
            obj.thickness(1)      = cParams.t1;
            obj.thickness(2)      = cParams.t2;
            obj.thickness(3)      = cParams.t3;
            obj.height(1)         = cParams.h1;
            obj.height(2)         = cParams.h2;
            obj.width             = cParams.a;
            obj.airfolDist        = cParams.d;
            obj.xcDist            = cParams.xc;
            obj.xsDist            = cParams.xs;
            obj.xeDist            = cParams.xe;
            obj.xmDist            = cParams.xm;
            obj.youngModulus      = cParams.E;
            obj.shearModulus      = cParams.G;
            obj.closedTorsInertia = cParams.Jc;
            obj.engineMass = cParams.We;
            obj.gravity = cParams.g;

            
            inertias = computeInertias(obj);
            obj.inertiaX = inertias.IxxT;
            
            obj.fixNodes = [1 1 0; 1 2 0; 1 3 0];
            obj.x = [0:cParams.b/obj.nElements:cParams.b]';  
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

         function computeMaterialProperties(obj)
             E   = obj.youngModulus;
             Ixx = obj.inertiaX;
             G   = obj.shearModulus;
             Jc  = obj.closedTorsInertia;

             obj.materialProperties = [E   Ixx*1e-12   G   Jc ;...
                                       1       1       1   1];
         end

         function computeConnectivites(obj)
            obj.Tn = zeros(size(obj.x,1)-1,2);  
            for j = 1:size(obj.x)-1
                obj.Tn(j,:) = [j j+1];
            end 

            obj.Tm = ones(size(obj.Tn,1),1);
         end

         function computePointLoads(obj)
             Nel = obj.nElements;
             Me  = obj.engineMass;
             g   = obj.gravity;
             xe  = obj.xeDist;
             xm  = obj.xmDist;

             obj.fData = [Nel/4+1 1 -Me*g; Nel/4+1 3 -Me*g*(xe-xm)/1000]; 
         end

         function computeDimensions(obj)
            obj.dim.nd   = size(obj.x,2);   
            obj.dim.nel  = size(obj.Tn,1);
            obj.dim.nnod = size(obj.x,1); 
            obj.dim.nne  = size(obj.Tn,2); 
            obj.dim.ni   = 3;          
            obj.dim.ndof = obj.dim.nnod*obj.dim.ni;  
         end


    end

end