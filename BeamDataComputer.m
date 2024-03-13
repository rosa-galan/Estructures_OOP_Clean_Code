classdef BeamDataComputer < handle
  
    properties (Access = public)
        Nel 
        defHist
        torsHist
        materialProperties 
        fixNodes
        x 
        Tm 
        Tn 
        fData 
        dim 
    end

    methods (Access = public)

        function obj = BeamDataComputer(cParams)
            obj.init(cParams);
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.Nel = 32;
            obj.defHist = zeros(size(obj.Nel,1), obj.Nel(end));
            obj.torsHist = zeros(size(obj.Nel,1), obj.Nel(end));
            obj.materialProperties = [cParams.E cParams.Ixx*1e-12 cParams.G cParams.Jc ; 1 1 1 1];
            obj.fixNodes = [1 1 0; 1 2 0;1 3 0];
            obj.x = [0:cParams.b/obj.Nel:cParams.b]';  
            obj.Tn = zeros(size(obj.x,1)-1,2);  
            for j = 1:size(obj.x)-1
                obj.Tn(j,:) = [j j+1];
            end 
            obj.Tm = ones(size(obj.Tn,1),1);
            obj.fData = [obj.Nel/4+1 1 -cParams.Me*cParams.g; obj.Nel/4+1 3 -cParams.Me*cParams.g*(cParams.xe-cParams.xm)/1000]; 
            obj.dim.nd = size(obj.x,2);   
            obj.dim.nel = size(obj.Tn,1);
            obj.dim.nnod = size(obj.x,1); 
            obj.dim.nne = size(obj.Tn,2); 
            obj.dim.ni = 3;          
            obj.dim.ndof = obj.dim.nnod*obj.dim.ni;  
        end

    end

end