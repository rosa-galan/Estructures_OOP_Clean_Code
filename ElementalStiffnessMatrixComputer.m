classdef ElementalStiffnessMatrixComputer < handle

    properties (Access = public)
        Kel
    end

    properties (Access = private)
        xDist
        matConnectivities
        matProp
        nodalConnectivities
        dimensions
        gProp
    end

    methods (Access = public)

        function obj = ElementalStiffnessMatrixComputer(cParams)
            obj.init(cParams);
        end

        function compute(obj)
            obj.computeElementalStiffnessMatrix();
        end
        
    end

    methods (Access = private)

        function init(obj,cParams)
            obj.xDist               = cParams.x;
            obj.matConnectivities   = cParams.Tm;
            obj.matProp             = cParams.matProp;
            obj.nodalConnectivities = cParams.Tn;
            obj.dimensions          = cParams.dim;
            obj.gProp               = computeGeneralProperties(obj);
        end

        function computeElementalStiffnessMatrix(obj)
            KelBending = computeKelBending(obj);
            KelTorsion = computeKelTorsion(obj);
            dim        = obj.dimensions;
            
            obj.Kel = zeros(dim.ni*dim.nne, dim.ni*dim.nne, dim.nel);
            for e=1:dim.nel
                obj.Kel(:,:,e) = KelBending + KelTorsion;
            end
        end

        function s = computeGeneralProperties(obj)
            dim = obj.dimensions;
            Tn  = obj.nodalConnectivities;
            Tm  = obj.matConnectivities;
            x   = obj.xDist;
            m   = obj.matProp;

            for e=1:dim.nel
                s.le    = abs(x(Tn(e,2),1)- x(Tn(e,1),1))/1000;
                s.matEl = m(Tm(e),:);
                s.Ee    = s.matEl(1);
                s.Ie    = s.matEl(2);
                s.Ge    = s.matEl(3);
                s.Je    = s.matEl(4);
            end
        end

        function s = computeKelBending(obj)
            Ee  = obj.gProp.Ee;
            Ie  = obj.gProp.Ie;
            le  = obj.gProp.le;
            dim = obj.dimensions;

            for e=1:dim.nel
                s = (Ee*Ie/(le)^3) * [ 12     6*le     0    -12     6*le      0;
                                      6*le  4*(le)^2   0   -6*le   2*(le)^2   0;
                                       0      0        0      0       0       0;
                                      -12   -6*le      0     12     -6*le     0;
                                      6*le  2*(le)^2   0    -6*le   4*(le)^2  0;
                                        0      0       0      0       0       0;];
            end
        end

        function s = computeKelTorsion(obj)
            Ge  = obj.gProp.Ge;
            Je  = obj.gProp.Je;
            le  = obj.gProp.le;
            dim = obj.dimensions;

             for e=1:dim.nel
                   s = (Ge*Je/le) * [0  0   0   0   0   0;
                                     0  0   0   0   0   0;
                                     0  0   1   0   0   -1;
                                     0  0   0   0   0   0;
                                     0  0   0   0   0   0;
                                     0  0   -1  0   0   1;];
             end
        end
        
    end

end