classdef ElementalStiffnessMatrixComputer < handle

    properties (Access = public)
        Kel
    end

    properties (Access = private)
        x
        Tm
        matProp
        Tn
        dimensions
        gProp
        KelBending
        KelTorsion
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
            obj.x = cParams.x;
            obj.Tm = cParams.Tm;
            obj.matProp = cParams.matProp;
            obj.Tn = cParams.Tn;
            obj.dimensions = cParams.dim;
            obj.gProp = computeGeneralProperties(obj);
            obj.KelBending = computeKelBending(obj);
            obj.KelTorsion = computeKelTorsion(obj);
        end

        function computeElementalStiffnessMatrix(obj)
            dim = obj.dimensions;
            obj.Kel = zeros(dim.ni*dim.nne, dim.ni*dim.nne, dim.nel);
            for e=1:dim.nel
                obj.Kel(:,:,e) = obj.KelBending + obj.KelTorsion;
            end
        end

        function s = computeGeneralProperties(obj)
            for e=1:obj.dimensions.nel
                s.le = abs(obj.x(obj.Tn(e,2),1)- obj.x(obj.Tn(e,1),1))/1000;
                s.matEl = obj.matProp(obj.Tm(e),:);
                s.Ee = s.matEl(1);
                s.Ie = s.matEl(2);
                s.Ge = s.matEl(3);
                s.Je = s.matEl(4);
            end
        end

        function s = computeKelBending(obj)
            Ee = obj.gProp.Ee;
            Ie = obj.gProp.Ie;
            le = obj.gProp.le;
            for e=1:obj.dimensions.nel
                s = (Ee*Ie/(le)^3) * [ 12     6*le     0    -12     6*le      0;
                                      6*le  4*(le)^2   0   -6*le   2*(le)^2   0;
                                       0      0        0      0       0       0;
                                      -12   -6*le      0     12     -6*le     0;
                                      6*le  2*(le)^2   0    -6*le   4*(le)^2  0;
                                        0      0       0      0       0       0;];
            end
        end

        function s = computeKelTorsion(obj)
             for e=1:obj.dimensions.nel
                 Ge = obj.gProp.Ge;
                 Je = obj.gProp.Je;
                 le = obj.gProp.le;
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