classdef GlobalStiffnessMatrixComputer < handle

      properties (Access = public)
          KG
      end

      properties (Access = private)
          Td
          dim
          Kel
          nneNi
      end

      methods (Access = public)

          function obj = GlobalStiffnessMatrixComputer(cParams)
              obj.init(cParams)
          end

          function compute(obj)
              obj.computeGlobalStiffnessMatrix();
          end

      end

      methods(Access = private)

          function init(obj, cParams)
              obj.Td = cParams.Td;
              obj.dim = cParams.dim;
              obj.Kel = cParams.Kel;
              obj.nneNi = obj.dim.ni*obj.dim.nne;
          end

          function computeGlobalStiffnessMatrix(obj)
              obj.KG = zeros(obj.dim.nnod*obj.dim.ni, obj.dim.nnod*obj.dim.ni);
              for e=1:obj.dim.nel
                for i=1:obj.nneNi
                    I=obj.Td(e,i);
                    for j=1:obj.nneNi
                        J=obj.Td(e, j);
                        obj.KG(I,J)=obj.KG(I,J)+obj.Kel(i,j,e);
                    end
                end
              end
          end

      end

end