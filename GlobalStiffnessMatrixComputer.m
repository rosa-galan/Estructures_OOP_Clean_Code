classdef GlobalStiffnessMatrixComputer < handle

      properties (Access = public)
          KG
      end

      properties (Access = private)
          nodalConnectivities
          dimensions
          xDist
          matConnectivities
          matProp
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
              obj.nodalConnectivities = cParams.Tn;
              obj.dimensions          = cParams.dim;
              obj.xDist               = cParams.x;
              obj.matConnectivities   = cParams.Tm;
              obj.matProp             = cParams.matProp;
          end

          function computeGlobalStiffnessMatrix(obj)

              connectMatrix = computeConnectMatrix(obj);     
              Td = connectMatrix.Td;

              elemStiffMatrix = computeElemMatrix(obj);
              Kel = elemStiffMatrix.Kel;

              dim = obj.dimensions;
              nneNi = dim.ni*dim.nne;

              obj.KG = zeros(dim.nnod*dim.ni, dim.nnod*dim.ni);
              for e=1:dim.nel
                for i=1:nneNi
                    I=Td(e,i);
                    for j=1:nneNi
                        J=Td(e, j);
                        obj.KG(I,J)=obj.KG(I,J)+Kel(i,j,e);
                    end
                end
              end
          end

          function connectMatrix = computeConnectMatrix(obj)
              s.dim = obj.dimensions;
              s.Tn  = obj.nodalConnectivities;

              connectMatrix = ConnectDOFComputer(s);
              connectMatrix.compute();
          end

          function elemStiffMatrix = computeElemMatrix(obj)

              s.dim     = obj.dimensions;
              s.Tn      = obj.nodalConnectivities;
              s.x       = obj.xDist;
              s.Tm      = obj.matConnectivities;
              s.matProp = obj.matProp;

              elemStiffMatrix = ElementalStiffnessMatrixComputer(s);
              elemStiffMatrix.compute()
          end

      end

end