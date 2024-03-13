classdef CleanCodeTests < matlab.unittest.TestCase

     methods (Test)
        
        % 1. TESTS TO VERIFY INERTIAS

        % 1.1. Test to verify Ixx
        function verifyInertiaX(testCase)
            load('OriginalData.mat','Ixx');
            originalSolution = Ixx;
            inertias = main();
            actualSolution = inertias.Ixx;
            testCase.verifyEqual(actualSolution,originalSolution);
        end

        % 1.2. Test to verify Izz
        function verifyInertiaZ(testCase)
            load('OriginalData.mat','Izz');
            originalSolution = Izz;
            inertias = main();
            actualSolution = inertias.Izz;
            testCase.verifyEqual(actualSolution,originalSolution);
        end

        % 2. TESTS TO VERIFY LIFT COMPUTATIONS
        
        % 2.1. Test to verify Lift Force Distribution
        function verifyLiftForceDist(testCase)
            load('OriginalData.mat','lift_dist');
            originalSolution = lift_dist;
            liftInfo = main();
            actualSolution = liftInfo.liftForceDist;
            testCase.verifyEqual(actualSolution,originalSolution);
        end 

        % 2.2. Test to verify Lift Moment Distribution
        function verifyLiftMomentDist(testCase)
            load('OriginalData.mat','liftM_dist');
            load('ParametersData.mat','s');
            cParams = s;
            originalSolution = liftM_dist;
            liftInfo = LiftWeightComputer(cParams);
            actualSolution = liftInfo.liftMomentDist;
            testCase.verifyEqual(actualSolution,originalSolution);
        end 

        % 3. TESTS TO VERIFY WEIGHT COMPUTATIONS
        
        % 3.1. Test to verify Weight Force Distribution
        function verifyWeightForceDist(testCase)
            load('OriginalData.mat','weight_dist');
            load('ParametersData.mat','s');
            cParams = s;
            originalSolution = weight_dist;
            weightInfo = LiftWeightComputer(cParams);
            actualSolution = weightInfo.weightForceDist;
            testCase.verifyEqual(actualSolution,originalSolution);
        end 

        % 3.2. Test to verify Weight Moment Distribution
        function verifyWeightMomentDist(testCase)
            load('OriginalData.mat','weightM_dist');
            load('ParametersData.mat','s');
            cParams = s;
            originalSolution = weightM_dist;
            weightInfo = LiftWeightComputer(cParams);
            actualSolution = weightInfo.weightMomentDist;
            testCase.verifyEqual(actualSolution,originalSolution);
        end 

        % 4. TEST TO VERIFY FORCE COMPUTATION

        % 4.1. Test to verify vertical force 
        function verifyForce(testCase)
            load('OriginalData.mat','Sz');
            load('ParametersData.mat','s');
            cParams = s;
            originalSolution = double(Sz);
            force = ForceZComputer(cParams);
            actualSolution = force.Sz;
            testCase.verifyEqual(actualSolution,originalSolution);
        end 

        % 5. TEST TO VERIFY FLUX AND SHEAR CENTER COMPUTATION

        % 5.1. Test to verify flux in open section 
        function verifyFlux(testCase)
            load('OriginalData.mat','flux');
            load('ParametersData.mat','s');
            cParams = s;
            originalSolution = double(flux);
            q = FluxComputer(cParams);
            actualSolution = q.q;
            testCase.verifyEqual(actualSolution,originalSolution, "AbsTol",1e-10);
        end 
        
        % 5.2. Test to verify shear center position
        function verifyShearCenter(testCase)
            load('OriginalData.mat','xsc');
            load('ParametersData.mat','s');
            cParams = s;
            originalSolution = xsc;
            xsc = ShearCenterComputer(cParams);
            actualSolution = xsc.shearCenter;
            testCase.verifyEqual(actualSolution,originalSolution, "AbsTol",1e-10);
        end 

        % 6. BEAM SOLVER TESTS

        % 6.1. Test to verify connectivity matrix
        function verifyConnectivityMatrix(testCase)
            load('OriginalData.mat','Td');
            load('ParametersData.mat','s');
            cParams = s;
            originalSolution = Td;
            TdMatrix = ConnectDOFComputer(cParams);
            actualSolution = TdMatrix.Td;
            testCase.verifyEqual(actualSolution,originalSolution, "AbsTol",1e-10);
        end 

        % 6.2. Test to verify elemental stiffness matrix
        function verifyStiffnessElementalMatrix(testCase)
            load('OriginalData.mat','Kel');
            load('ParametersData.mat','s');
            cParams = s;
            originalSolution = Kel;
            KelMatrix = ElementalStiffnessMatrixComputer(cParams);
            actualSolution = KelMatrix.Kel;
            testCase.verifyEqual(actualSolution,originalSolution, "AbsTol",1e-10);
        end 

        % 6.3. Test to verify global stiffness matrix
        function verifyStiffnessGlobalMatrix(testCase)
            load('OriginalData.mat','KG_original');
            load('ParametersData.mat','s');
            cParams = s;
            originalSolution = KG_original;
            KGMatrix = GlobalStiffnessMatrixComputer(cParams);
            actualSolution = KGMatrix.KG;
            testCase.verifyEqual(actualSolution,originalSolution, "AbsTol",1e-10);
        end 

         % 6.4. Test to verify elemental force vector
         function verifyElementalForce(testCase)
            load('OriginalData.mat','fel');
            load('ParametersData.mat','s');
            cParams = s;
            originalSolution = fel;
            elemForce = ElementalForceComputer(cParams);
            actualSolution = elemForce.fel;
            testCase.verifyEqual(actualSolution,originalSolution, "AbsTol",1e-10);
        end 

        % 6.5. Test to global force vector
        function verifyGlobalForce(testCase)
           load('OriginalData.mat','Fext');
           load('ParametersData.mat','s');
           cParams = s;
           originalSolution = Fext;
           extForce = GlobalExternalForceComputer(cParams);
           actualSolution = extForce.Fext;
           testCase.verifyEqual(actualSolution,originalSolution, "AbsTol",1e-10);
         end 

         % 6.6. Test to global force vector
         function verifyDOFArrays(testCase)
            load('OriginalData.mat','ur','vr','vf');
            load('ParametersData.mat','s');
            cParams = s;
            originalSolution{1} = ur;
            originalSolution{2} = vr;
            originalSolution{3} = vf;
            DOFarrays = FixedDOFComputer(cParams);
            actualSolution{1} = DOFarrays.ur;
            actualSolution{2} = DOFarrays.vr;
            actualSolution{3} = DOFarrays.vf;
            testCase.verifyEqual(actualSolution{1},originalSolution{1}, "AbsTol",1e-10);
            testCase.verifyEqual(actualSolution{2},originalSolution{2}, "AbsTol",1e-10);
            testCase.verifyEqual(actualSolution{3},originalSolution{3}, "AbsTol",1e-10);
         end 

        % % 6.7. Test to compare iterative and direct displacement
        % function compareIterativeDirectSolvers(testCase)
        %    load('IterativeDisplacement.mat','u');
        %    load('DirectDisplacement.mat','u');
        %    cParams = s;
        %    iterativeSolution = u;
        %    extForce = GlobalExternalForceComputer(cParams);
        %    actualSolution = extForce.Fext;
        %    testCase.verifyEqual(actualSolution,originalSolution, "AbsTol",1e-10);
        %  end 
       
     end

end