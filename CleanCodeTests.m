classdef CleanCodeTests < matlab.unittest.TestCase

     methods (Test)

         function verifyTangentialStress(testCase)
             load('OriginalData.mat','tau');
             originalSolution = double(tau);
             type = 'Iterative';
             stress = main(type);
             stress.compute();
             actualSolution = stress.tanStress;
             testCase.verifyEqual(actualSolution,originalSolution,"AbsTol",1e-10);
         end

         function verifyStiffnessMatrix(testCase)
             load('OriginalData.mat','KG_original');
             originalSolution = KG_original;
             type = 'Iterative';
             stiffMatrix = main(type);
             stiffMatrix.compute();
             actualSolution = stiffMatrix.StiffnessMatrix;
             testCase.verifyEqual(actualSolution,originalSolution,"AbsTol",1e-10);
         end
        
         function verifyIterativeCase(testCase)
             load('IterativeResults.mat','u','Fint');
             originalDisp = u;
             originalQ = Fint.Q;
             originalMb = Fint.Mb;
             originalMt = Fint.Mt;
             type = 'Iterative';
             solution = main(type);
             solution.compute();
             actualDisp = solution.displacement;
             actualQ = solution.internalForces.Q;
             actualMb = solution.internalForces.Mb;
             actualMt = solution.internalForces.Mt;
             testCase.verifyEqual(actualDisp,originalDisp,"AbsTol",1e-10);
             testCase.verifyEqual(actualQ,originalQ,"AbsTol",1e-10);
             testCase.verifyEqual(actualMb,originalMb,"AbsTol",1e-10);
             testCase.verifyEqual(actualMt,originalMt,"AbsTol",1e-10);
         end

         function verifyDirectCase(testCase)
             load('DirectResults.mat','u','Fint');
             originalDisp = u;
             originalQ = Fint.Q;
             originalMb = Fint.Mb;
             originalMt = Fint.Mt;
             type = 'Direct';
             solution = main(type);
             solution.compute();
             actualDisp = solution.displacement;
             actualQ = solution.internalForces.Q;
             actualMb = solution.internalForces.Mb;
             actualMt = solution.internalForces.Mt;
             testCase.verifyEqual(actualDisp,originalDisp,"AbsTol",1e-10);
             testCase.verifyEqual(actualQ,originalQ,"AbsTol",1e-10);
             testCase.verifyEqual(actualMb,originalMb,"AbsTol",1e-10);
             testCase.verifyEqual(actualMt,originalMt,"AbsTol",1e-10);
         end

         function verifyDeflexion(testCase)
             load('OriginalData.mat','Deflexion');
             originalDefX = Deflexion.x;
             originalDefY = Deflexion.y;
             originalDefTheta = Deflexion.theta;
             type = 'Iterative';
             solution = main(type);
             solution.compute();
             actualDefX = solution.deflexion.defX;
             actualDefY = solution.deflexion.defY;
             actualDefTheta = solution.deflexion.defTheta;
             testCase.verifyEqual(actualDefX,originalDefX,"AbsTol",1e-10);
             testCase.verifyEqual(actualDefY,originalDefY,"AbsTol",1e-10);
             testCase.verifyEqual(actualDefTheta,originalDefTheta,"AbsTol",1e-10);
         end
       
     end

end