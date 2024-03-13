% Clean code main

function [s] =  main%(solverType) 

clear all
close all
clc
%% 1. Take all constant data from the ConstantsComputer class

data = ConstantsComputer();

data = data.computeDerivedProperties();

%% 2. Compute intertias

% First, select the constants needed to compute the inertia
s.t1 = data.t1;
s.t2 = data.t2;
s.t3 = data.t3;
s.h1 = data.h1;
s.h2 = data.h2;
s.a = data.a;
s.d = data.d;
s.xc = data.xc;
s.xs = data.xs;
s.xe = data.xe;
s.xm = data.xm;
s.xa = data.xa;

% Now call the class 
Inertias = InertiaComputer(s);
Inertias.compute();

% Save to s the values of interest for future computations
s.Ixx = Inertias.IxxT;
s.Izz = Inertias.IzzT;
s.J = Inertias.J;

%% 3. Compute lift and weight forces and moments

% First add the constants needed to compute 
s.rho = data.rho;
s.V = data.V;
s.Cl = data.Cl;
s.lambda = data.lambda;
s.g = data.g;
s.b = data.b;
s.c = data.c;
s.Me = data.Me;

% Call the class
LiftWeight = LiftWeightComputer(s);
LiftWeight.compute();

% Save the values
s.liftForceDist = LiftWeight.liftForceDist; 
s.liftIntegral = LiftWeight.liftIntegral;
s.liftMomentDist = LiftWeight.liftMomentDist;
s.weightForceDist = LiftWeight.weightForceDist;
s.weightMomentDist = LiftWeight.weightMomentDist;

%% 3. Compute z force Sz

% Call the class
Sz = ForceZComputer(s);
Sz.compute();

% Save the value for future computations
s.Sz = Sz.Sz;

%% 4. Compute flux in open section

% Call the class
FluxOpenSection = FluxComputer(s);
FluxOpenSection.compute();

% Save values of flux to s
s.q = FluxOpenSection.q;
s.qSection = FluxOpenSection.qSection;

%% 5. Compute shear center

% Call the class
shearCenter = ShearCenterComputer(s);
shearCenter.compute();

% Save the value
s.xsc = shearCenter.shearCenter;


%% 6. Compute flux in closed section

% Save necessary values from data
s.be = data.be;

% Call the class
FluxClosedSection = ClosedSectionFluxComputer(s);
FluxClosedSection.compute();

% Save necessary values

s.Ma = FluxClosedSection.Ma;

%% 7. Compute moment in the shear center

% 7.1 Open section case 
OpenMoment = OpenSectionMomentComputer(s);
OpenMoment.compute();
s.openMsc = OpenMoment.openMsc;

% 7.2. Closed section case
ClosedMoment = ClosedSectionMomentComputer(s);
ClosedMoment.compute();
s.closedMsc = ClosedMoment.closedMsc;


%% 8. Compute tangential and normal stress

% Call the class 
TanNormStress = StressComputer(s);
TanNormStress.compute();

%% 9. Take all the constant data from the beamDataComputer class

% Add the necessary parameters 
s.Jc = data.Jc;
s.G = data.G;
s.E = data.E;

% Call the class
beamData = BeamDataComputer(s);

% Save necessary values
s.Tn = beamData.Tn;
s.Tm = beamData.Tm;
s.dim = beamData.dim;
s.x = beamData.x;
s.matProp = beamData.materialProperties;
s.Tm = beamData.Tm;
s.fdata = beamData.fData;
s.fixNodes = beamData.fixNodes;
s.Nel = beamData.Nel;

%% 10. Beam solver

% 10.1. Create degrees of freedom connectivities matrix
ConnectMatrix = ConnectDOFComputer(s);
ConnectMatrix.compute();
s.Td = ConnectMatrix.Td;

% 10.2. Compute element stiffness matrices
StiffnessElementalMatrix = ElementalStiffnessMatrixComputer(s);
StiffnessElementalMatrix.compute();
s.Kel = StiffnessElementalMatrix.Kel;

% 10.3. Assemble global stiffness matrix
KG = GlobalStiffnessMatrixComputer(s);
KG.compute();
s.KG = KG.KG;

% 10.4. Compute element force vector
ElementalForce = ElementalForceComputer(s);
ElementalForce.compute();
s.fel = ElementalForce.fel;

% 10.5. Assemble global external forces vector
GlobalFext = GlobalExternalForceComputer(s);
GlobalFext.compute();
s.Fext = GlobalFext.Fext;

% 10.6. Create arrays of fixed and free degrees of freedom
fixDOF = FixedDOFComputer(s);
fixDOF.compute();
s.ur = fixDOF.ur;
s.vr = fixDOF.vr;
s.vf = fixDOF.vf;

%save('ParametersData.mat','s');

% 10.7. Solve system and obtain reactions and displacements

SystemSolution = SystemSolver(s);
SystemSolution.compute();
s.u = SystemSolution.u;
s.R = SystemSolution.R;

% 10.8. Compute internal forces distribution
internalForces = InternalForcesComputer(s);
internalForces.compute();
s.Q = internalForces.Q;
s.Mb = internalForces.Mb;
s.Mt = internalForces.Mt;

% 10.9 Deflexion angles
deflexion = DeflexionComputer(s);
deflexion.compute();
s.defX = deflexion.defX;
s.defY = deflexion.defY;
s.defTheta = deflexion.defTheta;


%% 11. Postprocess

% 11.1 Distribution plots

% Force and moment plots
figure(1);
subplot(3, 1, 1);
plot(s.x(s.Tn'),s.Q','color','r');
title("Shear force distribution (Q)");
xlabel('Position along the beam (mm)');
ylabel('Shear force (N)');
subplot(3, 1, 2);
plot(s.x(s.Tn'),s.Mb','color','g');
title("Bending moment distribution");
xlabel('Position along the beam (mm)');
ylabel('Bending moment (Nm)');
subplot(3, 1, 3);
plot(s.x(s.Tn'),s.Mt','color','b');
title("Torsion moment distribution")
xlabel('Position along the beam (mm)');
ylabel('Torsion moment (Nm)');

% Deflexion plots
figure(2);
subplot(3, 1, 1);
plot(s.x(1:s.Nel),s.u(s.defX,1),'color','b');
title("Vertical deflexion distribution");
xlabel('Position along the beam (mm)');
ylabel('Vertical deflexion (m)');
subplot(3, 1, 2);
plot(s.x(1:s.Nel),s.u(s.defY,1),'color','b');
title("Deflexion angle distribution - Bending");
xlabel('Position along the beam (mm)');
ylabel('Deflexion angle (º)');
subplot(3, 1, 3);
plot(s.x(1:s.Nel),s.u(s.defTheta,1),'color','b');
title("Deflexion angle distribution - Torsion");
xlabel('Position along the beam (mm)');
ylabel('Deflexion angle (º)');


% 11.2. Colormaps

% Open case
%sectionAnalysis('open',s.d/1000,s.h1/1000,s.h2/1000,s.t1/1000,s.t2/1000,s.t3/1000,s.E,s.G,-s.xsc/1000 + (s.xc-(s.xs+s.d/2))/1000,double(s.Sz),s.Mb(1),double(s.openMsc));

% Closed case
%sectionAnalysis('closed',s.d/1000,s.h1/1000,s.h2/1000,s.t1/1000,s.t2/1000,s.t3/1000,s.E,s.G,-(s.xc-(s.xs+s.d/2))/1000,double(s.Sz),-s.Mb(1),double(s.closedMsc));


end