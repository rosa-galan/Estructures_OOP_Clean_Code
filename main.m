classdef main < handle

    properties (Access = public)
        internalForces
        externalForces
        deflexion
        StiffnessMatrix
        tanStress
        normStress
        displacement
    end

    properties (Access = private)
        problemData
        beamData
        forceZ
        liftForce
        liftMoment
        weightMoment
        weightForce
        shearCenter
        Ma
        solverType
    end

    methods (Access = public)

        function obj = main(type)
            obj.init(type);
        end

        function compute(obj)
            obj.computeIntermediateValues();
            obj.computeStress();
            obj.computeBeamData();
            obj.solveBeam();
        end

    end
   
    methods(Access = private)

        function init(obj,type)
            obj.problemData = ConstantsComputer();
            obj.solverType = type;
        end

        function computeIntermediateValues(obj)
            liftAndWeight    = computeLiftAndWeight(obj);
            obj.liftForce    = liftAndWeight.liftForceDist;
            obj.liftMoment   = liftAndWeight.liftMomentDist;
            obj.weightForce  = liftAndWeight.weightForceDist;
            obj.weightMoment = liftAndWeight.weightMomentDist;

            Fz = computeVerticalForce(obj);
            obj.forceZ = Fz.Sz;

            fluxClosedSection = computeFluxClosedSection(obj);
            obj.Ma = fluxClosedSection.Ma;

            xsc = computeShearCenter(obj);
            obj.shearCenter = xsc.shearCenter;
        end

        function computeStress(obj)
            cParams = obj.problemData;
            s.be = cParams.be;
            s.t1 = cParams.t1;
            s.t2 = cParams.t2;
            s.t3 = cParams.t3;
            s.h1 = cParams.h1;
            s.h2 = cParams.h2;
            s.d  = cParams.d;
            s.xs = cParams.xs;
            s.xc = cParams.xc;
            s.a  = cParams.a;
            s.b  = cParams.b;
            s.We = cParams.We;
            s.g  = cParams.g;
            s.Ma = obj.Ma;
            s.liftMomentDist = obj.liftMoment;
            s.weightMomentDist = obj.weightMoment;
            s.Sz = obj.forceZ;

            stress = StressComputer(s);
            stress.compute();

            obj.tanStress = stress.tanStress;
        end

        function computeBeamData(obj)
            data = obj.problemData;
            bData = BeamDataComputer(data);
            bData.compute();
            obj.beamData = bData;         
        end

        function solveBeam(obj)
            globalStiffMatrix = computeGlobalStiffMatrix(obj);
            obj.StiffnessMatrix = globalStiffMatrix.KG;
            
            Fext = computeExternalForce(obj);
            obj.externalForces = Fext.Fext;

            displ = computeDisplacements(obj);
            obj.displacement = displ.u;

            Fint = computeInternalForce(obj);
            obj.internalForces = Fint;

            defl = computeDeflexion(obj);
            obj.deflexion = defl;    
        end

        function forces = computeLiftAndWeight(obj)
            cParams     = obj.problemData;
            s.rho    = cParams.rho;
            s.V      = cParams.V;
            s.Cl     = cParams.Cl;
            s.lambda = cParams.lambda;
            s.g      = cParams.g;
            s.b      = cParams.b;
            s.c      = cParams.c;

            forces = LiftWeightComputer(s);
            forces.compute();
        end

        function forceZ = computeVerticalForce(obj)
            cParams     = obj.problemData;
            s.g      = cParams.g;
            s.b      = cParams.b;
            s.c      = cParams.c;
            s.We     = cParams.We;
            s.rho    = cParams.rho;
            s.V      = cParams.V;
            s.Cl     = cParams.Cl;
            s.lambda = cParams.lambda;

            forceZ = ForceZComputer(s);
            forceZ.compute();
        end

        function fluxClosedSection = computeFluxClosedSection(obj)
            cParams = obj.problemData;
            s.Sz = obj.forceZ;
            s.be = cParams.be;
            s.t1 = cParams.t1;
            s.t2 = cParams.t2;
            s.t3 = cParams.t3;
            s.h1 = cParams.h1;
            s.h2 = cParams.h2;
            s.d  = cParams.d;
            s.xs = cParams.xs;
            s.xc = cParams.xc;
            s.a  = cParams.a;

            fluxClosedSection = ClosedSectionFluxComputer(s);
            fluxClosedSection.compute();
        end

        function shearCenter = computeShearCenter(obj)
            cParams = obj.problemData;
            s.Sz = obj.forceZ;
            s.t1 = cParams.t1;
            s.t2 = cParams.t2;
            s.t3 = cParams.t3;
            s.h1 = cParams.h1;
            s.h2 = cParams.h2;
            s.d  = cParams.d;
            s.xs = cParams.xs;
            s.xc = cParams.xc;
            s.a  = cParams.a;

            shearCenter = ShearCenterComputer(s);
            shearCenter.compute();
        end

        function globalStiffMatrix = computeGlobalStiffMatrix(obj)
            cParams   = obj.beamData;
            s.x       = cParams.x;
            s.Tm      = cParams.Tm;
            s.matProp = cParams.materialProperties;
            s.Tn      = cParams.Tn;
            s.dim     = cParams.dim;

            globalStiffMatrix = GlobalStiffnessMatrixComputer(s);
            globalStiffMatrix.compute();
        end

        function Fext = computeExternalForce(obj)
            cParams           = obj.beamData;
            data              = obj.problemData;
            s.dim             = cParams.dim;
            s.x               = cParams.x;
            s.Tn              = cParams.Tn;
            s.fdata           = cParams.fData;
            s.xa              = data.xa;
            s.xc              = data.xc;
            s.xm              = data.xm;
            s.liftForceDist   = obj.liftForce;
            s.weightForceDist = obj.weightForce;

            Fext = GlobalExternalForceComputer(s);
            Fext.compute();
        end

        function u = computeDisplacements(obj)
            cParams           = obj.beamData;
            data              = obj.problemData;
            s.dim             = cParams.dim;
            s.x               = cParams.x;
            s.Tn              = cParams.Tn;
            s.Tm              = cParams.Tm;
            s.matProp         = cParams.materialProperties;
            s.fdata           = cParams.fData;
            s.xa              = data.xa;
            s.xc              = data.xc;
            s.xm              = data.xm;
            s.type            = obj.solverType;
            s.liftForceDist   = obj.liftForce;
            s.weightForceDist = obj.weightForce;
            s.fixNodes        = cParams.fixNodes;

            u = SystemSolver(s);
            u.compute();
        end

        function Fint = computeInternalForce(obj)
            cParams   = obj.beamData;
            s.dim     = cParams.dim;
            s.x       = cParams.x;
            s.Tn      = cParams.Tn;
            s.Tm      = cParams.Tm;
            s.matProp = cParams.materialProperties;
            s.u       = obj.displacement;

            Fint = InternalForcesComputer(s);
            Fint.compute();
        end

        function defl = computeDeflexion(obj)
            cParams = obj.beamData;
            s.Nel   = cParams.nElements;

            defl = DeflexionComputer(s);
            defl.compute();
        end

    end

end