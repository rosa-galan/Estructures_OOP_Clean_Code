classdef SystemSolver < handle

    properties (Access = public)
        u
        R
        uL
    end

    properties (Access = private)
        solverType
        nodalConnectivities
        dimensions
        xDist
        matConnectivities
        matProp
        xaDist
        xcDist
        xmDist
        fdata
        lift
        weight
        fixNodes
    end

    methods (Access = public)

        function obj = SystemSolver(cParams)
            obj.init(cParams);
        end

        function compute(obj)
            obj.solveSystem();
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.nodalConnectivities = cParams.Tn;
            obj.dimensions          = cParams.dim;
            obj.xDist               = cParams.x;
            obj.matConnectivities   = cParams.Tm;
            obj.matProp             = cParams.matProp;
            obj.xaDist              = cParams.xa;
            obj.xcDist              = cParams.xc;
            obj.xmDist              = cParams.xm;
            obj.fdata               = cParams.fdata;
            obj.lift                = cParams.liftForceDist;
            obj.weight              = cParams.weightForceDist;
            obj.fixNodes            = cParams.fixNodes;
            obj.solverType          = cParams.type;
            
        end

        function solveSystem(obj)

            stiffMatrix = computeStiffnessMatrix(obj);
            KG = stiffMatrix.KG;

            extForce = computeExternalForce(obj);
            Fext = extForce.Fext;

            fNodes = computeFixedNodes(obj);
            ur = (fNodes.ur)';
            vr = (fNodes.vr)';
            vf = (fNodes.vf)';

            Kff   = KG(vf,vf);
            Kfr   = KG(vf,vr);
            Krf   = KG(vr,vf);
            Krr   = KG(vr,vr);
            Fextf = Fext(vf,1);
            Fextr = Fext(vr,1);

            LHS = Kff;
            RHS = Fextf - Kfr*ur;

            disp = Solver.create(LHS,RHS,obj.solverType);
            obj.uL = disp.uL;

            obj.R = Krr*ur + Krf*obj.uL - Fextr;
            obj.u(vf, 1) = obj.uL';
            obj.u(vr, 1) = ur;
        end

        function stiffMatrix = computeStiffnessMatrix(obj)

            s.Tn      = obj.nodalConnectivities;
            s.dim     = obj.dimensions;
            s.x       = obj.xDist;
            s.Tm      = obj.matConnectivities;
            s.matProp = obj.matProp;

            stiffMatrix = GlobalStiffnessMatrixComputer(s);
            stiffMatrix.compute();
        end

        function extForce = computeExternalForce(obj)

            s.dim             = obj.dimensions;
            s.x               = obj.xDist;
            s.Tn              = obj.nodalConnectivities;
            s.xa              = obj.xaDist;
            s.xc              = obj.xcDist;
            s.xm              = obj.xmDist;
            s.fdata           = obj.fdata;
            s.liftForceDist   = obj.lift;
            s.weightForceDist = obj.weight;

            extForce = GlobalExternalForceComputer(s);
            extForce.compute();
        end

        function fNodes = computeFixedNodes(obj)
            
            s.dim = obj.dimensions;
            s.fixNodes = obj.fixNodes;

            fNodes = FixedDOFComputer(s);
            fNodes.compute();
        end
        
    end
    
end