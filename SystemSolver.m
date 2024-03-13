classdef SystemSolver < handle

    properties (Access = public)
        u
        R
        uL
        LHS
        RHS
    end

    properties (Access = private)
        dim
        KG
        Fext
        ur
        vr
        vf
        Kff
        Kfr
        Krf
        Krr
        Fextf
        Fextr
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
            obj.dim = cParams.dim;
            obj.KG = cParams.KG;
            obj.Fext = cParams.Fext;
            obj.ur = (cParams.ur)';
            obj.vr = (cParams.vr)';
            obj.vf = (cParams.vf)';
            obj.Kff = obj.KG(obj.vf, obj.vf);
            obj.Kfr = obj.KG(obj.vf, obj.vr);
            obj.Krf = obj.KG(obj.vr, obj.vf);
            obj.Krr = obj.KG(obj.vr, obj.vr);
            obj.Fextf = obj.Fext(obj.vf, 1);
            obj.Fextr = obj.Fext(obj.vr, 1);
            obj.LHS = obj.Kff;
            obj.RHS = obj.Fextf - obj.Kfr*obj.ur;
            disp = Solver.create(obj.LHS,obj.RHS);
            obj.uL = disp.uL;
        end

        function solveSystem(obj)
            obj.R = obj.Krr*obj.ur + obj.Krf*obj.uL - obj.Fextr;
            obj.u(obj.vf, 1) = obj.uL';
            obj.u(obj.vr, 1) = obj.ur;
        end

    end
    
end