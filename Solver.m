classdef Solver < handle

    properties (Access = public)
          uL
    end

    properties (Access = protected)
        LHS
        RHS
    end

    methods (Access = public, Static)

        function obj = create(L,R,type)
            switch type
                case 'Iterative'
                    obj = IterativeSolver(L,R);
                    obj.computeDisplacement();
                case 'Direct'
                    obj = DirectSolver(L,R);
                    obj.computeDisplacement();
                otherwise
                    error("The methodology chosen does not exist")
            end       
        end

    end
    
    methods (Access = protected)

        function init(obj,L,R)
            obj.LHS = L;
            obj.RHS = R; 
        end

    end

end