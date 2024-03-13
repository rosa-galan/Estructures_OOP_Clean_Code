classdef Solver < handle

    properties (Access = public)
          uL
    end

    properties (Access = protected)
        LHS
        RHS
    end

    methods (Access = public, Static)

        function obj = create(L,R)
            type = input('Enter a number - 1 [Iterative Solver] or 2 [Direct Solver]: ');
            switch type
                 case 1
                    obj = IterativeSolver(L,R);
                    obj.computeDisplacement();
                 case 2
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