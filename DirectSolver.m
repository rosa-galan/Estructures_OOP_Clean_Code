classdef DirectSolver < Solver

    methods (Access = public)
            
         function obj = DirectSolver(L,R)
            obj.init(L,R);
         end

         function computeDisplacement(obj)
             obj.uL = inv(obj.LHS)*obj.RHS;       
         end

    end

end