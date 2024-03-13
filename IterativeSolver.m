classdef IterativeSolver < Solver

    methods (Access = public)

        function obj = IterativeSolver(L,R)
            obj.init(L,R);
        end
            
        function computeDisplacement(obj)
             obj.uL = pcg(obj.LHS,obj.RHS);       
        end  

    end

end