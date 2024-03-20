classdef ConstantsComputer
    
    properties (Access = public)
        t1 
        t2 
        t3 
        c 
        h1
        h2 
        xs 
        d 
        a 
        b 
        rho 
        V 
        Cl
        lambda 
        g 
        We 
        be
        E 
        G 
        Jo 
        Jc 
        xa 
        xm 
        xe 
        xc
    end

    methods (Access = public)

        function obj = ConstantsComputer()
            obj.t1 = 25; 
            obj.t2 = 15; 
            obj.t3 = 5; 
            obj.c = 2000; 
            obj.b = 16e3; 
            obj.rho = 1.225; 
            obj.V = 800/3.6; 
            obj.Cl = 0.1;
            obj.lambda = 150; 
            obj.g = 9.81; 
            obj.We = 2250; 
            obj.E = 200e9; 
            obj.G = 75e9; 
            obj.Jo = 2.446e6*1e-12; 
            obj.Jc = 6.372e8*1e-12; 
            obj.h1 = 0.2 * obj.c;
            obj.h2 = 0.15 * obj.c;
            obj.xs = 0.3 * obj.c;
            obj.d = 0.3 * obj.c;
            obj.a = sqrt(((obj.h1 - obj.h2) / 2)^2 + obj.d^2);
            obj.be = 0.25 * obj.b;
            obj.xa = 0.25 * obj.c;
            obj.xm = 0.45 * obj.c; 
            obj.xe = 0.3 * obj.c;
            obj.xc = (obj.h1*obj.t1*obj.xs + obj.h2*obj.t2*(obj.d+obj.xs)+ 2*obj.t3*obj.a*(obj.xs+obj.d/2)) / (obj.h1*obj.t1 + obj.h2*obj.t2 + 2*obj.a*obj.t3);
        end
       
    end
    
end
