classdef DivisionAlgebra < replab.Str
% Division algebra over the real numbers
    properties (SetAccess = protected)
        description;
        
        % Description of the division algebra, see
        % E. de Klerk, C. Dobre, and D. V. Ṗasechnik, Math. Program. 129, 91 (2011).
        d; % Dimension of the algebra
           
        % Let now D(i) be a basis of the division algebra
        % such that every element is written sum_i w(i) D(i)
        gamma; % Multiplication parameters gamma(d,d,d) such that
               % D(i)*D(j) = sum_k gamma(i,j,k) D(k)
        
        % Matrix realization with (real) matrices of size m x m
        % such that M = sum_i w(i) basis(i)
        % and trace(basis(:,:,i)*dualBasis(:,:,j)) = delta(i,j)
        m;
        basis; % basis(m,m,d)
        dualBasis; % basis(m,m,d)
    end 
    
    methods
        
        function self = DivisionAlgebra(description, gamma, basis, dualBasis)
            self.description = description;
            self.d = size(gamma, 1);
            self.gamma = gamma;
            self.n = size(basis, 1);
            self.basis = basis;
            self.dualBasis = dualBasis;
        end
        
    end
   
    methods (Static)
        
        function D = real
            D = DivisionAlgebra('Real division algebra', 1, 1, 1);
        end
        
        function D = complex
            gamma = zeros(2,2,2);
            gamma(:,:,1) = [1  0
                            0 -1];
            gamma(:,:,2) = [0  1
                            1  0];
            basis = zeros(2,2,2);
            dualBasis = zeros(2,2,2);
            basis(:,:,1) = [1 0
                            0 1];
            basis(:,:,2) = [0 -1
                            1  0];
            dualBasis(:,:,1) = [1 0
                                0 1]/2;
            dualBasis(:,:,2) = [0 -1
                                1  0]/2;
            D = DivisionAlgebra('Complex division algebra', gamma, basis, dualBasis)
  
end
