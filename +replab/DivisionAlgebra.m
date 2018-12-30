classdef DivisionAlgebra < replab.Monoid
% Division algebra over the reals, whose elements 
% are represented by d x 1 real column vectors
    properties (SetAccess = protected)
        shortName;
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
        
        function self = DivisionAlgebra(description, shortName, gamma, identity, basis, dualBasis)
            self.description = description;
            self.shortName = shortName;
            self.d = size(gamma, 1);
            self.gamma = gamma;
            self.identity = identity;
            self.m = size(basis, 1);
            self.basis = basis;
            self.dualBasis = dualBasis;
        end
        
        function b = eqv(self, x, y)
            b = isequal(x, y);
        end
        
        function s = sample(self)
            s = randi([-10 10], self.d, 1);
        end
        
        function X = projectMatrix(self, X)
            m = self.m;
            d = self.d;
            R = size(X, 1);
            C = size(X, 2);
            bR = R/m;
            bC = C/m;
            X = reshape(X, [m bR m bC]);
            X = permute(X, [1 3 2 4]);
            X = reshape(X, [m*m bR*bC]);
            db = reshape(permute(self.dualBasis, [3 2 1]), [d m*m]);
            % ^ transpose as we did not write [3 1 2]
            b = reshape(self.basis, [m*m d]);
            X = b*(db*X);
            X = reshape(X, [m m bR bC]);
            X = permute(X, [1 3 2 4]);
            X = reshape(X, [m*bR m*bC]);
        end
        
        function X = toMatrix(self, x)
            d = self.d;
            m = self.m;
            X = reshape(self.basis, [m*m d]) * x;
            X = reshape(X, [m m]);
        end
        
        function x = fromMatrix(self, X)
            d = self.d;
            m = self.m;
            x = reshape(X.', [1 m*m]) * reshape(self.dualBasis, [m*m d]);
            x = x.';
        end
            
        function z = compose(self, x, y)
            d = self.d;
            z = x.'*reshape(self.gamma, [d d*d]);
            z = y.'*reshape(z, [d d]);
            z = z.';
        end
        
    end
    
    methods
        
        function law_matrix_isomorphism_compose_DD(self, x, y)
            X = self.toMatrix(x);
            Y = self.toMatrix(y);
            XY = X*Y;
            XY1 = self.toMatrix(self.compose(x, y));
            self.assertTrue(isequal(XY, XY1));
        end
        
        function law_matrix_isomorphism_D(self, x)
            X = self.toMatrix(x);
            x1 = self.fromMatrix(X);
            self.assertEqv(x, x1);
        end
        
    end
    
    methods (Static)
        
        function D = real
            D = replab.DivisionAlgebra('Real division algebra', 'R', 1, 1, 1, 1);
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
            dualBasis(:,:,2) = [0  1
                                -1 0]/2;
            identity = [1; 0];
            D = replab.DivisionAlgebra('Complex division algebra', 'C', gamma, identity, basis, dualBasis);
        end
        
        function D = quaternion
        % Quaternion division algebra
        %
        % Multiplication table
        %     1   i   j   k
        % 1   1   i   j   k
        % i   i  -1   k  -j
        % j   j  -k  -1   i
        % k   k   j  -i  -1
            gamma = zeros(4,4,4);
            gamma(:,:,1) = [ 1  0  0  0
                             0 -1  0  0
                             0  0 -1  0
                             0  0  0 -1];
            gamma(:,:,2) = [ 0  1  0  0
                             1  0  0  0
                             0  0  0  1
                             0  0 -1  0];
            gamma(:,:,3) = [ 0  0  1  0
                             0  0  0 -1
                             1  0  0  0
                             0  1  0  0];
            gamma(:,:,4) = [ 0  0  0  1
                             0  0  1  0
                             0 -1  0  0
                             1  0  0  0];
            % The matrix representation looks like
            % [a -b -c -d
            %  b  a -d  c
            %  c  d  a -b
            %  d -c  b  a]

            basis(:,:,1) = [ 1  0  0  0
                             0  1  0  0
                             0  0  1  0
                             0  0  0  1];
            basis(:,:,2) = [ 0 -1  0  0
                             1  0  0  0
                             0  0  0 -1
                             0  0  1  0];
            basis(:,:,3) = [ 0  0 -1  0
                             0  0  0  1
                             1  0  0  0
                             0 -1  0  0];
            basis(:,:,4) = [ 0  0  0 -1
                             0  0 -1  0
                             0  1  0  0
                             1  0  0  0];
            dualBasis = permute(basis, [2 1 3])/4;
            identity = [1;0;0;0];
            D = replab.DivisionAlgebra('Quaternion division algebra', 'H', gamma, identity, basis, dualBasis);
        end
        
    end

end
