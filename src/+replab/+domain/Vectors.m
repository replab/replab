classdef Vectors < replab.domain.VectorSpace
% Describes the vector space of d-dimensional real/complex column vectors

    properties
        d % integer: Dimension
    end

    methods

        %% Own methods

        function self = Vectors(field, d)
        % Constructs the space of d-dimensional real/complex column vectors
        %
        % Args:
        %   field ({'R', 'C'}): Real (R) or complex (C) coefficients
        %   d (integer): Dimension
            self.field = field;
            self.d = d;
        end

        %% Str methods

        function s = headerStr(self)
            if self.overR
                s = sprintf('%d x 1 real vectors', self.d);
            else
                s = sprintf('%d x 1 complex vectors', self.d);
            end
        end

        %% Domain methods

        function b = eqv(self, X, Y)
            b = ~replab.numerical.isNonZeroMatrix(X - Y, replab.globals.doubleEigTol);
        end

        function X = sample(self)
            if self.overR
                X = randn(self.d, 1);
            else
                realPart = randn(self.d, 1);
                imagPart = randn(self.d, 1);
                X = (realPart + 1i * imagPart)/sqrt(2);
            end
        end

    end

end
