classdef Matrices < replab.domain.VectorSpace
% Describes the vector space of nR x nC real/complex matrices

    properties
        nR % (integer): Row size
        nC % (integer): Column size
    end

    methods

        %% Own methods

        function self = Matrices(field, nR, nC)
        % Constructs the space of nR x nC real/complex matrices
        %
        % Args:
        %   field ({'R', 'C'}): Real (R) or complex (C) coefficients
        %   nR (integer): Number of rows
        %   nC (integer): Number of columns
            self.field = field;
            self.nR = nR;
            self.nC = nC;
        end

        %% Str methods

        function s = headerStr(self)
            if self.overR
                s = sprintf('%d x %d real matrices', self.nR, self.nC);
            else
                s = sprintf('%d x %d complex matrices', self.nR, self.nC);
            end
        end

        %% Domain methods

        function b = eqv(self, X, Y)
            b = ~replab.isNonZeroMatrix(X - Y, replab.Parameters.doubleEigTol);
        end

        function X = sample(self)
            if self.overR
                X = randn(self.nR, self.nC);
            else
                realPart = randn(self.nR, self.nC);
                imagPart = randn(self.nR, self.nC);
                X = (realPart + 1i * imagPart)/sqrt(2);
            end
        end

    end

end
