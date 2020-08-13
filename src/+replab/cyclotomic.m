classdef cyclotomic
% Matrix with elements in the cyclotomic field
%
% Requires Java, and the cyclolab support library.
%
% The implementation below is incomplete:
%
% * Only some binary operatoins convert one of the arguments from double to cyclotomic if necessary.
%
% * Broadcasting when using binary operations is incomplete (i.e. ``matrix op scalar`` or ``scalar op matrix``)
%
% * Matrix division (i.e. solving linear systems) is not implemented
%
% * The underlying library is pretty slow in general, the speed is mostly comparable with the symbolic toolbox
%
% * Only basic operations are implemented; additional features will be added as required
%
% Example:
%   >>> I = replab.cyclotomic.eye(3)
%       I =
%       1  0  0
%       0  1  0
%       0  0  1

    properties (SetAccess = protected)
        mat % (cell(\*,\*) of cyclo.Cyclo): Matrix of cyclotomic numbers
    end

    methods (Static)

        function c = eye(n)
        % Constructs the ``n x n`` identity matrix
        %
        % Args:
        %    n (integer): Matrix size
        %
        % Returns:
        %    `.cyclotomic`: The constructed matrix
            c = replab.cyclotomic.zeros(n, n);
            mat = c.mat;
            one = javaMethod('one', 'cyclo.Cyclo');
            for i = 1:n
                mat{i,i} = one;
            end
            c = replab.cyclotomic(mat);
        end

        function c = zeros(m, n)
        % Constructs a cyclotomic matrix filled with zeros
        %
        % Args:
        %    m (integer): Number of rows
        %    n (integer): Number of columns
        %
        % Returns:
        %    `.cyclotomic`: The constructed matrix
            mat = repmat({javaMethod('zero', 'cyclo.Cyclo')}, m, n);
            c = replab.cyclotomic(mat);
        end

        function c = ones(m, n)
        % Constructs a cyclotomic matrix filled with ones
        %
        % Args:
        %    m (integer): Number of rows
        %    n (integer): Number of columns
        %
        % Returns:
        %    `.cyclotomic`: The constructed matrix
            mat = repmat({javaMethod('one', 'cyclo.Cyclo')}, m, n);
            c = replab.cyclotomic(mat);
        end

        function c = E(n)
        % Returns a root of unity
        %
        % Args:
        %    n (integer): Root order
        %
        % Returns:
        %    `.cyclotomic`: The value ``exp(2i*pi/n)``
            c = replab.cyclotomic({javaMethod('e', 'cyclo.Cyclo')});
        end

        function c = fromStrings(strings)
        % Constructs a cyclotomic matrix from its string representation
        %
        % Example:
        %   >>> replab.cyclotomic.fromStrings({'1' '1/2'; '1/2' '1'})
        %         1   1/2
        %        1/2   1
        %
        % Args:
        %    strings (cell(\*,\*) of charstring): Expressions
        %
        % Returns:
        %    `.cyclotomic`: The constructed matrix
            ja = javaMethod('parse', 'cyclo.Lab', strings(:));
            c = replab.cyclotomic.fromJavaArray(ja, size(strings));
        end

        function c = fromDoubles(doubles)
        % Constructs a cyclotomic matrix from doubles
        %
        % Note that only fractions with a power-of-two denominator can be represented exactly.
        %
        % Example:
        %   >>> replab.cyclotomic.fromDoubles([1/2])
        %       1/2
        %   >>> replab.cyclotomic.fromDoubles([1/3])
        %       6004799503160661/18014398509481984
        %
        % Args:
        %    doubles (double(\*,\*)): Matrix to convert
        %
        % Returns:
        %    `.cyclotomic`: The constructed matrix
            c = replab.cyclotomic.fromJavaArray(javaMethod('fromDouble', 'cyclo.Lab', doubles(:)), size(doubles));
        end

        function c = fromRationals(numerators, denominators)
        % Constructs a cyclotomic matrix of rational numbers
        %
        % Example:
        %   >>> replab.cyclotomic.fromRationals([1 1; 1 1], [2 3; 3 2])
        %       1/2  1/3
        %       1/3  1/2
        %
        % Args:
        %    numerators (double(\*,\*)): Coefficient numerators
        %    denominators (double(\*,\*)): Coefficient denominators
        %
        % Returns:
        %    `.cyclotomic`: The constructed matrix
            assert(isequal(size(numerators), size(denominators)));
            c = replab.cyclotomic.fromJavaArray(javaMethod('fromRational', 'cyclo.Lab', numerators(:), denominators(:)), size(numerators));
        end

        function c = sqrtRational(num, den)
        % Returns the square root of a rational number
        %
        % Example:
        %   >>> replab.cyclotomic.sqrtRational(2)
        %       E(8) - E(8)^3
        %   >>> replab.cyclotomic.sqrtRational(1,2)
        %       E(8)/2 - E(8)^3/2
        %
        % Args:
        %    num (integer): Numerator
        %    den (integer, optional): Denominator, default value ``1``
        %
        % Returns:
        %    `.cyclotomic`: The value ``sqrt(num/den)``
            if nargin == 1
                c = replab.cyclotomic({javaMethod('sqrt', 'cyclo.Lab', num)});
            else
                c = replab.cyclotomic({javaMethod('sqrt', 'cyclo.Lab', num, den)});
            end
        end

        function c = fromJavaArray(ja, sz)
        % Creates a cyclotomic matrix from the data in a Java array
        %
        % Args:
        %   ja (1D Java array of ``cyclo.Cyclo``): Matrix elements
        %   sz (integer(1, 2)): Size of the cyclotomic matrix
        %
        % Returns:
        %   `.cyclotomic`: Cyclotomic matrix
            c = replab.cyclotomic(reshape(replab.compat.javaArrayToCell(ja), sz));
        end

    end

    methods (Access = protected)

        function res = matArray(self)
        % Returns the data of this cyclotomic matrix as a 1D Java array
            c = self.mat;
            c = c(:);
            if prod(size(self.mat)) == 1 || replab.compat.isOctave
                res = javaArray('cyclo.Cyclo', length(c));
                for i = 1:length(c)
                    res(i) = c{i};
                end
            else
                res = [c{:}];
            end
        end

        function self = cyclotomic(mat)
        % Constructs a cyclotomic matrix from a cell array of coefficients
        %
        % Args:
        %   mat (cell(\*,\*) of Java ``cyclo.Cyclo``): Coefficients
            assert(iscell(mat));
            self.mat = mat;
        end

    end

    methods

        function disp(self)
        % Standard display method
            t = replab.compat.javaArrayToCell(javaMethod('print', 'cyclo.Lab', self.matArray));
            t = replab.str.Table(reshape(t, self.size), 'uniform', 0);
            disp(t);
        end

        function s = num2str(self)
        % Conversion to string
            t = replab.compat.javaArrayToCell(javaMethod('print', 'cyclo.Lab', self.matArray));
            t = replab.str.Table(reshape(t, self.size), 'uniform', 0);
            s = t.format(1000, 1000);
        end

        function res = ne(self, rhs)
        % (Non-)equality test
            res = ~(self == rhs);
        end

        function res = eq(lhs, rhs)
        % Equality test
            if isa(lhs, 'double')
                lhs = replab.cyclotomic.fromDoubles(lhs);
            end
            if isa(rhs, 'double')
                rhs = replab.cyclotomic.fromDoubles(rhs);
            end
            res = reshape(javaMethod('eqv', 'cyclo.Lab', lhs.matArray, rhs.matArray), size(lhs.mat));
        end

        function res = conj(self)
        % Complex conjugation
            res = replab.cyclotomic.fromJavaArray(javaMethod('conjugate', 'cyclo.Lab', self.matArray), self.size);
        end

        function res = sqrt(self)
        % Square root
        %
        % Requires that the argument contains rational coefficients.
            res = replab.cyclotomic.fromJavaArray(javaMethod('sqrt', 'cyclo.Lab', self.matArray), self.size);
        end

        function res = plus(lhs, rhs)
        % Standard ``+`` operator
        %
        % Does not support broadcasting (i.e. ``M + 1`` when ``M`` is not a scalar)
            if isa(lhs, 'double')
                lhs = replab.cyclotomic.fromDoubles(lhs);
            end
            if isa(rhs, 'double')
                rhs = replab.cyclotomic.fromDoubles(rhs);
            end
            res = replab.cyclotomic.fromJavaArray(javaMethod('plus', 'cyclo.Lab', lhs.matArray, rhs.matArray), size(lhs));
        end

        function res = uminus(self)
            res = replab.cyclotomic.fromJavaArray(javaMethod('negate', 'cyclo.Lab', self.matArray), size(self.matArray));
        end

        function res = minus(lhs, rhs)
        % Standard ``-`` operator
        %
        % Does not support broadcasting (i.e. ``M + 1`` when ``M`` is not a scalar)
            if isa(lhs, 'double')
                lhs = replab.cyclotomic.fromDoubles(lhs);
            end
            if isa(rhs, 'double')
                rhs = replab.cyclotomic.fromDoubles(rhs);
            end
            res = replab.cyclotomic.fromJavaArray(javaMethod('minus', 'cyclo.Lab', lhs.matArray, rhs.matArray), size(lhs));
        end

        function res = times(lhs, rhs)
        % Pointwise ``*`` operator
        %
        % Does not support broadcasting (i.e. ``M * 2`` when ``M`` is not a scalar)
            if isa(lhs, 'double')
                lhs = replab.cyclotomic.fromDoubles(lhs);
            end
            if isa(rhs, 'double')
                rhs = replab.cyclotomic.fromDoubles(rhs);
            end
            mat = reshape(cell(com.faacets.gluon.Cyclotomic.pw_times(lhs.matArray, rhs.matArray)), size(lhs.mat));
            res = replab.cyclotomic(mat);
        end

        function res = mrdivide(lhs, rhs)
        % Division
        %
        % Only supports the case where the right hand side is scalar
            if isa(lhs, 'double')
                lhs = replab.cyclotomic.fromDoubles(lhs);
            end
            if isa(rhs, 'double')
                rhs = replab.cyclotomic.fromDoubles(rhs);
            end
            assert(isscalar(rhs.mat), 'The / operator is only implemented for scalar rhs');
            rm = rhs.mat;
            res = replab.cyclotomic.fromJavaArray(javaMethod('divideScalar', 'cyclo.Lab', lhs.matArray, rm{1}), lhs.size);
        end

        function res = transpose(self)
            res = replab.cyclotomic(self.mat.');
        end

        function res = ctranspose(self)
            res = conj(transpose(self));
        end

        function res = mtimes(lhs, rhs)
        % Matrix multiplication
        %
        % Support both ``m x n`` by ``n x p`` matrix multiplication, and the case where one of the arguments is a scalar
            if isa(lhs, 'double')
                lhs = replab.cyclotomic.fromDoubles(lhs);
            end
            if isa(rhs, 'double')
                rhs = replab.cyclotomic.fromDoubles(rhs);
            end
            if ~isscalar(rhs.mat) && isscalar(lhs.mat)
                res = rhs * lhs;
                return
            end
            if isscalar(rhs.mat)
                l = size(lhs.mat, 1);
                m = size(lhs.mat, 2);
                rm = rhs.mat;
                res = replab.cyclotomic.fromJavaArray(javaMethod('timesScalar', 'cyclo.Lab', lhs.matArray, rm{1,1}), size(lhs));
                return
            end
            l = size(lhs.mat, 1);
            m = size(lhs.mat, 2);
            assert(m == size(lhs.mat, 1));
            n = size(rhs.mat, 2);
            res = replab.cyclotomic.fromJavaArray(javaMethod('times', 'cyclo.Lab', l, m, n, lhs.matArray, rhs.matArray), [l n]);
        end

        function res = mpower(self, e)
        % Matrix power
            n = size(self.mat, 1);
            assert(size(self.mat, 2) == n);
            res = replab.cyclotomic.fromJavaArray(javaMethod('power', 'cyclo.Lab', n, self.matArray, e), [n n]);
        end

        function res = inv(self)
        % Matrix inverse
        %
        % Requires that the matrix is full rank
            n = size(self.mat, 1);
            assert(size(self.mat, 2) == n);
            res = replab.cyclotomic.fromJavaArray(javaMethod('inverse', 'cyclo.Lab', n, self.matArray), [n n]);
        end

        function s = size(self, varargin)
        % Matrix size
            s = size(self.mat, varargin{:});
        end

        function l = length(self)
        % Matrix length
            l = length(self.mat);
        end

        function varargout = subsref(self, s)
        % Indexing
            switch s(1).type
              case '.'
                [varargout{1:nargout}] = builtin('subsref', self, s);
              case '()'
                res = replab.cyclotomic(subsref(self.mat, s(1)));
                if length(s) > 1
                    res = subsref(res, s(2:end));
                end
                varargout{1} = res;
              otherwise
                error('Not a valid indexing expression')
            end
        end

        function self = subsasgn(self, s, val)
        % Assignment
            switch s(1).type
              case '.'
                self = builtin('subsasgn', self, s, val);
              case '()'
                if length(s) == 1
                    if isa(val, 'replab.cyclotomic')
                        self.mat = builtin('subsasgn', self.mat, s(1), val.mat);
                    elseif isa(val, 'double')
                        arg = replab.cyclotomic.fromDoubles(val);
                        self = subsasgn(self, s(1), arg);
                    end
                end
              otherwise
                error('Not a valid indexing expression')
            end
        end

        function res = double(self)
        % Conversion to floating-point double
            pairs = javaMethod('toDouble', 'cyclo.Lab', self.matArray);
            if replab.compat.isOctave
                R = pairs(1);
                I = pairs(2);
            else
                R = pairs(1,:);
                I = pairs(2,:);
            end
            res = reshape(R, size(self.mat)) + 1i * reshape(I, size(self.mat));
        end

        function res = horzcat(self, varargin)
        % Horizontal concatenation
            rhs = cellfun(@(a) a.mat, varargin, 'uniform', 0);
            res = horzcat(self.mat, rhs{:});
            res = replab.cyclotomic(res);
        end

        function res = vertcat(self, varargin)
        % Vertical concatenation
            rhs = cellfun(@(a) a.mat, varargin, 'uniform', 0);
            res = vertcat(self.mat, rhs{:});
            res = replab.cyclotomic(res);
        end

    end

end
