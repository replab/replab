classdef cyclotomic
% Multidimensional array with elements in the cyclotomic field
%
% Requires Java and the cyclolab support library.
%
% The implementation is pretty slow, but values are exact.
%
% Example:
%   >>> I = replab.cyclotomic.eye(3)
%       I =
%       1  0  0
%       0  1  0
%       0  0  1

    properties (Access = protected)
        size_ % (integer(1,\*)): Shape
        data_ % (``cyclo.Cyclo[]``): 1D array of cyclotomic numbers, elements are stored as in ``array(:)``
    end

    methods (Static)

        function c = sparse(I, J, K, nR, nC)
        % Constructs a 2D cyclotomic array using sparse matrix syntax (incomplete)
        %
        % Only supports the five argument syntax
        %
        % Returns:
        %    `.cyclotomic`: The constructed matrix
            if isa(K, 'double')
                K = replab.cyclotomic(K);
            end
            c = replab.cyclotomic.zeros(nR, nC);
            for i = 1:length(I)
                c(I(i), J(i)) = K(i);
            end
        end

        function c = eye(n)
        % Constructs the ``n x n`` identity matrix (incomplete)
        %
        % Does not support the syntax with two arguments.
        %
        % Args:
        %    n (integer): Matrix size
        %
        % Returns:
        %    `.cyclotomic`: The constructed matrix
            array = repmat({javaMethod('zero', 'cyclo.Cyclo')}, n, n);
            one = javaMethod('one', 'cyclo.Cyclo');
            for i = 1:n
                array{i,i} = one;
            end
            c = replab.cyclotomic(array);
        end

        function c = zeros(varargin)
        % Constructs a cyclotomic array filled with zeros
        %
        % Returns:
        %    `.cyclotomic`: The constructed matrix
            array = repmat({javaMethod('zero', 'cyclo.Cyclo')}, varargin{:});
            c = replab.cyclotomic(array);
        end

        function c = ones(varargin)
        % Constructs a cyclotomic array filled with ones
        %
        % Returns:
        %    `.cyclotomic`: The constructed matrix
            array = repmat({javaMethod('one', 'cyclo.Cyclo')}, varargin{:});
            c = replab.cyclotomic(array);
        end

        function c = E(orders)
        % Returns roots of unity
        %
        % Args:
        %    n (integer(\*,\*)): Root orders
        %
        % Returns:
        %    `.cyclotomic`: The value ``exp(2i*pi/n)``
            ja = javaMethod('E', 'cyclo.Lab', orders);
            c = replab.cyclotomic.fromJavaArray(ja, size(orders));
        end

% $$$         function c = approximate(lowerBounds, upperBounds)
% $$$         % Constructs the simplest rational approximation within an interval, for coefficients of a matrix
% $$$         %
% $$$         % Example:
% $$$         %   >>> replab.cyclotomic.approximate(3.141592, 3.141593)
% $$$         %     355/113
% $$$         %
% $$$         % Args:
% $$$         %   lowerBounds (double(\*,\*)): Interval lower bounds
% $$$         %   upperBounds (double(\*,\*)): Interval upper bounds
% $$$         %
% $$$         % Returns:
% $$$         %   `.cyclotomic`: The simplest rational approximation, coefficient by coefficient
% $$$             assert(isequal(size(lowerBounds), size(upperBounds)));
% $$$             c = replab.cyclotomic.fromJavaArray(javaMethod('approximate', 'cyclo.Lab', lowerBounds(:), upperBounds(:)), size(lowerBounds));
% $$$         end
% $$$
% $$$
% $$$         function c = fromRationals(numerators, denominators)
% $$$         % Constructs a cyclotomic matrix of rational numbers
% $$$         %
% $$$         % Example:
% $$$         %   >>> replab.cyclotomic.fromRationals([1 1; 1 1], [2 3; 3 2])
% $$$         %       1/2  1/3
% $$$         %       1/3  1/2
% $$$         %
% $$$         % Args:
% $$$         %    numerators (double(\*,\*)): Coefficient numerators
% $$$         %    denominators (double(\*,\*)): Coefficient denominators
% $$$         %
% $$$         % Returns:
% $$$         %    `.cyclotomic`: The constructed matrix
% $$$             assert(isequal(size(numerators), size(denominators)));
% $$$             c = replab.cyclotomic.fromJavaArray(javaMethod('fromRational', 'cyclo.Lab', numerators(:), denominators(:)), size(numerators));
% $$$         end
% $$$
% $$$         function c = sqrtRational(num, den)
% $$$         % Returns the square root of a rational number
% $$$         %
% $$$         % Example:
% $$$         %   >>> replab.cyclotomic.sqrtRational(2)
% $$$         %       E(8) - E(8)^3
% $$$         %   >>> replab.cyclotomic.sqrtRational(1,2)
% $$$         %       E(8)/2 - E(8)^3/2
% $$$         %
% $$$         % Args:
% $$$         %    num (integer): Numerator
% $$$         %    den (integer, optional): Denominator, default value ``1``
% $$$         %
% $$$         % Returns:
% $$$         %    `.cyclotomic`: The value ``sqrt(num/den)``
% $$$             if nargin == 1
% $$$                 c = replab.cyclotomic({javaMethod('sqrt', 'cyclo.Lab', num)});
% $$$             else
% $$$                 c = replab.cyclotomic({javaMethod('sqrt', 'cyclo.Lab', num, den)});
% $$$             end
% $$$         end
% $$$
% $$$         function c = fromJavaArray(ja, sz)
% $$$         % Creates a cyclotomic matrix from the data in a Java array
% $$$         %
% $$$         % Args:
% $$$         %   ja (1D Java array of ``cyclo.Cyclo``): Matrix elements
% $$$         %   sz (integer(1, 2)): Size of the cyclotomic matrix
% $$$         %
% $$$         % Returns:
% $$$         %   `.cyclotomic`: Cyclotomic matrix
% $$$             c = replab.cyclotomic(reshape(replab.compat.javaArrayToCell(ja), sz));
% $$$         end

    end

    methods (Static, Access = protected)

        function c = convertVpi(vec)
            vec = vec(:);
            s = arrayfun(@(i) strtrim(num2str(v(i))), 1:length(v), 'uniform', 0);
            c = javaMethod('parse', 'cyclo.Lab', s);
        end

        function c = convertCell(vec)
        % Constructs a cell vector to a Java array
        %
        % Args:
        %    vec (cell(\*,1) of double, charstring, vpi, cyclotomic): Vector to convert
        %
        % Returns:
        %    ``cyclo.Cyclo[]``: A Java array of cyclo.Cyclo instances
            vec = vec(:);
            realMask = cellfun(@(v) isa(v, 'double') && isreal(v), vec);
            complexMask = cellfun(@(v) isa(v, 'double') && ~isreal(v), vec);
            charMask = cellfun(@iscell, vec);
            vpiMask = cellfun(@(v) isa(v, 'vpi'), vec);
            cyclotomicMask = cellfun(@(v) isa(v, 'replab.cyclotomic'), vec);
            cycloMask = cellfun(@(v) isa(v, 'cyclo.Cyclo'), vec);
            restMask = ~(realMask | complexMask | charMask | vpiMask | cyclotomicMask | cycloMask);
            if any(restMask)
                ind = find(restMask);
                types = strjoin(unique(cellfun(@class, vec(ind), 'uniform', 0)), ', ');
                error('Array to convert contains unsupported types: %s', types);
            end
            c = javaArray('cyclo.Cyclo', length(vec));
            ind = find(realMask);
            if ~isempty(ind)
                c(ind) = javaMethod('fromDouble', 'cyclo.Lab', [vec{ind}]);
            end
            ind = find(complexMask);
            if ~isempty(ind)
                c(ind) = javaMethod('fromDouble', 'cyclo.Lab', real([vec{ind}]), imag([vec{ind}]));
            end
            ind = find(charMask);
            if ~isempty(ind)
                c(ind) = javaMethod('parse', 'cyclo.Lab', vec(ind));
            end
            ind = find(vpiMask);
            if ~isempty(ind)
                c(ind) = javaMethod('parse', 'cyclo.Lab', cellfun(@(v) strtrim(num2str(v)), vec(ind), 'uniform', 0));
            end
            ind = find(cyclotomicMask);
            if ~isempty(ind)
                for i = ind(:)'
                    c(i) = vec{i}.cyclo(1);
                end
            end
            ind = find(cycloMask);
            if ~isempty(ind)
                for i = ind(:)'
                    c(i) = vec{i};
                end
            end
        end

        function c = convertDouble(vec)
        % Converts a vector of floating-point doubles to a Java array
        %
        % Args:
        %    vec (double(\*,1)): Vector to convert
        %
        % Returns:
        %    ``cyclo.Cyclo[]``: A Java array of cyclo.Cyclo instances
            if isreal(vec)
                c = javaMethod('fromDouble', 'cyclo.Lab', vec(:));
            else
                c = javaMethod('fromDouble', 'cyclo.Lab', real(vec(:)), imag(vec(:)));
            end
        end

        function c = convert(vec)
        % Converts a vector to a Java array
        % Args:
        %    vec (double(\*,1) or vpi(\*,1) or cell(\*,1)): Vector to convert
        %
        % Returns:
        %    ``cyclo.Cyclo[]``: A Java array of cyclo.Cyclo instances
            if isa(vec, 'double')
                c = replab.cyclotomic.convertDouble(vec);
            elseif isa(vec, 'vpi')
                c = replab.cyclotomic.convertVpi(vec);
            elseif iscell(vec)
                c = replab.cyclotomic.convertCell(vec);
            end
        end

% $$$         function arg = shapeArg(arg)
% $$$             if isa(arg, 'double')
% $$$                 arg = replab.cyclotomic.fromDoubles(arg);
% $$$             end
% $$$         end
% $$$
% $$$         function [lhs rhs] = shapeArgs(lhs, rhs)
% $$$             if isa(lhs, 'double')
% $$$                 lhs = replab.cyclotomic.fromDoubles(lhs);
% $$$             end
% $$$             if isa(rhs, 'double')
% $$$                 rhs = replab.cyclotomic.fromDoubles(rhs);
% $$$             end
% $$$             if isscalar(lhs) && ~isscalar(rhs)
% $$$                 lhs = replab.cyclotomic(repmat(lhs.mat, size(rhs)));
% $$$             end
% $$$             if ~isscalar(lhs) && isscalar(rhs)
% $$$                 rhs = replab.cyclotomic(repmat(rhs.mat, size(lhs)));
% $$$             end
% $$$         end

    end

    methods

        function self = cyclotomic(array, size_)
        % Constructs a cyclotomic array from an array of coefficients
        %
        % We can construct cyclotomic from floating-point numbers. Note that only fractions with a power-of-two
        % denominator can be represented exactly.
        %
        % Example:
        %   >>> replab.cyclotomic(1/2)
        %       1/2
        %   >>> replab.cyclotomic(1/3)
        %       6004799503160661/18014398509481984
        %
        % The constructor also accepts heterogenous cell arrays:
        %
        % Example:
        %   >>> replab.cyclotomic({'1' '1/2'; '1/2' '1'})
        %         1   1/2
        %        1/2   1
        %
        % Args:
        %   array: Coefficients
            if isa(array, 'replab.cyclotomic')
                self.data_ = array.data;
                self.size_ = size(array);
            elseif isa(array, 'double') || isa(array, 'vpi') || iscell(array)
                assert(nargin == 1, 'The two argument call is reserved when the array is of type cyclo.Cyclo[]');
                self.data_ = replab.cyclotomic.convert(array(:));
                self.size_ = size(array);
            elseif isa(array, 'cyclo.Cyclo[]')
                self.data_ = array;
                self.size_ = size_;
            elseif isa(array, 'cyclo.Cyclo')
                data = javaArray('cyclo.Cyclo', 1);
                data(1) = array;
                self.data_ = data;
                self.size_ = size_;
            else
                error('Incorrect constructor call');
            end
        end

        function c = cyclo(self, varargin)
        % Returns the cyclo.Cyclo object at the given index in this multidimensional array
        %
        % Arguments are integer indices, as passed to subsref
            l = length(varargin)
            if l == 1
                c = self.data_(varargin{1});
            elseif l == length(self.size_)
                c = self.data_(sub2ind(self.size_, varargin{:}));
            else
                error('Wrong number of indices');
            end
        end

        function c = data(self)
        % Returns the java array containing the flat data of this array
            c = javaMethod('cloneArray', 'cyclo.Lab', self.data_);
        end

    end

    methods % Display

        function disp(self)
        % Standard display method
            t = replab.compat.javaArrayToCell(javaMethod('print', 'cyclo.Lab', self.data_));
            t = replab.str.Table(reshape(t, size(self)), 'uniform', 0);
            disp(t);
        end

        function s = num2str(self)
        % Conversion to string (incomplete)
        %
        % Does not support format specifiers or a specified precision.
        %
        % Returns:
        %   charstring: Possibly multiline string representation of the matrix
            t = replab.compat.javaArrayToCell(javaMethod('print', 'cyclo.Lab', self.data_));
            t = replab.str.Table(reshape(t, size(self)), 'uniform', 0);
            s = t.format(1000, 1000);
        end

    end

    methods % Standard MATLAB properties

% $$$         function res = isrational(self)
% $$$         % Returns which coefficients are rational
% $$$             res = reshape(javaMethod('isRational', 'cyclo.Lab', self.matArray), self.size);
% $$$         end
% $$$
% $$$         function res = isreal(self)
% $$$             res = all(all(self == conj(self)));
% $$$         end
% $$$
% $$$         function res = isscalar(self)
% $$$             res = isscalar(self.mat);
% $$$         end
% $$$
% $$$         function res = isvector(self)
% $$$             res = isvector(self.mat);
% $$$         end
% $$$
% $$$         function res = iswhole(self)
% $$$         % Returns which coefficients are integers
% $$$             res = reshape(javaMethod('isWhole', 'cyclo.Lab', self.matArray), self.size);
% $$$         end
% $$$
        function l = length(self)
        % Length of vector
        %
        % Equivalent to ``max(size(X))`` if ``isempty(X)`` is false, otherwise ``0``.
        %
        % Example:
        %   >>> length(replab.cyclotomic.zeros(10))
        %       10
        %
        % Returns:
        %   integer: Vector length
            if isempty(self)
                l = 0;
            else
                l = max(size(self));
            end
        end

        function s = size(self, dim)
        % Size of array (incomplete)
        %
        % Incomplete implementation: only supports the forms ``s = size(array)`` and ``s = size(array, i)`` where
        % ``i`` is an integer scalar.
        %
        % Example:
        %   >>> size(replab.cyclotomic.eye(2))
        %       2     2
        %
        % Args:
        %   dim (integer, optional): Dimension(s) to get the size of
        %
        % Returns:
        %   integer or integer(1,\*): Size
            if nargin < 2
                s = self.size_;
            else
                s = self.size_(dim);
            end
        end

    end

    methods % Binary operations

% $$$         function res = dot(lhs, rhs)
% $$$             res = sum(conj(lhs(:)).*rhs(:));
% $$$         end
% $$$
% $$$         function res = eq(lhs, rhs)
% $$$         % Equality test
% $$$             [lhs rhs] = replab.cyclotomic.shapeArgs(lhs, rhs);
% $$$             res = reshape(javaMethod('eqv', 'cyclo.Lab', lhs.matArray, rhs.matArray), size(lhs.mat));
% $$$         end
% $$$
% $$$         function res = minus(lhs, rhs)
% $$$         % Standard ``-`` operator
% $$$             [lhs rhs] = replab.cyclotomic.shapeArgs(lhs, rhs);
% $$$             res = replab.cyclotomic.fromJavaArray(javaMethod('minus', 'cyclo.Lab', lhs.matArray, rhs.matArray), size(lhs));
% $$$         end
% $$$
% $$$         function res = mpower(self, e)
% $$$         % Matrix power
% $$$             n = size(self.mat, 1);
% $$$             assert(size(self.mat, 2) == n);
% $$$             res = replab.cyclotomic.fromJavaArray(javaMethod('power', 'cyclo.Lab', n, self.matArray, e), [n n]);
% $$$         end
% $$$
% $$$         function res = mrdivide(lhs, rhs)
% $$$         % Division
% $$$         %
% $$$         % Only supports the case where the right hand side is scalar
% $$$             if isa(lhs, 'double')
% $$$                 lhs = replab.cyclotomic.fromDoubles(lhs);
% $$$             end
% $$$             if isa(rhs, 'double')
% $$$                 rhs = replab.cyclotomic.fromDoubles(rhs);
% $$$             end
% $$$             assert(isscalar(rhs.mat), 'The / operator is only implemented for scalar rhs');
% $$$             rm = rhs.mat;
% $$$             res = replab.cyclotomic.fromJavaArray(javaMethod('divideScalar', 'cyclo.Lab', lhs.matArray, rm{1}), lhs.size);
% $$$         end
% $$$
% $$$         function res = mtimes(lhs, rhs)
% $$$         % Matrix multiplication
% $$$         %
% $$$         % Support both ``m x n`` by ``n x p`` matrix multiplication, and the case where one of the arguments is a scalar
% $$$             if isa(lhs, 'double')
% $$$                 lhs = replab.cyclotomic.fromDoubles(lhs);
% $$$             end
% $$$             if isa(rhs, 'double')
% $$$                 rhs = replab.cyclotomic.fromDoubles(rhs);
% $$$             end
% $$$             if ~isscalar(rhs.mat) && isscalar(lhs.mat)
% $$$                 res = rhs * lhs;
% $$$                 return
% $$$             end
% $$$             if isscalar(rhs.mat)
% $$$                 l = size(lhs.mat, 1);
% $$$                 m = size(lhs.mat, 2);
% $$$                 rm = rhs.mat;
% $$$                 res = replab.cyclotomic.fromJavaArray(javaMethod('timesScalar', 'cyclo.Lab', lhs.matArray, rm{1,1}), size(lhs));
% $$$                 return
% $$$             end
% $$$             l = size(lhs.mat, 1);
% $$$             m = size(lhs.mat, 2);
% $$$             assert(m == size(rhs.mat, 1));
% $$$             n = size(rhs.mat, 2);
% $$$             res = replab.cyclotomic.fromJavaArray(javaMethod('times', 'cyclo.Lab', l, m, n, lhs.matArray, rhs.matArray), [l n]);
% $$$         end
% $$$
% $$$         function res = ne(self, rhs)
% $$$         % (Non-)equality test
% $$$             res = ~(self == rhs);
% $$$         end
% $$$
% $$$         function res = plus(lhs, rhs)
% $$$         % Standard ``+`` operator
% $$$         %
% $$$         % Does not support broadcasting (i.e. ``M + 1`` when ``M`` is not a scalar)
% $$$             [lhs rhs] = replab.cyclotomic.shapeArgs(lhs, rhs);
% $$$             res = replab.cyclotomic.fromJavaArray(javaMethod('plus', 'cyclo.Lab', lhs.matArray, rhs.matArray), size(lhs));
% $$$         end
% $$$
% $$$         function res = rdivide(lhs, rhs)
% $$$         % Pointwise ``/`` operator
% $$$             [lhs rhs] = replab.cyclotomic.shapeArgs(lhs, rhs);
% $$$             res = replab.cyclotomic.fromJavaArray(javaMethod('pw_divide', 'cyclo.Lab', lhs.matArray, rhs.matArray), size(lhs));
% $$$         end
% $$$
% $$$         function res = times(lhs, rhs)
% $$$         % Pointwise ``*`` operator
% $$$             [lhs rhs] = replab.cyclotomic.shapeArgs(lhs, rhs);
% $$$             res = replab.cyclotomic.fromJavaArray(javaMethod('pw_times', 'cyclo.Lab', lhs.matArray, rhs.matArray), size(lhs));
% $$$         end

    end

    methods % Shape

% $$$         function d = diag(self)
% $$$             m = min(self.size);
% $$$             matd = cell(m, 1);
% $$$             for i = 1:m
% $$$                 matd{i} = self.mat{i,i};
% $$$             end
% $$$             d = replab.cyclotomic(matd);
% $$$         end
% $$$
% $$$         function res = reshape(self, varargin)
% $$$             res = replab.cyclotomic(reshape(self.mat, varargin{:}));
% $$$         end

    end


    methods % Involutions

% $$$         function res = conj(self)
% $$$         % Complex conjugation
% $$$             res = replab.cyclotomic.fromJavaArray(javaMethod('conjugate', 'cyclo.Lab', self.matArray), self.size);
% $$$         end
% $$$
% $$$         function res = ctranspose(self)
% $$$             res = conj(transpose(self));
% $$$         end
% $$$
% $$$         function res = transpose(self)
% $$$             res = replab.cyclotomic(self.mat.');
% $$$         end
% $$$
% $$$         function res = uminus(self)
% $$$             res = replab.cyclotomic.fromJavaArray(javaMethod('negate', 'cyclo.Lab', self.matArray), size(self));
% $$$         end

    end

    methods % Matrix decompositions

% $$$         function [L, U, p] = lu(self)
% $$$             m = size(self, 1);
% $$$             n = size(self, 2);
% $$$             res = javaMethod('lu', 'cyclo.Lab', self.matArray, m, n);
% $$$             L = javaMethod('L', res);
% $$$             U = javaMethod('U', res);
% $$$             p = javaMethod('pivots', res);
% $$$             p = double(p(:)') + 1;
% $$$             if m >= n
% $$$                 L = replab.cyclotomic.fromJavaArray(L, [m n]);
% $$$                 U = replab.cyclotomic.fromJavaArray(U, [n n]);
% $$$             else
% $$$                 L = replab.cyclotomic.fromJavaArray(L, [m m]);
% $$$                 U = replab.cyclotomic.fromJavaArray(U, [m n]);
% $$$             end
% $$$         end
% $$$
% $$$         function [L, D, p] = ldl(self)
% $$$             m = size(self, 1);
% $$$             n = size(self, 2);
% $$$             assert(m == n);
% $$$             res = javaMethod('ldl', 'cyclo.Lab', self.matArray, m, n);
% $$$             L = javaMethod('L', res);
% $$$             D = javaMethod('D', res);
% $$$             p = javaMethod('pivots', res);
% $$$             p = double(p(:)') + 1;
% $$$             L = replab.cyclotomic.fromJavaArray(L, [n n]);
% $$$             D = replab.cyclotomic.fromJavaArray(D, [n n]);
% $$$         end

    end

    methods % Indexing

        function res = broadcast(self, shape)
        % Broadcast the current array to the given shape
            if length(self) == 1
                res = replab.cyclotomic(repmat({self}, shape));
            elseif isequal(size(self), shape)
                res = self;
            end
        end

        function varargout = subsref(self, s)
        % Indexing
            switch s(1).type
              case '.'
                [varargout{1:nargout}] = builtin('subsref', self, s);
              case '()'
                ind = reshape(1:prod(self.size_), self.size_);
                ind = subsref(ind, s(1));
                res = replab.cyclotomic(self.data_(ind(:)), size(ind));
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
                    ind = reshape(1:prod(self.size_), self.size_);
                    ind = subsref(ind, s(1));
                    c = self.data_;
                    val = replab.cyclotomic(val);
                    val = val.broadcast(size(ind));
                    ind = ind(:);
                    val = val.data;
                    for i = 1:length(ind)
                        c(ind(i)) = val(i);
                    end
                    self.data_ = c;
                else
                    error('Not supported');
                end
              otherwise
                error('Not a valid indexing expression')
            end
        end

    end

    methods

        function h = hash(self)
        % Returns a multidimensional array of hash codes
        %
        % Example:
        %   >>> num2str(replab.cyclotomic.zeros(1))
        %       '-1579025880'
        %
        % Returns:
        %   integer(...): Integer array of hash codes
            h = arrayfun(@(c) javaMethod('hashCode', c), self.data_);
            h = reshape(h, self.size_);
        end
% $$$
% $$$         function p = prod(self)
% $$$             assert(isscalar(self.mat) || isvector(self.mat));
% $$$             p = replab.cyclotomic({javaMethod('prod', 'cyclo.Lab', self.matArray)});
% $$$         end
% $$$
% $$$         function s = sum(self)
% $$$             assert(isscalar(self.mat) || isvector(self.mat));
% $$$             s = replab.cyclotomic({javaMethod('sum', 'cyclo.Lab', self.matArray)});
% $$$         end
% $$$
% $$$         function v = trace(self)
% $$$             v = sum(diag(self));
% $$$         end
% $$$
% $$$         function res = sqrt(self)
% $$$         % Square root
% $$$         %
% $$$         % Requires that the argument contains rational coefficients.
% $$$             res = replab.cyclotomic.fromJavaArray(javaMethod('sqrt', 'cyclo.Lab', self.matArray), self.size);
% $$$         end
% $$$
% $$$         function res = null(self)
% $$$         % Computes the null space of a cyclotomic matrix
% $$$         %
% $$$         % Example:
% $$$         %   >>> M = replab.cyclotomic.fromDoubles([1 3 0; -2 -6 0; 3 9 6]);
% $$$         %   >>> null(M)'
% $$$         %       -3  1  0
% $$$         %
% $$$         % Returns:
% $$$         %   `.cyclotomic`: The matrix null space
% $$$             rr = javaMethod('rref', 'cyclo.Lab', self.matArray, size(self, 1), size(self, 2));
% $$$             rank = double(rr.rank);
% $$$             rows = size(self, 2);
% $$$             cols = size(self, 2) - rank;
% $$$             res = replab.cyclotomic.fromJavaArray(javaMethod('nullSpace', rr), [rows cols]);
% $$$         end
% $$$
% $$$         function res = blkdiag(lhs, rhs)
% $$$             if isa(lhs, 'double')
% $$$                 lhs = replab.cyclotomic.fromDoubles(lhs);
% $$$             end
% $$$             if isa(rhs, 'double')
% $$$                 rhs = replab.cyclotomic.fromDoubles(rhs);
% $$$             end
% $$$             s1 = size(lhs);
% $$$             s2 = size(rhs);
% $$$             M = replab.cyclotomic.zeros(s1(1) + s2(1), s1(2) + s2(2));
% $$$             mat = M.mat;
% $$$             mat(1:s1(1),1:s1(2)) = lhs.mat;
% $$$             mat(s2(1)+1:end,s2(2)+1:end) = rhs.mat;
% $$$             res = replab.cyclotomic(mat);
% $$$         end
% $$$
% $$$         function res = kron(lhs, rhs)
% $$$             if isa(lhs, 'double')
% $$$                 lhs = replab.cyclotomic.fromDoubles(lhs);
% $$$             end
% $$$             if isa(rhs, 'double')
% $$$                 rhs = replab.cyclotomic.fromDoubles(rhs);
% $$$             end
% $$$             res = replab.cyclotomic.fromJavaArray(javaMethod('kron', 'cyclo.Lab', size(lhs, 1), size(lhs, 2), size(rhs, 1), size(rhs, 2), lhs.matArray, rhs.matArray), [size(lhs, 1)*size(rhs, 1) size(lhs, 2)*size(rhs, 2)]);
% $$$         end
% $$$
% $$$         function res = inv(self)
% $$$         % Matrix inverse
% $$$         %
% $$$         % Requires that the matrix is full rank
% $$$             n = size(self.mat, 1);
% $$$             assert(size(self.mat, 2) == n);
% $$$             res = replab.cyclotomic.fromJavaArray(javaMethod('inverse', 'cyclo.Lab', n, self.matArray), [n n]);
% $$$         end
% $$$
% $$$         function varargout = find(self, varargin)
% $$$             assert(isempty(varargin), 'Additional input arguments are not supported');
% $$$             mask = self ~= 0;
% $$$             switch nargout
% $$$               case 0
% $$$               case 1
% $$$                 I = find(mask);
% $$$                 varargout = {I};
% $$$               case 2
% $$$                 [I J] = find(mask);
% $$$                 varargout = {I J};
% $$$               case 3
% $$$                 [I J] = find(mask);
% $$$                 V = replab.cyclotomic(self.mat(mask));
% $$$                 varargout = {I J V};
% $$$               otherwise
% $$$                 error('Too many output arguments');
% $$$             end
% $$$         end
% $$$
% $$$
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
% $$$
% $$$         function [res error] = doubleApproximation(self)
% $$$         % Conversion to floating-point double
% $$$             data = javaMethod('toDoubleApproximation', 'cyclo.Lab', self.matArray);
% $$$             res = reshape(data.real, size(self.mat)) + 1i * reshape(data.imag, size(self.mat));
% $$$             error = reshape(data.error, size(self.mat));
% $$$         end
% $$$
% $$$         function res = horzcat(lhs, varargin)
% $$$         % Horizontal concatenation
% $$$             lhs = replab.cyclotomic.shapeArg(lhs);
% $$$             rhs = cellfun(@(a) replab.cyclotomic.shapeArg(a), varargin, 'uniform', 0);
% $$$             rhs = cellfun(@(a) a.mat, rhs, 'uniform', 0);
% $$$             res = horzcat(lhs.mat, rhs{:});
% $$$             res = replab.cyclotomic(res);
% $$$         end
% $$$
% $$$         function res = vertcat(lhs, varargin)
% $$$         % Vertical concatenation
% $$$             lhs = replab.cyclotomic.shapeArg(lhs);
% $$$             rhs = cellfun(@(a) replab.cyclotomic.shapeArg(a), varargin, 'uniform', 0);
% $$$             rhs = cellfun(@(a) a.mat, rhs, 'uniform', 0);
% $$$             res = vertcat(lhs.mat, rhs{:});
% $$$             res = replab.cyclotomic(res);
% $$$         end

    end

end
