classdef cyclotomic
% Multidimensional array with elements in the cyclotomic field
%
% Requires Java and the cyclolab support library.
%
% The implementation is pretty slow, but values are exact.
%
% We can construct cyclotomic arrays in various ways.
%
% They can be constructed from double floating-point arrays. Note that only fractions with a power-of-two
% denominator can be represented exactly.
%
% Example:
%   >>> replab.cyclotomic(1/2)
%       1/2
%   >>> replab.cyclotomic(1/3)
%       6004799503160661/18014398509481984
%
% Note that passing a cyclotomic array as a parameter is the identity operation. This is useful when converting
% number types.
%
% Example:
%   >>> c = replab.cyclotomic(1)
%       1
%   >>> replab.cyclotomic(c)
%       1
%
% Cyclotomics can also be constructed from vpis:
%
% Example:
%   >>> c = replab.cyclotomic(vpi('100000000000000'))
%       100000000000000
%
% They can also be constructed from strings. A vector/matrix syntax is also available, but row coefficients
% must always be separated by commas.
%
% Example:
%   >>> replab.cyclotomic('2/3')
%       2/3
%   >>> replab.cyclotomic('[1, 0; 0, 1]')
%       1  0
%       0  1
%
% The string syntax accepts:
%
% - integers,
% - the four operations ``+``, ``-``, ``*``, ``/``,
% - powers ``^`` with integer exponents,
% - parentheses,
% - roots of unity as in Gap System (``E(n)`` is ``exp(2*i*pi/n)``),
% - square roots of rational numbers,
% - sines and cosines of rational multiples of ``pi``.
%
% Example:
%   >>> replab.cyclotomic('sqrt(2)')
%       E(8)-E(8)^3
%   >>> replab.cyclotomic('cos(-pi/3)')
%       -1/2
%   >>> replab.cyclotomic('(3 + 2/3)^3')
%       1331/27
%   >>> replab.cyclotomic('E(3)^2')
%       E(3)^2
%
% Finally, the constructor also accepts heterogenous cell arrays, where the types above can be mixed and matched
% (except that strings must represent scalars, not matrices/vectors).
%
% Example:
%   >>> replab.cyclotomic({'1' vpi(2); '1/2' 1})
%         1    2
%        1/2   1
%
% A variety of static methods is available.
%
% Example:
%   >>> I = replab.cyclotomic.eye(3)
%       1  0  0
%       0  1  0
%       0  0  1

    properties (Access = protected)
        size_ % (integer(1,\*)): Shape
        data_ % (``cyclo.Cyclo[]``): 1D array of cyclotomic numbers, elements are stored as in ``array(:)``
    end

    methods (Static, Access = protected)

        function [lhs1, rhs1] = shape(lhs, rhs)
            lhs1 = replab.cyclotomic(lhs);
            rhs1 = replab.cyclotomic(rhs);
            if isscalar(lhs1)
                if ~isscalar(rhs1)
                    lhs1 = replab.cyclotomic.broadcast(lhs1, size(rhs1));
                end
            else
                if isscalar(rhs1)
                    rhs1 = replab.cyclotomic.broadcast(rhs1, size(lhs1));
                end
            end
            assert(isequal(size(lhs1), size(rhs1)));
        end

        function res = broadcast(array, shape)
        % Broadcast the given array to the given shape
            array = replab.cyclotomic(array);
            if length(array) == 1
                res = replab.cyclotomic(repmat({array}, shape));
            elseif isequal(size(array), shape)
                res = array;
            else
                error('Broadcast not implemented');
            end
        end

    end

    methods (Static) % Standard MATLAB methods

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
            c = c.data;
            for i = 1:length(I)
                c(I(i) + (J(i) - 1)*nR) = K.cyclo(i);
            end
            c = replab.cyclotomic(c, [nR nC]);
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
            if nargin < 1
                n = 1;
            end
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
            if nargin < 1
                varargin = {1};
            end
            array = repmat({javaMethod('zero', 'cyclo.Cyclo')}, varargin{:});
            c = replab.cyclotomic(array);
        end

        function c = ones(varargin)
        % Constructs a cyclotomic array filled with ones
        %
        % Returns:
        %    `.cyclotomic`: The constructed matrix
            if nargin < 1
                varargin = {1};
            end
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
            ja = javaMethod('E', 'cyclo.Lab', orders(:));
            c = replab.cyclotomic(ja, size(orders));
        end

        function M = rand(varargin)
        % Random cyclotomic matrix
        %
        % Keyword Args:
        %   maximalNumberOfTerms (integer): Maximal number of terms
        %   field ('R', 'C', 'Q'): Whether to sample from real or complex cyclotomics
        %   maximalOrder (integer): Maximal root-of-unity order
        %   maximalDenominator (integer): Largest denominator for rational numbers
        %   maximalAbsoluteValue (integer): Approximate largest absolute value
            args = struct('maximalNumberOfTerms', 5, 'field', 'C', 'maximalOrder', 12, 'maximalDenominator', 12, 'maximalAbsoluteValue', 5);
            mask = cellfun(@(x) isa(x, 'double'), varargin);
            startKWargs = find(~mask);
            if ~isempty(startKWargs)
                args = replab.util.populateStruct(args, varargin(startKWargs:end));
                dims = varargin(1:startKWargs-1);
            else
                dims = varargin;
            end
            if length(dims) == 0
                dims = [1 1];
            elseif length(dims) == 1
                dims = dims{1};
                if length(dims) == 1
                    dims = [dims dims];
                end
            else
                assert(all(cellfun(@length, dims) == 1))
                dims = cell2mat(dims);
            end
            if args.field == 'Q'
                args.field = 'R';
                args.maximalOrder = 1;
            end
            M = replab.cyclotomic.zeros(dims);
            maximalNumerator = args.maximalAbsoluteValue * args.maximalDenominator;
            for i = 1:args.maximalNumberOfTerms
                order = randi([1 args.maximalOrder]);
                exponent = randi([0 order]);
                root = replab.cyclotomic.E(order)^exponent;
                numerator = randi([-maximalNumerator maximalNumerator], dims);
                numerator = numerator .* (randi([1 args.maximalNumberOfTerms], dims) <= 2);
                numerator = replab.cyclotomic(numerator);
                denominator = randi([1 args.maximalDenominator]);
                M = M + (root*numerator)/denominator;
            end
            if args.field == 'R'
                M = (M + conj(M))/2;
            end
        end

    end

    methods (Static) % Cyclotomic constructions

        function c = approximate(lowerBounds, upperBounds)
        % Constructs the simplest rational approximation within an interval, for coefficients of a matrix
        %
        % Example:
        %   >>> replab.cyclotomic.approximate(3.141592, 3.141593)
        %     355/113
        %
        % Args:
        %   lowerBounds (double(...)): Interval lower bounds
        %   upperBounds (double(...)): Interval upper bounds
        %
        % Returns:
        %   `.cyclotomic`: The simplest rational approximation, coefficient by coefficient
            assert(isequal(size(lowerBounds), size(upperBounds)));
            c = replab.cyclotomic(javaMethod('approximate', 'cyclo.Lab', lowerBounds(:), upperBounds(:)), size(lowerBounds));
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
            c = replab.cyclotomic(javaMethod('fromRational', 'cyclo.Lab', numerators(:), denominators(:)), size(numerators));
        end

        function c = sqrtRational(num, den)
        % Returns the square root of a rational number
        %
        % Example:
        %   >>> replab.cyclotomic.sqrtRational(2)
        %       E(8)-E(8)^3
        %   >>> replab.cyclotomic.sqrtRational(1,2)
        %       E(8)/2-E(8)^3/2
        %
        % Args:
        %    num (integer): Numerator
        %    den (integer, optional): Denominator, default value ``1``
        %
        % Returns:
        %    `.cyclotomic`: The value ``sqrt(num/den)``
            if nargin < 2
                den = ones(size(num));
            end
            assert(isequal(size(num), size(den)));
            c = sqrt(replab.cyclotomic.fromRationals(num, den));
        end

    end

    methods (Static, Access = protected) % Private conversion methods

        function c = convertVpi(vec)
            vec = vec(:);
            s = arrayfun(@(i) strtrim(num2str(vec(i))), 1:length(vec), 'uniform', 0);
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
            charMask = cellfun(@ischar, vec);
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
                a = javaMethod('fromDouble', 'cyclo.Lab', [vec{ind}]);
                for i = 1:length(ind)
                    c(ind(i)) = a(i);
                end
            end
            ind = find(complexMask);
            if ~isempty(ind)
                a = javaMethod('fromDouble', 'cyclo.Lab', real([vec{ind}]), imag([vec{ind}]));
                for i = 1:length(ind)
                    c(ind(i)) = a(i);
                end
            end
            ind = find(charMask);
            if ~isempty(ind)
                a = javaMethod('parse', 'cyclo.Lab', vec(ind));
                for i = 1:length(ind)
                    c(ind(i)) = a(i);
                end
            end
            ind = find(vpiMask);
            if ~isempty(ind)
                a = javaMethod('parse', 'cyclo.Lab', cellfun(@(v) strtrim(num2str(v)), vec(ind), 'uniform', 0));
                for i = 1:length(ind)
                    c(ind(i)) = a(i);
                end
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
        %
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

    end

    methods

        function self = cyclotomic(array, size_)
        % Constructs a cyclotomic array from an array of coefficients
        %
        % Args:
        %   array: Coefficients
            if isa(array, 'char')
                res = javaMethod('parseMatrix', 'cyclo.Lab', array);
                assert(~res.isError, 'Parse error at character %d in ''%s''', res.errorIndex + 1, array);
                d = javaMethod('data', res);
                s = [javaMethod('nRows', res) javaMethod('nCols', res)];
                self.data_ = d;
                self.size_ = s;
            elseif isa(array, 'replab.cyclotomic')
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
            for i = 1:length(self.data_)
                assert(isa(self.data_(i), 'cyclo.Cyclo'));
            end
        end

        function c = cyclo(self, varargin)
        % Returns the cyclo.Cyclo object at the given index in this multidimensional array
        %
        % Arguments are integer indices, as passed to subsref
            l = length(varargin);
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

        function display(self, name)
            if nargin < 2
                name = 'ans';
            end
            n = ndims(self);
            if n > 2
                d = size(self);
                sub = cell(1, ndims(self) - 2);
                for ind = 1:prod(d(3:end))
                    [sub{:}] = ind2sub(d(3:end), ind);
                    s = struct('type', {'()'}, 'subs', {horzcat({':' ':'}, sub)});
                    disp([name '(:,:,' strjoin(cellfun(@num2str, sub, 'uniform', 0), ',') ') =']);
                    slice = self.subsref(s);
                    t = replab.compat.javaArrayToCell(javaMethod('print', 'cyclo.Lab', slice.data_));
                    t = replab.str.Table(reshape(t, size(slice)));
                    disp(t);
                end
            else
                disp([name ' =']);
                t = replab.compat.javaArrayToCell(javaMethod('print', 'cyclo.Lab', self.data_));
                t = replab.str.Table(reshape(t, size(self)));
                disp(t);
            end
        end

        function disp(self)
        % Standard display method
            self.display('ans');
        end

        function s = num2str(self)
        % Conversion to string (incomplete)
        %
        % Only supports matrices or vectors.
        %
        % Does not support format specifiers or a specified precision.
        %
        % Returns:
        %   charstring: Possibly multiline string representation of the matrix
            t = replab.compat.javaArrayToCell(javaMethod('print', 'cyclo.Lab', self.data_));
            if self.numel == 1
                s = t{1};
            else
                t = replab.str.Table(reshape(t, size(self)), 'uniform', 0);
                s = t.format(1000, 1000);
            end
        end

    end

    methods % Standard MATLAB properties

        function res = imag(self)
        % Returns the imaginary part
            res = replab.cyclotomic(-1i)*(self - conj(self))/2;
        end

        function res = isempty(self)
            res = prod(self.size_) == 0;
        end

        function res = isrational(self)
        % Returns which coefficients are rational
            res = reshape(javaMethod('isRational', 'cyclo.Lab', self.data_), self.size_);
        end

        function res = isreal(self)
            res = all(all(self == conj(self)));
        end

        function res = isscalar(self)
            res = length(self) == 1;
        end

        function res = isvector(self)
            res = sum(self.size_ ~= 1) <= 1;
        end

        function res = iswhole(self)
        % Returns which coefficients are integers
            res = reshape(javaMethod('isWhole', 'cyclo.Lab', self.data_), self.size_);
        end

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
                l = max(self.size_);
            end
        end

        function n = ndims(self)
        % Number of dimensions
        %
        % Returns:
        %   integer: Number of dimensions
            n = length(self.size_);
        end

        function n = numel(self)
        % Number of elements
        %
        % Returns:
        %   integer: Number of elements
            n = prod(self.size_);
        end

        function n = numArgumentsFromSubscript(self, s, indexingContext)
            n = 1;
        end

        function res = real(self)
        % Returns the real part
            res = (self + conj(self))/2;
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

        function res = dot(lhs, rhs)
        % Vector dot product (incomplete)
        %
        % Only supports the case where both ``lhs`` and ``rhs`` are vectors, and does not accept a ``dim`` argument.
            assert(isvector(lhs) && isvector(rhs));
            lhs = replab.cyclotomic(lhs);
            rhs = replab.cyclotomic(rhs);
            res = sum(conj(lhs(:)).*rhs(:));
        end

        function res = eq(lhs, rhs)
        % Equality test
            [lhs, rhs] = replab.cyclotomic.shape(lhs, rhs);
            res = reshape(javaMethod('eqv', 'cyclo.Lab', lhs.data, rhs.data), size(lhs));
        end

        function res = galois(self, ord)
        % Action of the Galois group
            res = replab.cyclotomic(javaMethod('galois', 'cyclo.Lab', self.data, ord), size(self));
        end

        function res = kron(lhs, rhs)
            lhs = replab.cyclotomic(lhs);
            rhs = replab.cyclotomic(rhs);
            res = replab.cyclotomic(javaMethod('kron', 'cyclo.Lab', size(lhs, 1), size(lhs, 2), size(rhs, 1), size(rhs, 2), lhs.data, rhs.data), [size(lhs, 1)*size(rhs, 1) size(lhs, 2)*size(rhs, 2)]);
        end

        function res = minus(lhs, rhs)
        % Standard ``-`` operator
            [lhs, rhs] = replab.cyclotomic.shape(lhs, rhs);
            res = replab.cyclotomic(javaMethod('minus', 'cyclo.Lab', lhs.data, rhs.data), size(lhs));
        end

        function res = mpower(self, e)
        % Matrix power
            assert(ndims(self) == 2);
            n = size(self, 1);
            assert(size(self, 2) == n);
            res = replab.cyclotomic(javaMethod('power', 'cyclo.Lab', n, self.data, e), [n n]);
        end

        function res = mrdivide(lhs, rhs)
        % Division
        %
        % Only supports the case where the right hand side is scalar
            lhs = replab.cyclotomic(lhs);
            rhs = replab.cyclotomic(rhs);
            assert(isscalar(rhs), 'The / operator is only implemented for scalar rhs');
            rd = rhs.data;
            res = replab.cyclotomic(javaMethod('divideScalar', 'cyclo.Lab', lhs.data, rd(1)), size(lhs));
        end

        function res = mtimes(lhs, rhs)
        % Matrix multiplication
        %
        % Support both ``m x n`` by ``n x p`` matrix multiplication, and the case where one of the arguments is a scalar
            lhs = replab.cyclotomic(lhs);
            rhs = replab.cyclotomic(rhs);
            if ~isscalar(rhs) && isscalar(lhs)
                res = rhs * lhs;
                return
            end
            if isscalar(rhs)
                l = size(lhs, 1);
                m = size(lhs, 2);
                rval = rhs.cyclo(1);
                res = replab.cyclotomic(javaMethod('timesScalar', 'cyclo.Lab', lhs.data_, rval), size(lhs));
                return
            end
            l = size(lhs, 1);
            m = size(lhs, 2);
            assert(m == size(rhs, 1));
            n = size(rhs, 2);
            res = replab.cyclotomic(javaMethod('times', 'cyclo.Lab', l, m, n, lhs.data, rhs.data), [l n]);
        end

        function res = ne(self, rhs)
        % (Non-)equality test
            res = ~(self == rhs);
        end

        function res = plus(lhs, rhs)
        % Standard ``+`` operator
            [lhs, rhs] = replab.cyclotomic.shape(lhs, rhs);
            res = replab.cyclotomic(javaMethod('plus', 'cyclo.Lab', lhs.data, rhs.data), size(lhs));
        end

        function res = rdivide(lhs0, rhs0)
        % Pointwise ``/`` operator
            [lhs, rhs] = replab.cyclotomic.shape(lhs0, rhs0);
            res = replab.cyclotomic(javaMethod('pw_divide', 'cyclo.Lab', lhs.data, rhs.data), size(lhs));
        end

        function res = times(lhs0, rhs0)
        % Pointwise ``*`` operator
            [lhs, rhs] = replab.cyclotomic.shape(lhs0, rhs0);
            res = replab.cyclotomic(javaMethod('pw_times', 'cyclo.Lab', lhs.data, rhs.data), size(lhs));
        end

    end

    methods % Shape

        function d = diag(self, varargin)
        % Diagonal matrices and diagonals of a matrix
        %
        % Fully compatible with MATLAB's diag
            ind = reshape(1:numel(self), size(self));
            ind = diag(ind, varargin{:});
            c = javaArray('cyclo.Cyclo', numel(ind));
            z = javaMethod('zero', 'cyclo.Cyclo');
            for i = 1:numel(ind)
                if ind(i) == 0
                    c(i) = z;
                else
                    c(i) = self.data_(ind(i));
                end
            end
            d = replab.cyclotomic(c, size(ind));
        end

        function res = reshape(self, varargin)
        % Reshape array
        %
        % Fully compatible with MATLAB's reshape
            ind = reshape(1:numel(self), size(self));
            ind = reshape(ind, varargin{:});
            res = replab.cyclotomic(self.data_, size(ind));
        end

        function res = permute(self, order)
        % Permute array dimensions
        %
        % Fully compatible with MATLAB's permute
            ind = reshape(1:numel(self), size(self));
            ind = permute(ind, order);
            c = self.data_;
            c = c(ind(:));
            res = replab.cyclotomic(c, size(ind));
        end

        function res = ipermute(self, order)
        % Inverse permute array dimensions
        %
        % Fully compatible with MATLAB's permute
            ind = reshape(1:numel(self), size(self));
            ind = ipermute(ind, order);
            c = self.data_;
            c = c(ind(:));
            res = replab.cyclotomic(c, size(ind));
        end

    end

    methods % Involutions

        function res = conj(self)
        % Complex conjugation
        %
        % Fully compatible with MATLAB's conj
            res = replab.cyclotomic(javaMethod('conjugate', 'cyclo.Lab', self.data), self.size);
        end

        function res = ctranspose(self)
        % Complex conjugate transpose
        %
        % Fully compatible with MATLAB's ctranspose
            res = conj(transpose(self));
        end

        function res = inv(self)
        % Matrix inverse
        %
        % Requires that the matrix is full rank
            assert(ndims(self) == 2);
            n = size(self, 1);
            assert(size(self, 2) == n);
            res = replab.cyclotomic(javaMethod('inverse', 'cyclo.Lab', n, self.data_), [n n]);
        end

        function res = transpose(self)
        % Transpose
        %
        % Fully compatible with MATLAB's transpose
            assert(ndims(self) == 2);
            ind = reshape(1:numel(self), size(self));
            ind = ind';
            ind = ind(:);
            c = self.data;
            c = c(ind);
            s = [size(self, 2) size(self, 1)];
            res = replab.cyclotomic(c, s);
        end

        function res = uminus(self)
        % Unary minus
            res = replab.cyclotomic(javaMethod('negate', 'cyclo.Lab', self.data), size(self));
        end

    end

    methods % Matrix decompositions

        function res = null(self)
        % Computes the null space of a cyclotomic matrix
        %
        % Example:
        %   >>> M = replab.cyclotomic([1 3 0; -2 -6 0; 3 9 6]);
        %   >>> null(M)'
        %       -3  1  0
        %
        % Returns:
        %   `.cyclotomic`: The matrix null space
            rr = javaMethod('rref', 'cyclo.Lab', self.data, size(self, 1), size(self, 2));
            rank = double(javaMethod('rank', rr));
            rows = size(self, 2);
            cols = size(self, 2) - rank;
            res = replab.cyclotomic(javaMethod('nullSpace', rr), [rows cols]);
        end

        function res = rref(self)
        % Computes the row reduced echelon form
        %
        % Returns:
        %   `.cyclotomic`: The row reduced echelon form
            rr = javaMethod('rref', 'cyclo.Lab', self.data, size(self, 1), size(self, 2));
            res = replab.cyclotomic(javaMethod('matrix', rr), size(self));
        end

        function [L, U, p] = lu(self)
        % Returns the LU decomposition of this matrix
        %
        % The ``m x n`` matrix must have ``m >= n``.
        %
        % It returns triangular matrices ``L``, ``U`` and a permutation vector ``p`` such that ``L*U == self(:,p)``
        %
        % Returns
        % -------
        %   L: `.cyclotomic`(\*,\*)
        %     Lower triangular matrix
        %   U: `.cyclotomic`(\*,\*)
        %     Upper triangular matrix
        %   p: integer(1,\*)
        %     Permutation vector
            m = size(self, 1);
            n = size(self, 2);
            assert(m >= n, 'LU decomposition only works for m x n matrices with m >= n');
            res = javaMethod('lu', 'cyclo.Lab', self.data, m, n);
            L = javaMethod('L', res);
            U = javaMethod('U', res);
            p = javaMethod('pivots', res);
            p = double(p(:)') + 1;
            L = replab.cyclotomic(L, [m n]);
            U = replab.cyclotomic(U, [n n]);
        end

        function [L, D, p] = ldl(self)
            m = size(self, 1);
            n = size(self, 2);
            assert(m == n);
            res = javaMethod('ldl', 'cyclo.Lab', self.data, m, n);
            L = javaMethod('L', res);
            D = javaMethod('D', res);
            p = javaMethod('pivots', res);
            p = double(p(:)') + 1;
            L = replab.cyclotomic(L, [n n]);
            D = replab.cyclotomic(D, [n n]);
        end

    end

    methods % Indexing

        function varargout = find(self, varargin)
        % Find indives of nonzero elements
        %
        % Compatible with MATLAB's find
            mask = self ~= 0;
            if nargout == 3
                [I, J] = find(mask, varargin{:});
                c = javaArray('cyclo.Cyclo', length(I));
                nr = size(self, 1);
                for i = 1:length(I)
                    c(i) = self.data_(I(i) + (J(i)-1)*nr);
                end
                V = replab.cyclotomic(c, size(I));
                varargout{1} = I;
                varargout{2} = J;
                varargout{3} = V;
            else
                [varargout{:}] = find(mask, varargin{:});
            end
        end

        function varargout = subsref(self, s)
        % Indexing
        %
        % Compatible with MATLAB's subsref
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
        %
        % Compatible with MATLAB's subsasgn
            switch s(1).type
              case '.'
                self = builtin('subsasgn', self, s, val);
              case '()'
                if length(s) == 1
                    ind = reshape(1:prod(size(self)), size(self));
                    ind = subsref(ind, s(1));
                    c = data(self);
                    val = replab.cyclotomic(val);
                    val = replab.cyclotomic.broadcast(val, size(ind));
                    ind = ind(:);
                    val = data(val);
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

    methods % Element-by-element computations

        function res = double(self)
        % Conversion to floating-point double
            pairs = javaMethod('toDouble', 'cyclo.Lab', self.data_);
            if replab.compat.isOctave
                R = pairs(1);
                I = pairs(2);
            else
                R = pairs(1,:);
                I = pairs(2,:);
            end
            res = reshape(R, self.size_) + 1i * reshape(I, self.size_);
        end

        function [res, error] = doubleApproximation(self)
        % Conversion to floating-point double with estimation of the error
            data = javaMethod('toDoubleApproximation', 'cyclo.Lab', self.data_);
            res = reshape(data.real, self.size_) + 1i * reshape(data.imag, self.size_);
            error = reshape(data.error, self.size_);
        end

        function h = hash(self)
        % Returns a multidimensional array of hash codes
        %
        % Example:
        %   >>> hash(replab.cyclotomic.zeros(1)) == -1579025880
        %       1
        %
        % Returns:
        %   integer(...): Integer array of hash codes
            h = arrayfun(@(c) javaMethod('hashCode', c), self.data_);
            h = reshape(h, self.size_);
        end

        function res = sqrt(self)
        % Square root
        %
        % Requires that the argument contains rational coefficients, otherwise compatible with MATLAB's sqrt
            res = replab.cyclotomic(javaMethod('sqrt', 'cyclo.Lab', self.data_), self.size_);
        end

    end

    methods % Contractions

        function p = prod(self, dim)
        % Product of elements (incomplete)
        %
        % Does not support arrays with more than two dimensions, and does not support
        % extra options (type, nanflag).
            if isempty(self)
                p = replab.cyclotomic.ones;
            elseif isscalar(self)
                p = self;
            else
                if nargin < 2 || isempty(dim)
                    dim = find(self.size_ ~= 1, 1); % Find first non-singleton dimension
                elseif isequal(dim, 'all')
                    dim = [1 2];
                end
                if isequal(dim, [1 2]) || isequal(dim, [2 1])
                    p = prod(self(:));
                elseif dim == 2
                    p = prod(self.').';
                elseif dim == 1
                    shift = 0;
                    data = javaArray('cyclo.Cyclo', size(self, 2));
                    for i = 1:size(self, 2)
                        data(i) = javaMethod('prod', 'cyclo.Lab', self.data_(shift+(1:size(self,1))));
                        shift = shift + size(self, 1);
                    end
                    p = replab.cyclotomic(data, [1 size(self, 2)]);
                else
                    error('Invalid dim argument');
                end
            end
        end

        function s = sum(self, dim)
        % Sum of elements (incomplete)
        %
        % Does not support arrays with more than two dimensions, and does not support
        % extra options (type, nanflag).
            if isempty(self)
                s = replab.cyclotomic.zeros;
            elseif isscalar(self)
                s = self;
            else
                if nargin < 2 || isempty(dim)
                    dim = find(self.size_ ~= 1, 1); % Find first non-singleton dimension
                elseif isequal(dim, 'all')
                    dim = [1 2];
                end
                if isequal(dim, [1 2]) || isequal(dim, [2 1])
                    s = sum(self(:));
                elseif dim == 2
                    s = sum(self.').';
                elseif dim == 1
                    shift = 0;
                    data = javaArray('cyclo.Cyclo', size(self, 2));
                    for i = 1:size(self, 2)
                        data(i) = javaMethod('sum', 'cyclo.Lab', self.data_(shift+(1:size(self,1))));
                        shift = shift + size(self, 1);
                    end
                    s = replab.cyclotomic(data, [1 size(self, 2)]);
                else
                    error('Invalid dim argument');
                end
            end
        end

        function v = trace(self)
        % Sum of diagonal elements
        %
        % Fully compatible with MATLAB's trace
            v = sum(diag(self));
        end

    end

    methods % Matrix concatenation

        function res = blkdiag(lhs, varargin)
        % Block diagonal concatenation of matrix input arguments
        %
        % Compatible with MATLAB's blkdiag
            lhs = replab.cyclotomic(lhs);
            rhs = cellfun(@(r) replab.cyclotomic(r), varargin, 'uniform', 0);
            assert(ndims(lhs) == 2);
            assert(all(cellfun(@(r) ndims(r) == 2, varargin)));
            dim1 = [size(lhs, 1) cellfun(@(r) size(r, 1), rhs)];
            dim2 = [size(lhs, 2) cellfun(@(r) size(r, 2), rhs)];
            res = replab.cyclotomic.zeros(sum(dim1), sum(dim2));
            ind = reshape(1:numel(res), size(res));
            res = res.data;
            I = ind(1:dim1(1), 1:dim2(1));
            data = lhs.data;
            for i = 1:length(I)
                res(I(i)) = data(i);
            end
            s1 = dim1(1);
            s2 = dim2(1);
            for i = 1:length(rhs)
                I = ind(s1+(1:dim1(i+1)), s2+(1:dim2(i+1)));
                data = rhs{i}.data;
                for j = 1:length(I)
                    res(I(i)) = data(j);
                end
                s1 = s1 + dim1(i+1);
                s2 = s2 + dim2(i+1);
            end
            res = replab.cyclotomic(res, size(ind));
        end

        function res = horzcat(lhs, varargin)
        % Horizontal concatenation (incomplete)
        %
        % Incomplete implementation: only works on vectors and matrices
            assert(ndims(lhs) == 2 && all(cellfun(@(a) ndims(a) == 2, varargin)), 'horzcat only works on matrices');
            dims = [size(lhs, 2) cellfun(@(a) size(a, 2), varargin)];
            lhs = replab.cyclotomic(lhs);
            rhs = cellfun(@(a) replab.cyclotomic(a), varargin, 'uniform', 0);
            rhs = cellfun(@(a) a.data, rhs, 'uniform', 0);
            m = size(lhs, 1);
            res = lhs.data;
            for i = 1:length(rhs)
                res = javaMethod('cat', 'cyclo.Lab', res, rhs{i});
            end
            res = replab.cyclotomic(res, [m sum(dims)]);
        end

        function res = vertcat(lhs, varargin)
        % Vertical concatenation
            lhs = replab.cyclotomic(lhs).';
            rhs = cellfun(@(a) replab.cyclotomic(a).', varargin, 'uniform', 0);
            res = [lhs rhs{:}];
            res = res.';
        end

    end

end
