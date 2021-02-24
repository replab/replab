classdef H
% Quaternion matrix type
%
% A quaternion matrix can be expressed using four real matrices A, B, C and D:
%
% $ Q = A + B i + C j + D k $ where $i^2 = j^2 = k^2 = i j k = -1$ and $i j = k$, $j k = i$ and $k i = j$.
%
% We choose to encode our quaternions using two matrices $X = A + B i$ and $Y = C + D i$, where
% we now identify $i$ with the complex imaginary unit.
%
% Thus $ Q = X + Y j $.
%
% The class supports sparse components for efficiency.

    properties (SetAccess = protected)
        X % First part
        Y % Second part
    end

    methods (Static)

        function res = i
        % Quaternion unit ``i``
        %
        % Returns:
        %   `.H`: Quaternion scalar matrix
            res = replab.H(0, 1);
        end

        function res = j
        % Quaternion unit ``j``
        %
        % Returns:
        %   `.H`: Quaternion scalar matrix
            res = replab.H(0, 0, 1);
        end

        function res = k
        % Quaternion unit ``k``
        %
        % Returns:
        %   `.H`: Quaternion scalar matrix
            res = replab.H(0,0,0,1);
        end

        function res = randn(varargin)
        % Returns a normally distributed quaternion multidimensional array
        %
        % Args:
        %   varargin: Arguments passed to ``randn``
        %
        % Returns:
        %   `.H`: Quaternion multidimensional array
            res = replab.H(randn(varargin{:}), randn(varargin{:}), randn(varargin{:}), randn(varargin{:}));
        end

        function [X, Y] = decompose(Q)
        % Decomposes the first and second part from the given matrix, interpreted as a quaternion matrix
        %
        % The result of this call is equivalent to ``tmp = replab.H(Q); X = tmp.X; Y = tmp.Y``, but it
        % avoids unnecessary construction of temporary objects.
        %
        % Args:
        %   Q (`replab.H` or double multidimensional array): Quaternion to decompose
        %
        % Returns
        % -------
        %   X: double multidimensional array
        %     First part of the quaternion matrix
        %   Y: double multidimensional array
        %     Second part of the quaternion matrix
            if isa(Q, 'replab.H')
                X = Q.X;
                Y = Q.Y;
            elseif isa(Q, 'double')
                X = Q;
                if ndims(Q) <= 2
                    Y = sparse(size(Q, 1), size(Q, 2));
                else
                    Y = zeros(size(Q));
                end
            else
                error('Unsupported type %s', class(Q));
            end
        end

        function Q = decode(Z)
        % Decodes the given complex matrix into a quaternion matrix
        %
        % Args:
        %   Z (double(\*,\*)): Complex matrix encoding a quaternion algebra, see `.encode`
        %
        % Returns:
        %   `+replab.H`: Quaternion matrix
            X = (Z(1:2:end,1:2:end)+conj(Z(2:2:end,2:2:end)))/2;
            Y = (Z(1:2:end,2:2:end)-conj(Z(2:2:end,1:2:end)))/2;
            Q = replab.H(X, [], Y);
        end

        function Z = encode(Q)
        % Encodes the given matrix, interpreted as a quaternion, into a complex matrix
        %
        % Uses the encoding of a quaternion ``q = a + b i + c j + d k`` as ``[a + i*b, c + i*d; -c + i*d, a - i*b]``
        %
        % Args:
        %   Q (`.H`): Quaternion matrix
        %
        % Returns:
        %   double(\*,\*): Complex matrix encoding the given quaternion matrix
            [X, Y] = replab.H.decompose(Q);
            Z = kron(real(X), [1 0; 0 1]) + kron(imag(X), [1i 0; 0 -1i]) + ...
                kron(real(Y), [0 1; -1 0]) + kron(imag(Y), [0 1i; 1i 0]);
        end

    end

    methods

        function self = H(M1, Mi, Mj, Mk)
        % Quaternion constructor
        %
        % Constructs a quaternion as:
        % `` Q = M1 + Mi * i + Mj * j + Mk * k `` where ``i``, ``j``, ``k`` are the quaternion units.
        %
        % Some matrices can be complex, in which case the complex imaginary unit is identified with ``i``.
        %
        % Args:
        %   M1 (double(d1,d2)): Must be provided
        %   Mi (double(d1,d2) or ``[]``): Matrix to be multiplied by ``i``
        %   Mj (double(d1,d2) or ``[]``): Matrix to be multiplied by ``j``
        %   Mk (double(d1,d2) or ``[]``): Matrix to be multiplied by ``k``

            if isa(M1, 'replab.H')
                self.X = M1.X;
                self.Y = M1.Y;
                return
            end

            if nargin < 3 || isempty(Mj)
                Mj = 0 * M1;
            end

            % first we move everything into M1 and Mj
            if nargin >= 2 && ~isempty(Mi)
                M1 = M1 + Mi * 1i;
            end

            if nargin == 4 && ~isempty(Mk)
                Mj = Mj + real(Mk) * 1i; % i j = k
                if ~isreal(Mk)
                    Mj = Mj - imag(Mk); % i k = - k i = -j
                end
            end

            assert(isequal(size(M1), size(Mj)));

            self.X = M1;
            self.Y = Mj;
        end

        function res = part1(self)
        % Returns the real part of the quaternion
        %
        % Returns:
        %   double(...): Real part
            res = real(self.X);
        end

        function res = parti(self)
        % Returns the ``i`` part of the quaternion
        %
        % Returns:
        %   double(...): Part
            res = imag(self.X);
        end

        function res = partj(self)
        % Returns the ``j`` part of the quaternion
        %
        % Returns:
        %   double(...): Part
            res = real(self.Y);
        end

        function res = partk(self)
        % Returns the ``k`` part of the quaternion
        %
        % Returns:
        %   double(...): Part
            res = imag(self.Y);
        end

        function C = toStringCell(self)
        % Computes a string representation of the elements of the array
        %
        % Returns:
        %   cell(...) of charstring: String representation
            X = self.X(:);
            Y = self.Y(:);
            C = cell(length(X), 1);
            for i = 1:length(C)
                x = X(i);
                y = Y(i);
                suffix = {'' 'i' 'j' 'k'};
                parts = full([real(x) imag(x) real(y) imag(y)]);
                s = '';
                plus = '';
                minus = '-';
                for j = 1:4
                    if parts(j) ~= 0
                        if parts(j) > 0
                            s = sprintf('%s%s%.5g%s', s, plus, parts(j), suffix{j});
                        else
                            s = sprintf('%s%s%.5g%s', s, minus, -parts(j), suffix{j});
                        end
                        plus = ' + ';
                        minus = ' - ';
                    end
                end
                if isempty(s)
                    s = '0';
                end
                C{i} = s;
            end
            C = reshape(C, self.size);
        end

    end

    methods % Implementations

        function r1 = eps(self)
        % Computes the spacing of floating point numbers
            r1 = eps(abs(self));
        end

        function r1 = abs(self)
        % Computes the absolute value, array coefficient by coefficient
        %
        % Returns:
        %   double: Magnitude of each quaternion coefficient
            r1 = sqrt(abs2(self));
        end

        function r1 = abs2(self)
        % Computes the square of the absolute value, array coefficient by coefficient
        %
        % Returns:
        %   double: Square of the magnitude of each quaternion coefficient

            r1 = real(self.X).^2 + imag(self.X).^2 + real(self.Y).^2 + imag(self.Y).^2;
        end

        function res = blkdiag(varargin)
            n = length(varargin);
            X = cell(1, n);
            Y = cell(1, n);
            for i = 1:n
                [X{i}, Y{i}] = replab.H.decompose(varargin{i});
            end
            X = blkdiag(X{:});
            Y = blkdiag(Y{:});
            res = replab.H(X, [], Y);
        end

        function res = conj(self)
            res = replab.H(conj(self.X), [], -self.Y);
        end

        function res = ctranspose(self)
            res = self.conj.transpose;
        end

        function e = end(self, k, n)
            e = builtin('end', self.X, k, n);
        end

        function res = horzcat(varargin)
            [X, Y] = replab.H.decompose(varargin{1});
            for i = 2:nargin
                [Xi, Yi] = replab.H.decompose(varargin{i});
                X = horzcat(X, Xi);
                Y = horzcat(Y, Yi);
            end
            res = replab.H(X, [], Y);
        end

        function res = diag(self)
            [X, Y] = replab.H.decompose(self);
            res = replab.H(diag(X), [], diag(Y));
        end

        function disp(self)
        % Standard display method
            t = replab.str.Table(self.toStringCell);
            disp(t);
        end

        function b = eq(Q1, Q2)
            [X1, Y1] = replab.H.decompose(Q1);
            [X2, Y2] = replab.H.decompose(Q2);
            b = (X1 == X2) & (Y1 == Y2);
        end

        function res = kron(Q1, Q2)
            [X1, Y1] = replab.H.decompose(Q1);
            [X2, Y2] = replab.H.decompose(Q2);
            X3 = kron(X1, X2) - kron(Y1, conj(Y2));
            Y3 = kron(Y1, conj(X2)) + kron(X1, Y2);
            res = replab.H(X3, [], Y3);
        end

        function res = inv(self)
            X = self.X;
            Y = self.Y;
            if nnz(Y) == 0
                res = replab.H(inv(X));
            elseif isequal(self.size, [1 1])
                res = self.conj * (1 ./ self.abs2);
            else
                Z = replab.H.encode(self);
                Zi = inv(Z);
                res = replab.H.decode(Zi);
            end
        end

        function b = isscalar(self)
            b = prod(self.size) == 1;
        end

        function res = minus(Q1, Q2)
            [X1, Y1] = replab.H.decompose(Q1);
            [X2, Y2] = replab.H.decompose(Q2);
            res = replab.H(X1 - X2, [], Y1 - Y2);
        end

        function res = mtimes(Q1, Q2)
            [X1, Y1] = replab.H.decompose(Q1);
            [X2, Y2] = replab.H.decompose(Q2);
            X = X1 * X2 - Y1 * conj(Y2);
            Y = Y1 * conj(X2) + X1 * Y2;
            res = replab.H(X, [], Y);
        end

        function Q3 = mldivide(Q1, Q2)
        % Backslash or left matrix divide
            Z1 = replab.H.encode(Q1);
            Z2 = replab.H.encode(Q2);
            Z3 = Z1\Z2;
            Q3 = replab.H.decode(Z3);
        end

        function Q3 = mrdivide(Q1, Q2)
        % Slash or right matrix divide
            if isscalar(Q2)
                Q3 = Q1 ./ Q2;
                return
            end
            Z1 = replab.H.encode(Q1);
            Z2 = replab.H.encode(Q2);
            Z3 = Z1/Z2;
            Q3 = replab.H.decode(Z3);
        end


        function b = ne(Q1, Q2)
            b = ~(Q1 == Q2);
        end

        function n = norm(self, p)
            assert(length(self.size) == 2, 'Input must be 2-D.');
            if isvector(self.X)
                assert(nargin < 2 || isequal(p, 2) || isequal(p, 'fro'));
                n = self.abs2;
                n = sqrt(sum(n(:)));
            elseif ismatrix(self.X)
                if nargin < 2
                    p = 2;
                end
                switch p
                  case 'fro'
                    n = self.abs2;
                    n = sqrt(sum(n(:)));
                  case 1
                    n = max(sum(self.abs, 1));
                  case inf
                    n = max(sum(self.abs, 2));
                  otherwise
                    s = replab.shortStr(p);
                    error(sprintf('Norm type %s undefined', s));
                end
            end
        end

        function res = plus(Q1, Q2)
            [X1, Y1] = replab.H.decompose(Q1);
            [X2, Y2] = replab.H.decompose(Q2);
            res = replab.H(X1 + X2, [], Y1 + Y2);
        end

        function res = real(self)
            res = real(self.X);
        end

        function res = rdivide(Q1, Q2)
            [X, Y] = replab.H.decompose(Q1 .* conj(Q2));
            if isa(Q2, 'double')
                f = real(Q2).^2 + imag(Q2).^2;
            else
                f = abs2(Q2);
            end
            res = replab.H(X./f, [], Y./f);
        end

        function res = sign(self)
            res = self ./ abs(self);
        end

        function s = size(self, varargin)
            s = size(self.X, varargin{:});
        end

        function res = sqrt(self)
        % Quaternion element-by-element square root
        %
        % See `<https://math.stackexchange.com/questions/382431/square-roots-of-quaternions>_`
            T = self.part1;
            X = self.parti;
            Y = self.partj;
            Z = self.partk;
            T = T(:);
            X = X(:);
            Y = Y(:);
            Z = Z(:);
            for i = 1:length(T)
                % q = t + v u
                % where t, v are real and u is a unit quaternion vector
                % so that u^2 = -1
                t = T(i);
                u = [X(i) Y(i) Z(i)];
                v = norm(u);
                if v ~= 0
                    u = u/v;
                else
                    u = [1 0 0];
                end
                s = sqrt(t + 1i*v);
                t = real(s);
                v = imag(s);
                T(i) = t;
                X(i) = v*u(1);
                Y(i) = v*u(2);
                Z(i) = v*u(3);
            end
            T = reshape(T, self.size);
            X = reshape(X, self.size);
            Y = reshape(Y, self.size);
            Z = reshape(Z, self.size);
            res = replab.H(T, X, Y, Z);
        end

        function varargout = subsasgn(self, s, p)
            if length(s) == 1 && isequal(s.type, '()')
                [qX, qY] = replab.H.decompose(self);
                [pX, pY] = replab.H.decompose(p);
                qX(s.subs{:}) = pX;
                qY(s.subs{:}) = pY;
                varargout{1} = replab.H(qX, [], qY);
            else
                [varargout{1:nargout}] = builtin('subsasgn', self, s, p);
            end
        end

        function varargout = subsref(self, s)
            if length(s) == 1 && isequal(s.type, '()')
                X = self.X(s.subs{:});
                Y = self.Y(s.subs{:});
                varargout{1} = replab.H(X, [], Y);
            else
                [varargout{1:nargout}] = builtin('subsref', self, s);
            end
        end

        function res = times(Q1, Q2)
            [X1, Y1] = replab.H.decompose(Q1);
            [X2, Y2] = replab.H.decompose(Q2);
            X = X1 .* X2 - Y1 .* conj(Y2);
            Y = Y1 .* conj(X2) + X1 .* Y2;
            res = replab.H(X, [], Y);
        end

        function res = trace(self)
            [X, Y] = replab.H.decompose(self);
            res = replab.H(trace(X), [], trace(Y));
        end

        function res = transpose(self)
            res = replab.H(self.X.', [], self.Y.');
        end

        function res = uminus(self)
            res = replab.H(-self.X, [], -self.Y);
        end

        function res = uplus(self)
            res = self;
        end

        function res = vertcat(varargin)
            [X, Y] = replab.H.decompose(varargin{1});
            for i = 2:nargin
                [Xi, Yi] = replab.H.decompose(varargin{i});
                X = vertcat(X, Xi);
                Y = vertcat(Y, Yi);
            end
            res = replab.H(X, [], Y);
        end

    end

end
