classdef rational
% Represents a matrix with rational coefficients
    properties
        N = []; % integer part
        d = 1; % denominator
    end
    methods
        function R = rational(M, d)
            switch nargin
              case 1
                if isa(M, 'rational')
                    R.N = M.N;
                    R.d = M.d;
                    return
                end
                N = [];
                D = [];
                switch class(M)
                  case 'double'
                      [N, D] = rat(M);
                  case 'sym'
                    [N, D] = numden(M);
                    if any(round(N)~=N)
                        error('Not a rational matrix.');
                    end
                    if any(round(D)~=D)
                        error('Not a rational matrix.');
                    end
                    if any(round(N)./round(D) ~= M)
                        error('Not a rational matrix.');
                    end
                    N = double(N);
                    D = double(D);
                end
                d = rational.matlcm(D);
                R.N = N .* (d./D);
                R.d = d;
              case 2
                if isequal(size(d), [1 1])
                    if isa(M, 'rational')
                        R.N = M.N;
                        R.d = M.d * d;
                    else
                        R.N = M;
                        R.d = d;
                    end
                else
                    D = d;
                    d = rational.matlcm(D);
                    R.N = N .* (d./D);
                    R.d = d;
                end
            end
        end
        function S = sym(R)
            S = sym(R.N)/sym(R.d);
        end
        function d = den(R)
            d = R.d;
        end
        function N = num(R)
            N = R.N;
        end
        function R = pinv(R)
            R1 = R.N(any(R.N ~= 0, 2), any(R.N ~= 0, 1));
            isdiag = @(X) isequal(diag(diag(X)), X);
            if ~isdiag(R1)
                error(['Does not know how to pseudo-inverse nondiagonal ' ...
                       'matrices']);
            end
            INV = R.d./rational(R.N(R.N~=0));
            R.N(R.N~=0) = INV.N;
            R.N = R.N';
            R.d = INV.d;
        end
        function R1 = inv(R)
            isdiag = @(X) isequal(diag(diag(X)), X);
            if ~isdiag(R)
                error(['Does not know how to inverse nondiagonal ' ...
                       'matrices']);
            end
            R1 = diag(1./diag(R));
        end
        function R = simplify(R)
            Ngcd = rational.matgcd(R.N);
            g = gcd(Ngcd, R.d);
            if g ~= 1
                R.N = R.N ./ g;
                R.d = R.d / g;
            end
        end
        function M = double(R)
            M = R.N/R.d;
        end
        function R1 = diag(R)
            R1 = rational(diag(R.N), R.d);
        end
        function R = mtimes(A, B)
            if isa(A, 'sym') || isa(B, 'sym')
                R = sym(A) * sym(B);
            else
                A = rational(A);
                B = rational(B);
                R = simplify(rational(A.N * B.N, A.d * B.d));
            end
        end
        function R = times(A, B)
            if isa(A, 'sym') || isa(B, 'sym')
                R = sym(A) .* sym(B);
            else
                A = rational(A);
                B = rational(B);
                R = simplify(rational(A.N .* B.N, A.d * B.d));
            end
        end
        function [root, res] = sqrt(R)
        % rational square root. root.^2 * res == R
            R1 = simplify(R);
            droot = 1;
            dres = 1;
            d = R.d;
            f = unique(factor(d));
            for i = f
                while gcd(d, i*i) > 1
                    d = d / (i*i);
                    droot = droot * i;
                end
                while gcd(d, i) > 1
                    d = d / i;
                    dres = dres * i;
                end
            end
            N = R.N;
            Nroot = zeros(size(N));
            Nres = zeros(size(N));
            Nroot(N ~= 0) = 1;
            Nres(N > 0) = 1;
            Nres(N < 0) = -1;
            N = abs(N);
            f = unique(factor(abs(rational.matlcm(N))));
            for i = f
                while gcd(rational.matlcm(N), i*i) == i*i && i > 1
                    ind = (mod(N, i*i) == 0) & (N ~= 0);
                    N(ind) = N(ind) / (i*i);
                    Nroot(ind) = Nroot(ind) * i;
                end
                while gcd(rational.matlcm(N), i) == i && i > 1
                    ind = (mod(N, i) == 0) & (N ~= 0);
                    N(ind) = N(ind) / i;
                    Nres(ind) = Nres(ind) * i;
                end
            end
            root = rational(Nroot, droot);
            res = rational(Nres, dres);
        end
        function R = power(R1, n)
            R = rational(R1.N.^n, R1.d.^n);
        end
        function R = transpose(R)
            R.N = R.N';
        end
        function R = ctranspose(R)
            R.N = R.N';
        end
        function R = rref(R)
            R.N = rref(R.N);
            R.d = 1;
        end
        function R = permute(R, order)
            R.N = permute(R.N, order);
        end
        function R = squeeze(R)
            R.N = squeeze(R.N);
        end
        function c = char(R)
            assert(isequal(size(R), [1 1]));
            [N, D] = rat(R);
            if D ~= 1
                c = [num2str(N) '/' num2str(D)];
            else
                c = num2str(N);
            end
        end
        function R = abs(R)
            R.N = abs(R.N);
        end
        function R = uminus(R)
            R.N = -R.N;
        end
        function [N, D] = rat(R)
            N = R.N;
            D = ones(size(R.N)) * R.d;
            g = gcd(N, D);
            N = N./g;
            D = D./g;
        end
        function s = size(R, varargin)
            s = size(R.N, varargin{:});
        end
        function R1 = sum(R, varargin)
            R1 = rational(sum(R.N, varargin{:}), R.d);
        end
        function R = horzcat(A, B, varargin)
            if nargin == 1
                R = A;
                return
            end
            if nargin > 2
                B = horzcat(B, varargin{:});
            end
            if isa(A, 'sym') || isa(B, 'sym')
                R = [sym(A) sym(B)];
                return
            end
            [A, B, d] = common(A, B);
            R = simplify(rational([A.N B.N],d));
        end
        function R = reshape(R, varargin)
            R.N = reshape(R.N, varargin{:});
        end
        function R = vertcat(A, B, varargin)
            if nargin > 2
                B = vertcat(B, varargin{:});
            end
            if isa(A, 'sym') || isa(B, 'sym')
                R = [sym(A);sym(B)];
                return
            end
            [A, B, d] = common(A, B);
            R = simplify(rational([A.N;B.N],d));
        end
        function R = kron(A, B)
            A = rational(A);
            B = rational(B);
            R = simplify(rational(kron(A.N, B.N), A.d * B.d));
        end
        function R = plus(A, B)
            if isa(A, 'sym') || isa(B, 'sym')
                R = sym(A) + sym(B);
            else
                [A, B, d] = common(A, B);
                R = simplify(rational([A.N+B.N],d));
            end
        end
        function R = sort(R, varargin)
            R.N = sort(R.N, varargin{:});
        end
        function R = diff(A, varargin)
            R = diff(A.N, varargin{:});
        end
        function R = minus(A, B)
            if isa(A, 'sym') || isa(B, 'sym')
                R = sym(A) - sym(B);
            else
                [A, B, d] = common(A, B);
                R = simplify(rational([A.N-B.N],d));
            end
        end
        function R = unique(R, varargin)
            R.N = unique(R.N, varargin{:});
        end
        function bl = ne(A, B)
            [A, B] = common(A, B);
            bl = ne(A.N, B.N);
        end
        function bl = eq(A, B)
            [A, B] = common(A, B);
            bl = eq(A.N, B.N);
        end
        function bl = lt(A, B)
            [A, B] = common(A, B);
            bl = lt(A.N, B.N);
        end
        function bl = gt(A, B)
            [A, B] = common(A, B);
            bl = gt(A.N, B.N);
        end
        function bl = le(A, B)
            [A, B] = common(A, B);
            bl = le(A.N, B.N);
        end
        function bl = ge(A, B)
            [A, B] = common(A, B);
            bl = ge(A.N, B.N);
        end
        function R = mrdivide(A, b)
            if isa(A, 'sym') || isa(b, 'sym')
                R = sym(A)/sym(b);
                return
            end
            A = rational(A);
            b = rational(b);
            assert(isequal(size(b.N), [1 1]));
            assert(b.N ~= 0);
            [b.N b.d] = deal(b.d,b.N);
            if b.d < 0
                b.N = -b.N;
                b.d = -b.d;
            end
            R = A*b;
        end
        function [A, B, d] = common(A, B)
            A = rational(A);
            B = rational(B);
            d = lcm(A.d, B.d);
            A.N = A.N * (d/A.d);
            B.N = B.N * (d/B.d);
            A.d = d;
            B.d = d;
        end
        function R = rdivide(A, B)
            A = rational(A);
            B = rational(B);
            assert(all(all(B.N ~= 0)));
            m = rational.matlcm(B.N);
            B.N = B.d*m./B.N;
            B.d = m;
            B = simplify(B);
            R = A.*B;
        end
        function R1 = subsref(R, S)
            switch S(1).type
              case '.'
                R1 = builtin('subsref', R, S);
              case '()'
                R1 = simplify(rational(subsref(R.N, S), R.d));
              otherwise
                error('');
            end
        end
        function A = subsasgn(A, S, B)
            [A, B] = common(A, B);
            A.N = subsasgn(A.N, S, B.N);
            A = simplify(A);
        end
        function n = numel(R, varargin)
            n = numel(R.N, varargin{:});
        end
        function ind = end(R, k, n)
            szd = size(R.N);
            if k < n
                ind = szd(k);
            else
            ind = prod(szd(k:end));
            end
        end
        function disp(R)
            if R.d ~= 1
                disp(['1/' num2str(R.d) ' * ']);
            end
            disp(R.N);
        end
        function str = num2str(R)
            if R.d ~= 1;
                str = ['1/' num2str(R.d) ' * ' num2str(R.N)];
            else
                str = num2str(R.N);
            end
        end
    end
    methods(Static)
        function g = matgcd(M, filtered)
            if nargin == 1 || ~filtered
                M = abs(M(find(M ~= 0)));
            end
            M = M(:);
            n = size(M, 1);
            if mod(n,2) == 1
                g = M(n);
                n = n - 1;
            else
                g = 0;
            end
            switch n
              case 0
              case 2
                g = gcd(g, gcd(M(1), M(2)));
              otherwise
                g = gcd(g, ...
                        rational.matgcd(gcd(M(1:n/2), M(n/2+1:n)), true));
            end
        end
        function l = matlcm(M, filtered)
            if nargin == 1 || ~filtered
                M = abs(M(find(M ~= 0)));
            end
            M = M(:);
            n = size(M, 1);
            if mod(n,2) == 1
                l = M(n);
                n = n - 1;
            else
                l = 1;
            end
            switch n
              case 0
              case 2
                l = lcm(l, lcm(M(1), M(2)));
              otherwise
                l = lcm(l, ...
                        rational.matlcm(lcm(M(1:n/2), M(n/2+1:n)), true));
            end
        end
    end
end
