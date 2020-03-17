classdef H < replab.Str
% Barebones quaternion class, defining just enough for our purposes

    properties
        r % (double(\*,\*)): Real part
        i % (double(\*,\*)): Coeffient in the first quaternion unit
        j % (double(\*,\*)): Coeffient in the second quaternion unit
        k % (double(\*,\*)): Coeffient in the third quaternion unit
    end

    methods
        function self = H(r, i, j, k)
        % Constructor from the expansion as a d=4 real vector space
        %
        % The real part should always be provided. Otherwise, undefined parameters
        % are replaced by zeros.
        %
        % Special case: if only one argument is given:
        %   - if the argument is complex,
        %     it fills the real and 'i' parts of the quaternion.
        %   - if the argument is a quaternion, the quaternion is copied
            if nargin < 4 || isempty(k)
                k = zeros(size(r));
            end
            if nargin < 3 || isempty(j)
                j = zeros(size(r));
            end
            if nargin == 1 && isa(r, 'replab.quaternion.H')
                i = r.i;
                j = r.j;
                k = r.k;
                r = r.r;
            elseif nargin == 1 && ~isreal(r)
                i = imag(r);
                r = real(r);
            elseif nargin < 2 || isempty(i)
                i = zeros(size(r));
            end
            assert(isreal(r));
            assert(isreal(i));
            assert(isreal(j));
            assert(isreal(k));
            assert(isequal(size(r), size(i)));
            assert(isequal(size(r), size(j)));
            assert(isequal(size(r), size(k)));
            self.r = r;
            self.i = i;
            self.j = j;
            self.k = k;
        end
        function r1 = abs(q)
            r1 = sqrt(abs2(q));
        end
        function r1 = abs2(q)
            r1 = q.r.^2 + q.i.^2 + q.j.^2 + q.k.^2;
        end
        function qc = conj(q)
            qc = replab.quaternion.H(q.r, -q.i, -q.j, -q.k);
        end
        function qt = ctranspose(q)
            qt = q.conj.transpose;
        end
        function e = end(q, k, n)
            e = builtin('end', q.r, k, n);
        end
        function b = eq(q1, q2)
            q1 = replab.quaternion.H(q1);
            q2 = replab.quaternion.H(q2);
            b = (q1.r == q2.r) & (q1.i == q2.i) & (q1.j == q2.j) & (q1.k == q2.k);
        end
        function q = horzcat(varargin)
            arg1 = replab.quaternion.H(varargin{1});
            qr = arg1.r;
            qi = arg1.i;
            qj = arg1.j;
            qk = arg1.k;
            for i = 2:nargin
                arg = replab.quaternion.H(varargin{i});
                qr = horzcat(qr, arg.r);
                qi = horzcat(qi, arg.i);
                qj = horzcat(qj, arg.j);
                qk = horzcat(qk, arg.k);
            end
            q = replab.quaternion.H(qr, qi, qj, qk);
        end
        function q1 = inv(q)
            assert(isscalar(q.r));
            q1 = q.conj * (1 ./ q.abs2);
        end
        function q1 = inverse(q)
        % Scalar inverse
            q1 = q.conj * (1 ./ q.abs2);
        end
        function q3 = minus(q1, q2)
            q1 = replab.quaternion.H(q1);
            q2 = replab.quaternion.H(q2);
            q3 = replab.quaternion.H(q1.r - q2.r, q1.i - q2.i, q1.j - q2.j, q1.k - q2.k);
        end
        function q3 = mtimes(q1, q2)
            q1 = replab.quaternion.H(q1);
            q2 = replab.quaternion.H(q2);
            % 1 = -ii = -jj = -kk
            % i = jk = -kj = 1i = i1
            % j = ki = -ik = 1j = j1
            % k = ij = -ji = 1k = k1
            q3r = q1.r*q2.r - q1.i*q2.i - q1.j*q2.j - q1.k*q2.k;
            q3i = q1.j*q2.k - q1.k*q2.j + q1.r*q2.i + q1.i*q2.r;
            q3j = q1.k*q2.i - q1.i*q2.k + q1.r*q2.j + q1.j*q2.r;
            q3k = q1.i*q2.j - q1.j*q2.i + q1.r*q2.k + q1.k*q2.r;
            q3 = replab.quaternion.H(q3r, q3i, q3j, q3k);
        end
        function q3 = mrdivide(q1, q2)
            assert(isscalar(q2));
            q3 = rdivide(q1, q2);
        end
        function n = norm(q, p)
            if isvector(q.r)
                assert(nargin < 2 || p == 2);
                n = q.abs2;
                n = sqrt(sum(n(:)));
            elseif ismatrix(q.r)
                if nargin < 2
                    p = 2;
                end
                switch nargin
                  case 'fro'
                    n = q.abs2;
                    n = sqrt(sum(n(:)));
                  case 1
                    n = max(sum(q.abs, 1));
                  case inf
                    n = max(sum(q.abs, 2));
                  otherwise
                    s = replab.shortStr(p);
                    error(sprintf('Norm type %s undefined', s));
                end
            end
        end
        function q3 = plus(q1, q2)
            q1 = replab.quaternion.H(q1);
            q2 = replab.quaternion.H(q2);
            q3 = replab.quaternion.H(q1.r+q2.r, q1.i+q2.i, q1.j+q2.j, q1.k+q2.k);
        end
        function r = real(q)
            r = q.r;
        end
        function q3 = rdivide(q1, q2)
            q1 = replab.quaternion.H(q1);
            q2 = replab.quaternion.H(q2);
            q3 = q1 .* q2.inverse;
        end
        function q1 = subsasgn(q, s, p)
            q = replab.quaternion.H(q);
            p = replab.quaternion.H(p);
            qr = q.r;
            qi = q.i;
            qj = q.j;
            qk = q.k;
            qr(s.subs{:}) = p.r;
            qi(s.subs{:}) = p.i;
            qj(s.subs{:}) = p.j;
            qk(s.subs{:}) = p.k;
            q1 = replab.quaternion.H(qr, qi, qj, qk);
        end
        function varargout = subsref(q, s)
            assert(length(s) == 1);
            if isequal(s.type, '()')
                q1r = q.r(s.subs{:});
                q1i = q.i(s.subs{:});
                q1j = q.j(s.subs{:});
                q1k = q.k(s.subs{:});
                varargout{1} = replab.quaternion.H(q1r, q1i, q1j, q1k);
            else
                [varargout{1:nargout}] = builtin('subsref', q, s);
            end
        end
        function q3 = times(q1, q2)
            q1 = replab.quaternion.H(q1);
            q2 = replab.quaternion.H(q2);
            % 1 = -ii = -jj = -kk
            % i = jk = -kj = 1i = i1
            % j = ki = -ik = 1j = j1
            % k = ij = -ji = 1k = k1
            q3r = q1.r.*q2.r - q1.i.*q2.i - q1.j.*q2.j - q1.k.*q2.k;
            q3i = q1.j.*q2.k - q1.k.*q2.j + q1.r.*q2.i + q1.i.*q2.r;
            q3j = q1.k.*q2.i - q1.i.*q2.k + q1.r.*q2.j + q1.j.*q2.r;
            q3k = q1.i.*q2.j - q1.j.*q2.i + q1.r.*q2.k + q1.k.*q2.r;
            q3 = replab.quaternion.H(q3r, q3i, q3j, q3k);
        end
        function q1 = transpose(q)
            q1 = replab.quaternion.H(q.r.', q.i.', q.j.', q.k.');
        end
        function qm = uminus(q)
            qm = replab.quaternion.H(-q.r, -q.i, -q.j, -q.k);
        end
        function qp = uplus(q)
            qp = q;
        end
        function q = vertcat(varargin)
            arg1 = replab.quaternion.H(varargin{1});
            qr = arg1.r;
            qi = arg1.i;
            qj = arg1.j;
            qk = arg1.k;
            for i = 2:nargin
                arg = replab.quaternion.H(varargin{i});
                qr = vertcat(qr, arg.r);
                qi = vertcat(qi, arg.i);
                qj = vertcat(qj, arg.j);
                qk = vertcat(qk, arg.k);
            end
            q = replab.quaternion.H(qr, qi, qj, qk);
        end
    end
    methods (Static)
        function q = zeros(varargin)
            r = zeros(varargin{:});
            q = replab.quaternion.H(r);
        end
        function q = eye(varargin)
            r = eye(varargin{:});
            q = replab.quaternion.H(r);
        end
    end
end
