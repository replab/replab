classdef cyclotomic
% Matrix with elements in the cyclotomic field
%
% Requires Java
%
% Example:
%   >>> I = replab.cyclotomic.eye(3)
%       1  0  0
%       0  1  0
%       0  0  1

    properties %(Access = protected)
        mat % (cell(\*,\*) of cyclo.Cyclo): Matrix of cyclotomic numbers
    end

    methods (Static)

        function c = eye(n)
            c = replab.cyclotomic.zeros(n, n);
            mat = c.mat;
            one = javaMethod('one', 'cyclo.Cyclo');
            for i = 1:n
                mat{i,i} = one;
            end
            c = replab.cyclotomic(mat);
        end

        function c = zeros(m, n)
            mat = repmat({javaMethod('zero', 'cyclo.Cyclo')}, m, n);
            c = replab.cyclotomic(mat);
        end

        function c = ones(m, n)
            mat = repmat({javaMethod('one', 'cyclo.Cyclo')}, m, n);
            c = replab.cyclotomic(mat);
        end

        function c = E(n)
            c = replab.cyclotomic({javaMethod('e', 'cyclo.Cyclo')});
        end

        function c = fromStrings(strings)
            ja = javaMethod('parse', 'cyclo.Lab', strings(:));
            c = replab.cyclotomic.fromJavaArray(ja, size(strings));
        end

        function c = fromDoubles(doubles)
            c = replab.cyclotomic.fromJavaArray(javaMethod('fromDouble', 'cyclo.Lab', doubles(:)), size(doubles));
        end

        function c = sqrtRational(num, den)
            if nargin == 1
                c = replab.cyclotomic({com.faacets.gluon.Cyclotomic.sqrt(num)});
            else
                c = replab.cyclotomic({com.faacets.gluon.Cyclotomic.sqrt(num, den)});
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

    methods

        function res = matArray(self)
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
            assert(iscell(mat));
            self.mat = mat;
        end

        function disp(self)
            t = replab.compat.javaArrayToCell(javaMethod('print', 'cyclo.Lab', self.matArray));
            t = replab.str.Table(reshape(t, self.size), 'uniform', 0);
            disp(t);
        end

        function s = num2str(self)
            t = replab.compat.javaArrayToCell(javaMethod('print', 'cyclo.Lab', self.matArray));
            t = replab.str.Table(reshape(t, self.size), 'uniform', 0);
            s = t.format(1000, 1000);
        end

        function res = eq(self, rhs)
            res = reshape(javaMethod('eqv', 'cyclo.Lab', self.matArray, rhs.matArray), size(self.mat));
        end

        function res = conj(self)
            res = replab.cyclotomic.fromJavaArray(javaMethod('conjugate', 'cyclo.Lab', self.matArray), self.size);
        end

        function res = sqrt(self)
            res = replab.cyclotomic.fromJavaArray(javaMethod('sqrt', 'cyclo.Lab', self.matArray), self.size);
        end

        function res = plus(lhs, rhs)
            if isa(lhs, 'double')
                lhs = replab.cyclotomic.fromDoubles(lhs);
            end
            if isa(rhs, 'double')
                rhs = replab.cyclotomic.fromDoubles(rhs);
            end
            res = replab.cyclotomic.fromJavaArray(javaMethod('plus', 'cyclo.Lab', lhs.matArray, rhs.matArray), size(lhs));
        end

        function res = minus(lhs, rhs)
            if isa(lhs, 'double')
                lhs = replab.cyclotomic.fromDoubles(lhs);
            end
            if isa(rhs, 'double')
                rhs = replab.cyclotomic.fromDoubles(rhs);
            end
            res = replab.cyclotomic.fromJavaArray(javaMethod('minus', 'cyclo.Lab', lhs.matArray, rhs.matArray), size(lhs));
        end

        function res = times(lhs, rhs)
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

        function res = mtimes(lhs, rhs)
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
            n = size(self.mat, 1);
            assert(size(self.mat, 2) == n);
            res = replab.cyclotomic.fromJavaArray(javaMethod('power', 'cyclo.Lab', n, self.matArray, e), [n n]);
        end

        function res = inv(self)
            n = size(self.mat, 1);
            assert(size(self.mat, 2) == n);
            res = replab.cyclotomic.fromJavaArray(javaMethod('inverse', 'cyclo.Lab', n, self.matArray), [n n]);
        end

        function s = size(self, varargin)
            s = size(self.mat, varargin{:});
        end

        function l = length(self)
            l = length(self.mat);
        end

        function varargout = subsref(self, s)
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
            rhs = cellfun(@(a) a.mat, varargin, 'uniform', 0);
            res = horzcat(self.mat, rhs{:});
            res = replab.cyclotomic(res);
        end

        function res = vertcat(self, varargin)
            rhs = cellfun(@(a) a.mat, varargin, 'uniform', 0);
            res = vertcat(self.mat, rhs{:});
            res = replab.cyclotomic(res);
        end

    end

end
