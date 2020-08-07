classdef cyclotomic
% Matrix of cyclotomic numbers
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
            for i = 1:n
                mat{i,i} = cyclo.Cyclo.one;
            end
            c = replab.cyclotomic(mat);
        end

        function c = zeros(m, n)
            mat = repmat({cyclo.Cyclo.zero}, m, n);
            c = replab.cyclotomic(mat);
        end

        function c = ones(m, n)
            mat = repmat({cyclo.Cyclo.one}, m, n);
            c = replab.cyclotomic(mat);
        end

        function c = E(n)
            c = replab.cyclotomic({cyclo.Cyclo.e(n)});
        end

        function c = fromString(strings)
            mat = cellfun(@(s) com.faacets.gluon.Cyclotomic.parse(s), strings, 'uniform', 0);
            c = replab.cyclotomic(mat);
        end

        function c = fromDouble(doubles)
            mat = arrayfun(@(d) com.faacets.gluon.Cyclotomic.fromDouble(d), doubles, 'uniform', 0);
            c = replab.cyclotomic(mat);
        end

    end

    methods

        function res = matArray(self)
            m = self.mat;
            res = [m{:}];
        end

        function self = cyclotomic(mat)
            assert(iscell(mat));
            self.mat = mat;
        end

        function disp(self)
            t = replab.str.Table(cellfun(@(c) char(com.faacets.gluon.Cyclotomic.printCyclo(c)), self.mat, 'uniform', 0));
            disp(t);
        end

        function res = conj(self)
            mat = reshape(cell(com.faacets.gluon.Cyclotomic.conjugate([self.matArray])), size(self.mat));
            res = replab.cyclotomic(mat);
        end

        function res = plus(self, rhs)
            mat = reshape(cell(com.faacets.gluon.Cyclotomic.plus([self.matArray], [rhs.matArray])), size(self.mat));
            res = replab.cyclotomic(mat);
        end

        function res = minus(self, rhs)
            mat = reshape(cell(com.faacets.gluon.Cyclotomic.minus([self.matArray], [rhs.matArray])), size(self.mat));
            res = replab.cyclotomic(mat);
        end

        function res = times(self, rhs)
            mat = reshape(cell(com.faacets.gluon.Cyclotomic.pw_times([self.matArray], [rhs.matArray])), size(self.mat));
            res = replab.cyclotomic(mat);
        end

        function res = mtimes(self, rhs)
            l = size(self.mat, 1);
            m = size(self.mat, 2);
            assert(m == size(self.mat, 1));
            n = size(rhs.mat, 2);
            mat = com.faacets.gluon.Cyclotomic.times(l, m, n, self.matArray, rhs.matArray);
            res = replab.cyclotomic(reshape(cell(mat), [l n]));
        end

        function res = inv(self)
            n = size(self.mat, 1);
            assert(size(self.mat, 2) == n);
            mat = com.faacets.gluon.Cyclotomic.inverse(n, self.matArray);
            res = replab.cyclotomic(reshape(cell(mat), [n n]));
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
                        arg = replab.cyclotomic.fromDouble(val);
                        self = subsasgn(self, s(1), arg);
                    end
                end
              otherwise
                error('Not a valid indexing expression')
            end
        end

        function res = double(self)
            pairs = com.faacets.gluon.Cyclotomic.toDouble([self.matArray]);
            r = pairs(1,:);
            i = pairs(2,:);
            res = reshape(r, size(self.mat)) + 1i * reshape(i, size(self.mat));
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
