classdef cyclotomic

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
            mat = cellfun(@(s) com.faacets.gluon.CycloUtils.parse(s), strings, 'uniform', 0);
            c = replab.cyclotomic(mat);
        end

        function c = fromDouble(doubles)
            mat = arrayfun(@(d) com.faacets.gluon.CycloUtils.fromDouble(d), doubles, 'uniform', 0);
            c = replab.cyclotomic(mat);
        end

    end

    methods

        function self = cyclotomic(mat)
            assert(iscell(mat));
            self.mat = mat;
        end

        function disp(self)
            t = replab.str.Table(cellfun(@(c) char(com.faacets.gluon.CycloUtils.printCyclo(c)), self.mat, 'uniform', 0));
            disp(t);
        end

        function res = conj(self)
            res = replab.cyclotomic(cellfun(@(c) c.conjugate, self.mat, 'uniform', 0));
        end

        function res = plus(self, rhs)
            I = com.faacets.gluon.Interface.compile({
                'import com.faacets.gluon._; import scalin.immutable.{Mat, DenseMat}; import cyclo.Cyclo; import scalin.immutable.dense._'
                'implicit val forCyclo: ArgAdapter[Array[Cyclo]] = ArgAdapter[Array[Cyclo]]("Array[Cyclo]")(_.asInstanceOf[Array[_]].map(_.asInstanceOf[Cyclo]))'
                'Interface[Int, Int, Array[Cyclo], Array[Cyclo], Array[Cyclo]]("plus", "nR", "nC", "lhs", "rhs") {'
                '(nR, nC, lhs, rhs) =>'
                'val l = DenseMat.tabulate(nR, nC)( (r, c) => lhs(r + nR*c) )'
                'val r = DenseMat.tabulate(nR, nC)( (r, c) => rhs(r + nR*c) )'
                'val z = l + r'
                'Array.tabulate(nR, nC)( (r, c) => z(r, c) ).flatten'
                '}'
                   });
            nRows = size(self.mat, 1);
            nCols = size(self.mat, 2);
            assert(size(rhs.mat, 1) == nRows);
            assert(size(rhs.mat, 2) == nCols);
            res1 = I.call4(nRows, nCols, self.mat(:), rhs.mat(:));
            res = replab.cyclotomic(reshape(cell(res1), nRows, nCols));
        end

        function res = minus(self, rhs)
            I = com.faacets.gluon.Interface.compile({
                'import com.faacets.gluon._; import scalin.immutable.{Mat, DenseMat}; import cyclo.Cyclo; import scalin.immutable.dense._'
                'implicit val forCyclo: ArgAdapter[Array[Cyclo]] = ArgAdapter[Array[Cyclo]]("Array[Cyclo]")(_.asInstanceOf[Array[_]].map(_.asInstanceOf[Cyclo]))'
                'Interface[Int, Int, Array[Cyclo], Array[Cyclo], Array[Cyclo]]("plus", "nR", "nC", "lhs", "rhs") {'
                '(nR, nC, lhs, rhs) =>'
                'val l = DenseMat.tabulate(nR, nC)( (r, c) => lhs(r + nR*c) )'
                'val r = DenseMat.tabulate(nR, nC)( (r, c) => rhs(r + nR*c) )'
                'val z = l - r'
                'Array.tabulate(nR, nC)( (r, c) => z(r, c) ).flatten'
                '}'
                   });
            nRows = size(self.mat, 1);
            nCols = size(self.mat, 2);
            assert(size(rhs.mat, 1) == nRows);
            assert(size(rhs.mat, 2) == nCols);
            res1 = I.call4(nRows, nCols, self.mat(:), rhs.mat(:));
            res = replab.cyclotomic(reshape(cell(res1), nRows, nCols));
        end

        function res = mtimes(self, rhs)
            I = com.faacets.gluon.Interface.compile({
                'import com.faacets.gluon._; import scalin.immutable.{Mat, DenseMat}; import cyclo.Cyclo; import scalin.immutable.dense._'
                'implicit val forCyclo: ArgAdapter[Array[Cyclo]] = ArgAdapter[Array[Cyclo]]("Array[Cyclo]")(_.asInstanceOf[Array[_]].map(_.asInstanceOf[Cyclo]))'
                'Interface[Int, Int, Int, Array[Cyclo], Array[Cyclo], Array[Cyclo]]("plus", "nR1", "nC1", "nC2", "lhs", "rhs") {'
                '(nR1, nC1, nC2, lhs, rhs) =>'
                'val l = DenseMat.tabulate(nR1, nC1)( (r, c) => lhs(r + nR1*c) )'
                'val r = DenseMat.tabulate(nC1, nC2)( (r, c) => rhs(r + nC1*c) )'
                'val z = l * r'
                'Array.tabulate(nR1, nC2)( (r, c) => z(r, c) ).flatten'
                '}'
                   });
            nR1 = size(self.mat, 1);
            nC1 = size(self.mat, 2);
            nR2 = size(rhs.mat, 1);
            assert(nC1 == nR2);
            nC2 = size(rhs.mat, 2);
            res1 = I.call5(nR1, nC1, nC2, self.mat(:), rhs.mat(:));
            res = replab.cyclotomic(reshape(cell(res1), nR1, nC2));
        end

        function s = size(self)
            s = size(self.mat);
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
            els = self.mat(:);
            res = zeros(length(els), 1);
            for i = 1:length(els)
                cd = com.faacets.gluon.CycloUtils.toDouble(els{i});
                res(i) = cd(1) + 1i*cd(2);
            end
            res = reshape(res, size(self.mat));
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
