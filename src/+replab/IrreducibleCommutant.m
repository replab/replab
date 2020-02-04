classdef IrreducibleCommutant < replab.Equivariant
% Algebra of matrices that commute with an irreducible decomposition
%
% Note that the ``rep`` property must be of type `+replab.Irreducible`.

    methods

        function self = IrreducibleCommutant(irreducible)
            self = self@replab.Equivariant(irreducible, irreducible, 'commutant');
        end

        function X = projectAndReduceFromParent(self, X)
            n = self.rep.nComponents;
            blocks = cell(1, n);
            for i = 1:n
                iso = self.rep.component(i);
                blocks{i} = iso.commutant.projectAndReduceFromParent(X);
            end
            X = blkdiag(blocks{:});
        end

        function X = projectAndReduce(self, X)
            n = self.rep.nComponents;
            blocks = cell(1, n);
            shift = 0;
            for i = 1:n
                iso = self.rep.component(i);
                r = shift + (1:iso.dimension);
                blocks{i} = iso.commutant.projectAndReduce(X(r, r));
                shift = shift + iso.dimension;
            end
            X = blkdiag(blocks{:});
        end

        function [X err] = project(self, X)
            n = self.rep.nComponents;
            blocks = cell(1, n);
            shift = 0;
            for i = 1:n
                iso = self.rep.component(i);
                d = iso.dimension;
                r = shift + (1:d);
                blocks{i} = iso.commutant.project(X(r, r));
                shift = shift + iso.dimension;
            end
            X = blkdiag(blocks{:});
            err = NaN;
        end

    end

end
