classdef IrreducibleCommutant < replab.Equivariant
% Commutant algebra of an irreducible decomposition whose isotypic components have been harmonized
%
%
% Algebra of matrices that commute with an irreducible decomposition
%
% Note that the ``rep`` property must be of type `+replab.Irreducible`.

    methods

        function self = IrreducibleCommutant(irreducible)
            self = self@replab.Equivariant(irreducible, irreducible, 'commutant');
            assert(all(cellfun(@(c) c.isHarmonized, irreducible.components)));
        end

    end

    methods (Access = protected)

        function X1 = project_exact(self, X)
            n = self.repR.nComponents;
            blocks = cell(1, n);
            shift = 0;
            for i = 1:n
                iso = self.repR.component(i);
                d = iso.dimension;
                r = shift + (1:d);
                blocks{i} = iso.commutant.project(X(r, r));
                shift = shift + iso.dimension;
            end
            X1 = blkdiag(blocks{:});
        end

        function [X1, err] = project_double_sparse(self, X)
            n = self.repR.nComponents;
            blocks = cell(1, n);
            shift = 0;
            for i = 1:n
                iso = self.repR.component(i);
                d = iso.dimension;
                r = shift + (1:d);
                blocks{i} = iso.commutant.project(X(r, r));
                shift = shift + iso.dimension;
            end
            X1 = blkdiag(blocks{:});
            err = inf;
        end

    end

end
