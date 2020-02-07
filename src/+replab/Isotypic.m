classdef Isotypic < replab.SubRep
% Describes an isotypic component in the decomposition of a representation
%
% It is expressed as a subrepresentation of the representation being decomposed, however
% key methods are implemented more efficiently as more structure is available. In particular
% the computation of images is done in a way that minimizes numerical error and returns
% true block diagonal matrices.
%
% An isotypic component regroups equivalent irreducible representations, expressed in the same basis.
% Note that if the multiplicity is not one, there is a degeneracy in the basis of the copies, and
% the particular basis choosen is not deterministic.
%
% However the subspace spanned by an isotypic component as a whole is unique.

    properties
        irreps % (cell(1,*) of `.SubRep`): Equivalent irreducible subrepresentations in this isotypic component
        multiplicity % (integer): Number of equivalent irreducible representations in this isotypic component
        irrepDimension % (integer): Dimension of each irreducible representation in this component
    end

    methods

        function self = Isotypic(parent, irreps)
            assert(length(irreps) >= 1, 'Isotypic component cannot be empty');
            assert(isa(parent, 'replab.Rep'));
            for i = 1:length(irreps)
                ci = irreps{i};
                assert(isa(ci, 'replab.SubRep'));
                assert(ci.isKnownIrreducible);
            end
            Us = cellfun(@(sr) sr.U, irreps, 'uniform', 0);
            U = vertcat(Us{:});
            nbs = cellfun(@(sr) sr.niceBasis, irreps, 'uniform', 0);
            niceBasis = replab.NiceBasis.vertcat(nbs);
            self = self@replab.SubRep(parent, U, niceBasis);
            self.irreps = irreps;
            self.multiplicity = length(irreps);
            self.irrepDimension = irreps{1}.dimension;
        end

        function n = nIrreps(self)
        % Returns the number of irreps = the multiplicity
            n = self.multiplicity;
        end

        function c = irrep(self, i)
        % Returns the i-th copy of the irreducible representation
            c = self.irreps{i};
        end

        function P = projector(self)
        % Returns the projector on this isotypic component
            P = full(self.U'*self.U);
        end

        function P = projectorOnIrrep(self, i)
        % Returns the projector on the i-th irreducible representation in this component
            Ui = self.irrep(i).U;
            P = full(Ui'*Ui);
        end

        %% Str methods

        function names = hiddenFields(self)
            names = hiddenFields@replab.SubRep(self);
            names{1, end+1} = 'irreps';
        end

        function [names values] = additionalFields(self)
            [names values] = additionalFields@replab.SubRep(self);
            for i = 1:self.nIrreps
                names{1, end+1} = sprintf('irrep(%d)', i);
                values{1, end+1} = self.irrep(i);
            end
        end

        function s = headerStr(self)
            if isequal(self.irrep(1).field, 'C')
                rt = 'C';
            else
                rt = self.irrep(1).irrepInfo.divisionAlgebra;
                assert(~isempty(rt));
            end
            if self.multiplicity > 1
                s = sprintf('Isotypic component I(%d)x%s(%d)', self.multiplicity, rt, self.irrepDimension);
            else
                s = sprintf('Isotypic component %s(%d)', rt, self.irrepDimension);
            end
        end

        %% Rep methods

        function rho = image(self, g)
            p = self.parent.image(g);
            U = self.irrep(1).U;
            rho = U*p*U';
            for i = 2:self.nIrreps
                U = self.irrep(i).U;
                rho = rho + U*p*U';
            end
            rho = rho / self.nIrreps;
            rho = kron(eye(self.nIrreps), rho);
        end

        function c = commutant(self)
            if isempty(self.commutant_)
                if self.overC
                    self.commutant_ = replab.IsotypicSimpleCommutant(self);
                else
                    switch self.irrep(1).irrepInfo.divisionAlgebra
                      case 'R'
                        self.commutant_ = replab.IsotypicSimpleCommutant(self);
                      case 'C'
                        self.commutant_ = replab.IsotypicComplexCommutant(self);
                      case 'H'
                        self.commutant_ = replab.IsotypicQuaternionCommutant(self);
                    end
                end
            end
            c = self.commutant_;
        end

    end

end
