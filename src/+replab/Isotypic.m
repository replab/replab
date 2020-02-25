classdef Isotypic < replab.SubRep
% Describes an isotypic component in the decomposition of a representation
%
% It is expressed as a subrepresentation of the representation being decomposed, however
% key methods are implemented more efficiently as more structure is available. In particular
% the computation of images is done in a way that minimizes numerical error and returns
% true block diagonal matrices.
%
% An isotypic component regroups equivalent irreducible representations, however not necessarily
% expressed in the same basis.
% Note that if the multiplicity is not one, there is a degeneracy in the basis of the copies, and
% the particular basis chosen is not deterministic.
%
% However the subspace spanned by an isotypic component as a whole is unique.
%
% We require that the embedding maps of the irreps do not overlap.

    properties (SetAccess = protected)
        irreps % (cell(1,\*) of `.SubRep`): Equivalent irreducible subrepresentations in this isotypic component
    end

    methods

        function self = Isotypic(parent, irreps)
        % Constructor
        %
        % Requirement: the embedding maps of the given irreps do not overlap
            m = length(irreps);
            assert(m >= 1, 'Isotypic component cannot be empty');
            assert(isa(parent, 'replab.Rep'));
            Bs = cell(1, m);
            Es = cell(m, 1);
            for i = 1:m
                ci = irreps{i};
                assert(isa(ci, 'replab.SubRep'));
                assert(isequal(ci.isIrreducible, true));
                Bs{1,i} = ci.B_internal;
                Es{i,1} = ci.E_internal;
            end
            B_internal = [Bs{:}];
            E_internal = vertcat(Es{:});
            self = self@replab.SubRep(parent, B_internal, E_internal);
            self.irreps = irreps;
            % mutable Rep properties
            areUnitary = cellfun(@(irr) irr.isUnitary, irreps, 'uniform', 0);
            if replab.util.areAllEqual(areUnitary)
                self.isUnitary = irreps{1}.isUnitary;
            end
            trivialDimensions = cellfun(@(irr) irr.trivialDimension, irreps, 'uniform', 0);
            if replab.util.areAllEqual(trivialDimensions)
                self.trivialDimension = m * irreps{1}.trivialDimension;
            end
            isDAC = cellfun(@(irr) irr.isDivisionAlgebraCanonical, irreps, 'uniform', 0);
            if replab.util.areAllEqual(isDAC)
                self.isDivisionAlgebraCanonical = irreps{1}.isDivisionAlgebraCanonical;
            end
        end

        function m = multiplicity(self)
        % Number of equivalent irreducible representations in this isotypic component
            m = length(self.irreps);
        end

        function d = irrepDimension(self)
        % Dimension of every single irreducible representation in this component
            d = self.irrep(1).dimension;
        end

        function n = nIrreps(self)
        % Returns the number of irreps in this isotypic component, which is their multiplicity
            n = length(self.irreps);
        end

        function c = irrep(self, i)
        % Returns the i-th copy of the irreducible representation
            c = self.irreps{i};
        end

        function P = projector(self)
        % Returns the projector on this isotypic component
            P = full(self.B_internal*self.E_internal);
        end

        function P = projectorOnIrrep(self, i)
        % Returns the projector on the i-th irreducible representation in this component
            Bi = self.irrep(i).B_internal;
            Ei = self.irrep(i).E_internal;
            P = full(Bi*Ei);
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
            if isa(self, 'replab.HarmonizedIsotypic')
                s = 'Isotypic component (harmonized)';
            else
                s = 'Isotypic component';
            end
            if self.overC
                rt = 'C';
            else
                fbsi = self.irrep(1).frobeniusSchurIndicator;
                if isempty(fbsi)
                    rt = '?';
                else
                    switch fbsi
                      case 0
                        rt = 'R';
                      case 1
                        rt = 'C';
                      case -1
                        rt = 'H';
                      otherwise
                        rt = '?';
                    end
                end
            end
            if self.multiplicity > 1
                s = sprintf('%s I(%d)x%s(%d)', s, self.multiplicity, rt, self.irrepDimension);
            else
                s = sprintf('%s %s(%d)', s, rt, self.irrepDimension);
            end
        end

        %% SubRep methods

        function iso = refine(self)
            iso = replab.rep.refineIsotypic(self, replab.Context.make);
        end

    end

end
