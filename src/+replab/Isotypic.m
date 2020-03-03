classdef Isotypic < replab.SubRep
% Describes an isotypic component in the decomposition of a representation
%
% It is expressed as a subrepresentation of the representation being decomposed, however
% key methods are implemented more efficiently as more structure is available. In particular
% the computation of images is done in a way that minimizes numerical error and returns
% true block diagonal matrices.
%
% An isotypic component regroups equivalent irreducible representations, however not necessarily
% expressed in the same basis (that would be `.HarmonizedIsotypic`).
% Note that if the multiplicity is not one, there is a degeneracy in the basis of the copies, and
% the particular basis chosen is not deterministic.
%
% However the subspace spanned by an isotypic component as a whole is unique.

    properties (SetAccess = protected)
        irreps % (cell(1,\*) of `.SubRep`): Equivalent irreducible subrepresentations in this isotypic component
    end

    methods (Static)

        function iso = fromIrreps(parent, irreps)
        % Builds an isotypic canonical component from equivalent subrepresentations
        %
        % Args:
        %   parent (`+replab.Rep`): Representation being decomposed
        %   irreps (cell(1,\*) of `+replab.SubRep`): Equivalent irreducible subrepresentations of ``rep``
        %
        % Returns:
        %   `+replab.Isotypic`: The corresponding isotypic component
            assert(length(irreps) >= 1, 'Isotypic component cannot be empty');
            m = length(irreps);
            for i = 1:m
                s = irreps{i};
                assert(isa(s, 'replab.SubRep'));
                assert(s.parent == parent);
                assert(isequal(s.isIrreducible, true));
            end
            d = irreps{1}.dimension;
            if isequal(parent.isUnitary, true)
                b = cellfun(@(s) isequal(s.B_internal, s.E_internal'), irreps);
                if all(b)
                    % all bases are unitary, parent is unitary, we use orthogonality
                    iso = replab.Isotypic(rep, sub);
                    Es = cell(m, 1);
                    for i = 1:m
                        Es{i,1} = irreps{i}.E_internal;
                    end
                    E_internal = vertcat(Es{:});
                    iso = replab.Isotypic(parent, irreps, E_internal);
                    return
                end
            end
            Bs = cell(1, m);
            for i = 1:m
                Bs{1,i} = irreps{i}.B_internal;
            end
            Biso = [Bs{:}];
            subiso = parent.subRep(Biso);
            E_internal = subiso.E_internal;
            iso = replab.Isotypic(parent, irreps, E_internal);
        end

    end

    methods

        function self = Isotypic(parent, irreps, E_internal)
        % Constructs an isotypic component of a parent representation
        %
        % The basis of the isotypic component comes from the concatenation of the bases of the irreps,
        % while its embedding map is supplied as an argument.
        %
        % The static method `.fromIrreps` computes this embedding map if necessary.
        %
        % Args:
        %   parent (`+replab.Rep`): Parent representation of which we construct a subrepresentation
        %   irreps (cell(1,\*) of `+replab.SubRep`): Irreducible subrepresentations
        %   E_internal (double(\*,\*), may be sparse): Embedding map of the isotypic component
            m = length(irreps);
            assert(m >= 1, 'Isotypic component cannot be empty');
            assert(isa(parent, 'replab.Rep'));
            Bs = cell(1, m);
            for i = 1:m
                ci = irreps{i};
                assert(isa(ci, 'replab.SubRep'));
                assert(ci.parent == parent);
                assert(isequal(ci.isIrreducible, true));
                Bs{1,i} = ci.B_internal;
            end
            B_internal = [Bs{:}];
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
        %
        % Returns:
        %   integer: Multiplicity
            m = length(self.irreps);
        end

        function d = irrepDimension(self)
        % Dimension of every single irreducible representation in this component
        %
        % Returns:
        %   integer: Irreducible representation dimension
            d = self.irrep(1).dimension;
        end

        function n = nIrreps(self)
        % Returns the number of irreps in this isotypic component, which is their multiplicity
        %
        % Returns:
        %   integer: Number of irreducible representations
            n = length(self.irreps);
        end

        function c = irrep(self, i)
        % Returns the i-th copy of the irreducible representation
        %
        % Returns:
        %   `.SubRep`: Irreducible subrepresentation of `.parent`
            c = self.irreps{i};
        end

        function P = projector(self)
        % Returns the projector on this isotypic component
        %
        % Returns:
        %   double(\*,\*): Projector matrix on the isotypic component
            P = full(self.B_internal*self.E_internal);
        end

        function P = projectorOnIrrep(self, i)
        % Returns the projector on the i-th irreducible representation in this component
        %
        % Returns:
        %   double(\*,\*): Projector matrix on the irreducible representation
            range = (i-1)*self.irrepDimension+(1:self.irrepDimension);
            Bi = self.irrep(i).B_internal;
            Ei = self.E_internal(range,:);
            P = full(Bi*Ei);
        end

        function iso = harmonize(self, context)
        % Harmonizes the isotypic component
        %
        % Returns:
        %   `.HarmonizedIsotypic`: Isotypic component with bases harmonized
            if isa(self, 'replab.HarmonizedIsotypic')
                iso = self;
            else
                if isempty(context)
                    c = replab.Context.make;
                else
                    c = context;
                end
                iso = replab.irreducible.harmonizeIsotypic(self, c);
                if isempty(context)
                    c.close;
                end
            end
        end

        function A = changeOfBasis(self, i, j, context)
        % Returns change of basis matrices that relate two irreducible representations
        %
        % Args:
        %   i (integer): Index of an irreducible representation
        %   j (integer): Index of an irreducible representation
        %   context (`+replab.Context`, optional): Sampling context
        % Returns:
        %   double(\*,\*): ``A`` such that ``A * self.irrep(j).image(g) * inv(A) = self.irrep(i).image(g)``
            if nargin < 4
                context = replab.Context.make;
            end
            C = self.parent.commutant.sampleInContext(context, 1);
            A = full(self.irrep(i).E_internal * C * self.irrep(j).B_internal);
            % TODO: is this choice of normalized sensical for nonunitary representations?
            A = A * sqrt(self.irrepDimension/real(trace(A*A'))) * sign(A(1,1));
            if nargin < 4
                context.close;
            end
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
                      case 1
                        rt = 'R';
                      case 0
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
            if isequal(self.trivialDimension, self.dimension)
                s = [s ' (trivial)'];
            end
        end

        %% SubRep methods

        function iso = refine(self)
            iso = replab.rep.refineIsotypic(self, replab.Context.make);
        end

    end

end
