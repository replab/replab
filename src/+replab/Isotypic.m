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
        irrepDimension % (integer): Dimension of every single irreducible representation in this component
        isHarmonized % (logical): Whether all irreducible subrepresentations are expressed in the same basis
    end

    methods (Static)

        function iso = fromTrivialSubRep(trivial)
        % Builds an isotypic component from the (full) trivial subrepresentation
        %
        % Args:
        %   trivial (`+replab.SubRep`): Maximal trivial subrepresentation of ``parent``
        %
        % Returns:
        %   `+replab.Isotypic`: The corresponding trivial isotypic component
            parent = trivial.parent;
            dT = trivial.dimension;
            irreps = cell(1, dT);
            I = trivial.injection_internal;
            P = trivial.projection_internal;
            for i = 1:dT
                irreps{i} = parent.subRep(I(:,i), 'projection', P(i,:), 'isIrreducible', true, 'frobeniusSchurIndicator', 1, 'trivialDimension', 1, 'isUnitary', true);
            end
            isHarmonized = true;
            irrepDimension = 1;
            iso = replab.Isotypic(parent, irreps, P, irrepDimension, isHarmonized);
        end

        function iso = fromIrreps(parent, irreps, irrepDimension, isHarmonized)
        % Builds an isotypic component from equivalent subrepresentations
        %
        % All irreps of the isotypic component must be provided, and their injection maps must be linearly
        % independent.
        %
        % Args:
        %   parent (`+replab.Rep`): Representation being decomposed
        %   irreps (cell(1,\*) of `+replab.SubRep`): Equivalent irreducible subrepresentations of ``parent``
        %   irrepDimension (integer): Dimension of a single irreducible representation
        %   isHarmonized (logical): Whether all irreps are expressed in the same basis
        %
        % Returns:
        %   `+replab.Isotypic`: The corresponding isotypic component
            if isempty(irreps)
                iso = replab.Isotypic(parent, irreps, zeros(0, parent.dimension), irrepDimension, true);
                return
            end
            m = length(irreps);
            if m == 1
                % Single multiplicity? Reuse projection
                iso = replab.Isotypic(parent, irreps, irreps{1}.projection_internal, irrepDimension, isHarmonized);
                return
            end
            mapsAreUnitary = cellfun(@(s) all(all(s.injection_internal == s.projection_internal')), irreps);
            if parent.knownUnitary && all(mapsAreUnitary)
                % all maps are unitary, parent is unitary, we use orthogonality
                projections = cell(1, m);
                for i = 1:m
                    projections{i} = irreps{i}.projection_internal;
                end
                projection = vertcat(projections{:});
                iso = replab.Isotypic(parent, irreps, projection, irrepDimension, isHarmonized);
                return
            end
            % non unitary case
            injections = cell(1, m);
            for i = 1:m
                injections{i} = irreps{i}.injection_internal;
            end
            injection = horzcat(injections{:});
            sub = parent.subRep(injection);
            iso = replab.Isotypic(parent, irreps, sub.projection, irrepDimension, isHarmonized);
        end

    end

    methods

        function self = Isotypic(parent, irreps, projection, irrepDimension, isHarmonized)
        % Constructs an isotypic component of a parent representation
        %
        % The injection map of the isotypic component comes from the sum of the injection maps of the irreps,
        % while its projection map is supplied as an argument.
        %
        % The static method `.fromIrreps` computes this embedding map if necessary.
        %
        % Additional keyword arguments are passed to the `.Rep` constructor.
        %
        % Args:
        %   parent (`+replab.Rep`): Parent representation of which we construct a subrepresentation
        %   irreps (cell(1,\*) of `+replab.SubRep`): Irreducible subrepresentations
        %   projection (double(\*,\*), may be sparse): Embedding map of the isotypic component
        %   irrepDimension (integer): Dimension of a single irreducible representation
        %   isHarmonized (logical): Whether all irreps are expressed in the same basis
            m = length(irreps);
            assert(isa(parent, 'replab.Rep'));
            assert(all(cellfun(@(irrep) irrep.parent == parent, irreps)));
            assert(all(cellfun(@(irrep) isa(irrep, 'replab.SubRep'), irreps)));
            assert(all(cellfun(@(irrep) irrep.isIrreducible, irreps)));
            args = {};
            if all(cellfun(@(irrep) irrep.knownUnitary, irreps))
                args = horzcat(args, {'isUnitary', true});
            end
            if all(cellfun(@(irrep) irrep.inCache('trivialDimension'), irreps))
                td = sum(cellfun(@(irrep) irrep.trivialDimension, irreps));
                args = horzcat(args, {'trivialDimension', td});
            end
            if all(cellfun(@(irrep) irrep.inCache('isDivisionAlgebraCanonical'), irreps))
                dac = all(cellfun(@(irrep) irrep.isDivisionAlgebraCanonical, irreps));
            end
            args = horzcat(args, {'isIrreducible', m > 1});
            injections = cellfun(@(irrep) irrep.injection_internal, irreps, 'uniform', 0);
            if m == 0
                injection = zeros(parent.dimension, 0);
            else
                injection = [injections{:}];
            end
            self = self@replab.SubRep(parent, injection, projection, args{:});
            self.irreps = irreps;
            self.irrepDimension = irrepDimension;
            self.isHarmonized = isHarmonized;
        end

        function m = multiplicity(self)
        % Number of equivalent irreducible representations in this isotypic component
        %
        % Returns:
        %   integer: Multiplicity
            m = length(self.irreps);
        end

        function r = irrepRange(self, i)
        % Returns the coefficient range where the i-th irrep lies
        %
        % Args:
        %   i (integer): Irreducible representation index
        %
        % Returns:
        %   integer(1,\*): Range of the rows/columns block of the i-th irrep
            id = self.irrepDimension;
            r = (1:id)+id*(i-1);
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

        function iso = harmonize(self, context)
        % Harmonizes the isotypic component
        %
        % Returns:
        %   `.Isotypic`: Isotypic component with `.isHarmonized` true
            if self.isHarmonized
                iso = self;
            elseif self.trivialDimension == self.dimension
                isHarmonized = true;
                iso = replab.Isotypic(self.parent, self.irreps, self.projection, self.irrepDimension, isHarmonized);
            else
                if nargin < 2 || isempty(context)
                    close = true;
                    c = replab.Context.make;
                else
                    close = false;
                    c = context;
                end
                iso = replab.irreducible.Isotypic_harmonize(self, c);
                if close
                    c.close;
                end
            end
        end

    end

    methods (Access = protected) % Implementations

        function rho = image_double_sparse(self, g)
            if self.isHarmonized && self.multiplicity > 0
                p = self.parent.image(g, 'double/sparse');
                I = self.irrep(1).injection;
                P = self.irrep(1).projection;
                rho = P*p*I;
                for i = 2:self.nIrreps
                    I = self.irrep(i).injection;
                    P = self.irrep(i).projection;
                    rho = rho + P*p*I;
                end
                rho = rho / self.nIrreps;
                rho = kron(eye(self.nIrreps), rho);
            else
                rho = image_double_sparse@replab.SubRep(self, g);
            end
        end

        function rho = image_exact(self, g)
            if self.isHarmonized && self.multiplicity > 0
                rho = self.irrep(1).image(g, 'exact');
                rho = kron(replab.cyclotomic.eye(self.nIrreps), rho);
            else
                rho = image_exact@replab.SubRep(self, g);
            end
        end

    end

    methods % Implementations

        % Str

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
                      case -2
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

        % Obj

        function l = laws(self)
            l = replab.laws.IsotypicLaws(self);
        end

        % SubRep

        function iso1 = refine(self, varargin)
        % Refines this isotypic component
        %
        % Assumes that the isotypic component is already close to an exact isotypic component, and refines its subspace
        % by an iterative procedure applied on its `.injection` and `.projection` maps.
        %
        % This procedure preserves biorthogonality of the irreducible representations, but not necessarily harmonization.
        %
        % Keyword Args:
        %   largeScale (logical or ``[]``, optional): Whether to use the large-scale version of the algorithm, default ``[]`` (automatic selection)
        %   numNonImproving (integer, optional): Number of non-improving steps before stopping the large-scale algorithm, default ``20``
        %   nSamples (integer, optional): Number of samples to use in the large-scale version of the algorithm, default ``5``
        %   nInnerIterations (integer, optional): Number of inner iterations in the medium-scale version of the algorithm, default ``10``
        %   maxIterations (integer, optional): Maximum number of (outer) iterations, default ``1000``
        %
        % Returns:
        %   `.Isotypic`: Isotypic component with refined subspace (injection/projection maps)
            replab.log(1, 'Refining isotypic component');
            args = struct('numNonImproving', 20, 'largeScale', self.parent.dimension > 1000, 'nSamples', 5, 'nInnerIterations', 10, 'maxIterations', 1000);
            args = replab.util.populateStruct(args, varargin);
            D = self.parent.dimension;
            m = self.multiplicity;
            irreps = cell(1, m);
            if self.parent.knownUnitary
                Q0 = self.injection;
                if ~all(all(Q0 == self.projection'))
                    [Q0, ~] = qr(Q0, 0);
                end
                Qo = zeros(D, 0);
                for i = 1:m
                    Q = Q0(:, self.irrepRange(i));
                    if args.largeScale
                        Q1 = replab.rep.refine_unitaryLargeScale(self.parent, Q, args.numNonImproving, args.nSamples, args.maxIterations, Qo);
                    else
                        Q1 = replab.rep.refine_unitaryMediumScale(self.parent, Q, args.nInnerIterations, args.maxIterations, Qo);
                    end
                    irreps{i} = self.parent.subRep(Q1, 'projection', Q1', 'isUnitary', true, 'isIrreducible', true);
                    Qo = [Qo Q1];
                end
                isHarmonized = false;
                iso1 = replab.Isotypic(self.parent, irreps, Qo', self.irrepDimension, isHarmonized);
            else
                I0 = self.injection;
                P0 = self.projection;
                Io = zeros(D, 0);
                Po = zeros(0, D);
                for i = 1:m
                    I = I(:, self.irrepRange(i));
                    P = P(self.irrepRange(i), :);
                    if args.largeScale
                        [I1, P1] = replab.rep.refine_nonUnitaryLargeScale(self.parent, I, P, args.numNonImproving, args.nSamples, args.maxIterations, Io, Po);
                    else
                        [I1, P1] = replab.rep.refine_nonUnitaryMediumScale(self.parent, I, P, args.nInnerIterations, args.maxIterations, Io, Po);
                    end
                    irreps{i} = self.parent.subRep(I1, 'projection', P1, 'isIrreducible', true);
                    Io = [Io I1];
                    Po = [Po; P1];
                end
                isHarmonized = false;
                iso1 = replab.Isotypic(self.parent, irreps, Po, self.irrepDimension, isHarmonized);
            end
        end

    end

end
