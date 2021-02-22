classdef Isotypic < replab.SubRep
% Describes an isotypic component in the decomposition of a representation
%
% It is expressed as a subrepresentation of the representation being decomposed, however key methods are implemented
% more efficiently as more structure is available. In particular the computation of images is done in a way that
% minimizes numerical error and returns true block diagonal matrices.
%
% An isotypic component regroups equivalent irreducible representations expressed in the same basis
% Note that if the multiplicity is not one, there is a degeneracy in the picking of a basis of of the multiplicity
% space, and the basis is not necessarily chosen in a deterministic way.
%
% However the subspace spanned by an isotypic component as a whole is unique.
%
% To cater for empty isotypic components, the isotypic component stores a "model" of the irreducible representation.
%
% The isotypic stores its own projection map, as the following equation is not necessarily satisfied for individual
% irreps:
%
% `` irreps{i}.projection * irrep{j}.injection = (i == j) * eye(d) ``

    properties (SetAccess = protected)
        irreps % (cell(1,\*) of `.SubRep`): Equivalent irreducible subrepresentations in this isotypic component
        modelIrrep % (`.Rep`): Irreducible representation identical to the subrepresentations present in this isotypic component
    end

    methods (Static)

        function iso = fromTrivialSubRep(trivial)
        % Builds an isotypic component from the trivial subrepresentation
        %
        % The trivial subrepresentation must be complete, i.e. represent the space of all vectors that are fixed by the
        % action of the representation.
        %
        % Args:
        %   trivial (`.SubRep`): Maximal trivial subrepresentation of ``parent``
        %
        % Returns:
        %   `.Isotypic`: The corresponding trivial isotypic component
            parent = trivial.parent;
            dT = trivial.dimension;
            if dT == 0
                group = trivial.group;
                irreps = cell(1, 0);
                modelIrrep = group.trivialRep(trivial.field, 1);
                iso = replab.Isotypic.fromIrreps(parent, irreps, modelIrrep, 'irrepsAreBiorthogonal', true, 'irrepsAreHarmonized', true);
            else
                irreps = cell(1, dT);
                I = trivial.injection_internal;
                P = trivial.projection_internal;
                for i = 1:dT
                    irreps{i} = parent.subRep(I(:,i), 'projection', P(i,:), 'isIrreducible', true, 'frobeniusSchurIndicator', 1, 'trivialDimension', 1, 'isUnitary', true);
                end
                iso = replab.Isotypic.fromIrreps(parent, irreps, [], 'irrepsAreBiorthogonal', true, 'irrepsAreHarmonized', true);
            end
        end

        function iso = fromIrreps(parent, irreps, modelIrrep, varargin)
        % Creates an isotypic component from linearly independent equivalent irreducible subrepresentations
        %
        % Args:
        %   parent (`.Rep`): Parent representation
        %   irreps (cell(1,\*) of `.SubRep`): Irreducible equivalent subrepresentations of ``parent``
        %   modelIrrep (`.Rep` or ``[]``): Model of the irreducible representation in this component (can be empty, and is then set to ``irreps{1}``)
        %
        % Keyword Args:
        %   irrepsAreBiorthogonal (logical, optional): True if the irreps injection/projection are biorthogonal, default: false
        %   irrepsAreHarmonized (logical, optional): True if the irreps are in the same basis as ``modelIrrep``, default: false
        %
        % Returns:
        %   `.Isotypic`: Isotypic component
            if isempty(irreps)
                assert(~isempty(modelIrrep), 'modelIrrep must be provided if the isotypic component is empty');
                iso = replab.Isotypic(parent, cell(1, 0), modelIrrep, zeros(0, parent.dimension));
                return
            end
            n = length(irreps);
            % from now on, assume that the isotypic component is not empty
            args = struct('irrepsAreBiorthogonal', false, 'irrepsAreHarmonized', false);
            args = replab.util.populateStruct(args, varargin);
            if isempty(modelIrrep)
                harmonizeFrom = 2; % if harmonization is needed, skip first irrep
                modelIrrep = irreps{1};
            else
                harmonizeFrom = 1;
            end
            if ~args.irrepsAreHarmonized
                genModel = replab.rep.GenSubRep.fromRep(modelIrrep);
                for i = harmonizeFrom:n
                    gen = replab.rep.GenSubRep.fromSubRep(irreps{i});
                    gen1 = gen.harmonize(genModel);
                    sub1 = gen1.toSubRep;
                    irreps{i} = irreps{i}.withUpdatedMaps(sub1.injection, sub1.projection);
                end
            end
            if ~args.irrepsAreBiorthogonal || ~args.irrepsAreHarmonized
                injections = cellfun(@(irrep) irrep.injection_internal, irreps, 'uniform', 0);
                I = horzcat(injections{:});
                P = replab.rep.findProjection_largeScale(parent, I, 5, replab.rep.Tolerances, [], []);
            else
                projections = cellfun(@(irrep) irrep.projection_internal, irreps, 'uniform', 0);
                P = vertcat(projections{:});
            end
            iso = replab.Isotypic(parent, irreps, modelIrrep, P);
        end

    end

    methods

        function self = Isotypic(parent, irreps, modelIrrep, projection)
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
        %   modelIrrep (`+replab.Rep`): Canonical model of the irreducible representation contained in this component
        %   projection (double(\*,\*), may be sparse): Embedding map of the isotypic component
            m = length(irreps);
            assert(isa(parent, 'replab.Rep'));
            assert(all(cellfun(@(irrep) isa(irrep, 'replab.SubRep'), irreps)));
            assert(all(cellfun(@(irrep) irrep.isIrreducibleAndCanonical, irreps)));
            injections = cellfun(@(irrep) irrep.injection_internal, irreps, 'uniform', 0);
            if m == 0
                injection = zeros(parent.dimension, 0);
            else
                injection = horzcat(injections{:});
            end
            self = self@replab.SubRep(parent, injection, projection, ...
                                      'isUnitary', modelIrrep.isUnitary, ...
                                      'isIrreducible', m <= 1, ...
                                      'trivialDimension', modelIrrep.trivialDimension * m, ...
                                      'divisionAlgebraName', modelIrrep.divisionAlgebraName);
            self.irreps = irreps;
            self.modelIrrep = modelIrrep;
        end

    end

    methods

        function m = multiplicity(self)
        % Number of equivalent irreducible representations in this isotypic component
        %
        % Returns:
        %   integer: Multiplicity
            m = length(self.irreps);
        end

        function d = irrepDimension(self)
        % Returns the dimension of a single irreducible representation contained in this component
        %
        % Returns:
        %   integer: Irrep dimension
            d = self.modelIrrep.dimension;
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

    end

    methods % Equivariant spaces

        function E = isotypicEquivariantFrom(self, repC, varargin)
        % Returns the space of equivariant linear maps from another isotypic component to this isotypic component
        %
        % The equivariant vector space contains the matrices X such that
        %
        % ``X * repC.image(g) = self.image(g) * X``
        %
        % Both isotypic components must be harmonized.
        %
        % Args:
        %   repC (`+replab.Isotypic`): Isotypic component, representation on the source/column space
        %
        % Keyword Args:
        %   special ('commutant', 'hermitian', 'trivialRows', 'trivialCols' or '', optional): Special structure if applicable, see `.Equivariant`, default: ''
        %   type ('exact', 'double' or 'double/sparse', optional): Whether to obtain an exact equivariant space, default 'double' ('double' and 'double/sparse' are equivalent)
        %
        % Returns:
        %   `+replab.IsotypicEquivariant` or ``[]``: The equivariant vector space, or ``[]`` if the space has dimension zero or contains only the zero matrix
            E = replab.IsotypicEquivariant.make(self, repC, varargin{:});
        end

        function E = isotypicEquivariantTo(self, repR, varargin)
        % Returns the space of equivariant linear maps from this isotypic component to another isotypic component
        %
        % The equivariant vector space contains the matrices X such that
        %
        % ``X * self.image(g) = repR.image(g) * X``
        %
        % Both isotypic components must be harmonized.
        %
        % Args:
        %   repR (`+replab.Isotypic`): Isotypic component, representation on the target/row space
        %
        % Keyword Args:
        %   special ('commutant', 'hermitian', 'trivialRows', 'trivialCols' or '', optional): Special structure if applicable, see `.Equivariant`, default: ''
        %   type ('exact', 'double' or 'double/sparse', optional): Whether to obtain an exact equivariant space, default 'double' ('double' and 'double/sparse' are equivalent)
        %
        % Returns:
        %   `+replab.IsotypicEquivariant` or ``[]``: The equivariant vector space, or ``[]`` if the space has dimension zero or contains only the zero matrix
            E = replab.IsotypicEquivariant.make(repR, self, varargin{:});
        end

    end

    methods (Access = protected) % Implementations

        function rho = image_double_sparse(self, g)
            if self.multiplicity > 1
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
            if self.multiplicity > 1
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
            s = 'Isotypic component';
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
            if self.inCache('trivialDimension')
                if self.trivialDimension == self.dimension
                    s = [s ' (trivial)'];
                elseif self.trivialDimension == 0
                    s = [s ' (nontrivial)'];
                end
            end
        end

        % Obj

        function l = laws(self)
            l = replab.laws.IsotypicLaws(self);
        end

        % Rep

        function c = commutant(self, type)
            if nargin < 2 || isempty(type) || strcmp(type, 'double/sparse')
                type = 'double';
            end
            c = self.cached(['commutant_' type], @() self.isotypicEquivariantFrom(self,  'special', 'commutant', 'type', type));
        end

        % SubRep

        function iso = refine(self, varargin)
        % Refines this isotypic component
        %
        % Assumes that the isotypic component is already close to an exact isotypic component, and refines its subspace
        % by an iterative procedure applied on its `.injection` and `.projection` maps.
        %
        % Keyword Args:
        %   numNonImproving (integer, optional): Number of non-improving steps before stopping the large-scale algorithm, default ``20``
        %   nSamples (integer, optional): Number of samples to use in the large-scale version of the algorithm, default ``5``
        %   maxIterations (integer, optional): Maximum number of (outer) iterations, default ``1000``
        %
        % Returns:
        %   `.Isotypic`: Isotypic component with refined subspace (injection/projection maps)
            if self.multiplicity == 0
                iso1 = self;
                return
            end
            replab.msg(1, 'Refining isotypic component');
            args = struct('numNonImproving', 20, 'largeScale', self.parent.dimension > 1000, 'nSamples', 5, 'nInnerIterations', 10, 'maxIterations', 1000);
            args = replab.util.populateStruct(args, varargin);
            D = self.parent.dimension;
            m = self.multiplicity;
            irreps = self.irreps;
            parent = self.parent;
            modelIrrep = self.modelIrrep;
            if ~modelIrrep.isExact
                modelIrrep = modelIrrep.refine;
            end
            genM = replab.rep.GenSubRep.fromSubRep(modelIrrep);
            if parent.isUnitary && modelIrrep.isUnitary
                Qo = zeros(D, 0);
                for i = 1:m
                    irr = irreps{i};
                    gen = replab.rep.GenSubRep.fromSubRep(irr);
                    gen = replab.rep.harmonize_unitary_largeScale(gen, genM, 20, 5, 100, Qo);
                    sub = gen.toSubRep;
                    irr = irr.withUpdatedMaps(sub.injection, sub.projection);
                    Qo = horzcat(Qo, irr.injection);
                    irreps{i} = irr;
                end
                iso = replab.Isotypic(parent, irreps, modelIrrep, Qo');
            else
                I = self.injection;
                P = self.projection;
                Io = zeros(D, 0);
                Po = zeros(0, D);
                for i = 1:m
                    irr = irreps{i};
                    range = self.irrepRange(i);
                    irr = irr.withUpdatedMaps(I(:, range), P(range, :));
                    gen = replab.rep.GenSubRep.fromSubRep(irr);
                    gen = replab.rep.harmonize_nonUnitary_largeScale(gen, genM, 20, 5, 100, Io, Po);
                    sub = gen.toSubRep;
                    irr = irr.withUpdatedMaps(sub.injection, sub.projection);
                    Io = [Io, irr.injection];
                    Po = [Po; irr.projection];
                    irreps{i} = irr;
                end
                % apply one Newton-Raphson iteration for increased precision
                Po = 2*Po - Po*Io*Po;
                iso = replab.Isotypic(parent, irreps, modelIrrep, Po);
            end
        end

    end

end
