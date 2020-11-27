classdef Equivariant < replab.Domain
% Describes a vector space of group-equivariant matrices
%
% Let ``repC`` and ``repR`` be two representations of the same group ``G``.
%
% This describes the set of matrices ``X`` such that ``repR.image(g) * X = X * repC.image(g)``
%
% See Proposition 4 of
% J.-P. Serre, Linear Representations of Finite Groups (Springer, 1977).
%
% There are two special cases of equivariant spaces.
%
% - ``commutant`` equivariant spaces have ``repC == repR``,
% - ``hermitian`` equivariant spaces have ``repC == repR.conjugate.dual``.
% - ``trivial`` equivariant spaces have ``repC == repR.group.trivialRep(field, repR.dimension)``
%
% When ``repR`` is unitary, the ``commutant`` and ``hermitian`` cases are identical.

    properties (SetAccess = protected)
        field % ({'R', 'C'}): Field of the vector space real (R) or complex x(C)
        nR % (integer): Row size
        nC % (integer): Column size
        group % (`+replab.CompactGroup`): Group being represented
        repR % (`+replab.Rep`): Representation of row space
        repC % (`+replab.Rep`): Representation of column space
        special % ({'hermitian', 'commutant', 'trivialRows', 'trivialCols', []}): Whether the equivariant space has special structure
    end

    properties (Access = protected)
        domain % (`+replab.Domain`): Domain, real or complex matrices
        cachedSamples_ % (struct of cell(1,\*) of domain elements): Samples
        cachedErrors_ % (struct of double(1,\*)): Error information
    end

    methods % Projection

        function [X1 err] = project(self, X, type)
        % Projects any ``nR x nC`` matrix in the equivariant subspace
        %
        % Performs the integration
        %
        % `` X1 = int{g in G} dg rhoR.image(g) * X * rhoC.inverseImage(g) ``
        %
        % Args:
        %   X (double(\*,\*) or `.cyclotomic`(\*,\*), may be sparse): Matrix to project
        %
        % Returns
        % -------
        %   X1: double(\*,\*) or `.cyclotomic`(\*,\*)
        %     Projected matrix
        %   err: double
        %     Estimation of the numerical error, expressed as the distance of the returned ``X1`` to the invariant subspace in Frobenius norm
            error('Abstract');
        end

    end

    methods

        function self = Equivariant(repC, repR, special)
        % Constructor; use `+replab.Equivariant.make` in user code
        %
        % That function selects an optimized implementation depending on the use case.
            self.repR = repR;
            self.nR = repR.dimension;
            self.repC = repC;
            self.nC = repC.dimension;
            assert(isequal(repR.field, repC.field), ...
                   'Both representations must have be defined on the same field');
            assert(isempty(special) || ismember(special, {'hermitian', 'commutant', 'trivialRows', 'trivialCols'}));
            self.field = repR.field;
            assert(repR.group == repC.group, ...
                   'Both representations must be defined on the same group');
            self.group = repR.group;
            self.domain = replab.domain.Matrices(self.field, self.nR, self.nC);
            self.special = special;
            self.cachedSamples_ = struct;
            self.cachedErrors_ = struct;
        end

        function [X err] = sampleWithError(self)
        % Returns an approximate sample from this equivariant space along with estimated numerical error
        %
        % The samples are cached in a context.
        %
        % Args:
        %   context (`+replab.Context`): Context in which samples are cached
        %   i (double): 1-based index of the sample
        %
        % Returns
        % -------
        % X:
        %   double(\*,\*): A sample from this equivariant space
        % err:
        %   double: Estimation of the numerical error, expressed as the distance of the returned ``X`` to
        %           the invariant subspace in Frobenius norm
            [X err] = self.project(self.domain.sample);
        end
% $$$
% $$$         function clearCache(self, context)
% $$$         % Clears the samples cached for the given context
% $$$         %
% $$$         % Args:
% $$$         %   context (`+replab.Context`): Context to clear
% $$$             self.cachedSamples_ = rmfield(self.cachedSamples_, context.id);
% $$$             self.cachedErrors_ = rmfield(self.cachedErrors_, context.id);
% $$$         end
% $$$
% $$$         function [X err] = sampleInContext(self, context, ind)
% $$$         % Returns an approximate sample from this equivariant space along with estimated numerical error
% $$$         %
% $$$         % The samples are cached in a context.
% $$$         %
% $$$         % Args:
% $$$         %   context (`+replab.Context`): Context in which samples are cached
% $$$         %   ind (double): 1-based index of the sample
% $$$         %
% $$$         % Returns
% $$$         % -------
% $$$         % X:
% $$$         %   double(\*,\*): A sample from this equivariant space
% $$$         % err:
% $$$         %   double: Estimation of the numerical error, expressed as the distance of the returned ``X`` to
% $$$         %           the invariant subspace in Frobenius norm
% $$$             assert(~context.closed);
% $$$             id = context.id;
% $$$             if ~isfield(self.cachedSamples_, id)
% $$$                 context.register(self);
% $$$                 self.cachedSamples_.(id) = cell(1, 0);
% $$$                 self.cachedErrors_.(id) = zeros(1, 0);
% $$$             end
% $$$             n = length(self.cachedSamples_.(id));
% $$$             samples = self.cachedSamples_.(id);
% $$$             errors = self.cachedErrors_.(id);
% $$$             if ind > n
% $$$                 for i = n+1:ind
% $$$                     [X err] = self.sampleWithError;
% $$$                     samples{1, i} = X;
% $$$                     errors{1, i} = err;
% $$$                 end
% $$$                 self.cachedSamples_.(id) = samples;
% $$$                 self.cachedErrors_.(id) = errors;
% $$$             end
% $$$             X = samples{ind};
% $$$             err = errors(ind);
% $$$         end

    end

    methods % Implementations

% $$$         function E1 = subEquivariant(self, subC, subR, special)
% $$$         % Constructs a invariant subspace of an equivariant space
% $$$         %
% $$$         % Args:
% $$$         %   subC (`+replab.SubRep`): A subrepresentation of ``self.repC``
% $$$         %   subR (`+replab.SubRep`): A subrepresentation of ``self.repR``
% $$$         %   special (charstring): Whether the equivariant subspace has special structure
% $$$             assert(isa(subC, 'replab.SubRep'));
% $$$             assert(isa(subR, 'replab.SubRep'));
% $$$             assert(subC.parent == self.repC);
% $$$             assert(subR.parent == self.repR);
% $$$             E1 = replab.equi.ForSubReps(subC, subR, special, self);
% $$$         end

        % Str

        function s = headerStr(self)
            if isequal(self.special, 'hermitian')
                s = sprintf('%d x %d invariant Hermitian matrices over %s', ...
                            self.nR, self.nC, self.field);
            elseif isequal(self.special, 'commutant')
                s = sprintf('%d x %d commutant matrices over %s', ...
                            self.nR, self.nC, self.field);
            else
                s = sprintf('%d x %d equivariant matrices over %s', ...
                            self.nR, self.nC, self.field);
            end
        end

        % Domain

        function b = eqv(self, X, Y)
            b = self.domain.eqv(X, Y);
        end

        function X = sample(self)
            X = self.sampleWithError;
        end

    end

    methods (Static)

        function E = make(repC, repR, special)
        % Returns the space of equivariant linear maps between two representations
        %
        % The equivariant vector space contains the matrices X such that
        %
        % ``repC.image(g) * X = X * repR.image(g)``
        %
        % Args:
        %   repC (`+replab.Rep`): Representation on the source/column space
        %   repR (`+replab.Rep`): Representation on the target/row space
        %   special (charstring): Special structure see help on `+replab.Equivariant.special`
        %
        % Returns:
        %   `+replab.Equivariant`: The equivariant vector space
            if isa(repR.group, 'replab.FiniteGroup')
                E = replab.equi.Equivariant_forFiniteGroup(repC, repR, special);
            else
                error('Unimplemented');
            end
        end

    end

end
