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
%
% When ``repR`` is unitary, the ``commutant`` and ``hermitian`` cases are identical.

    properties (SetAccess = protected)
        field % ({'R', 'C'}): Field of the vector space real (R) or complex x(C)
        nR % (integer): Row size
        nC % (integer): Column size
        group % (`+replab.CompactGroup`): Group being represented
        repR % (`+replab.Rep`): Representation of row space
        repC % (`+replab.Rep`): Representation of column space
        special % ({'hermitian', 'commutant', []}): Whether the equivariant space has special structure
    end

    properties (Access = protected)
        domain % (`+replab.Domain`): Domain, real or complex matrices
        cachedSamples_ % (struct of cell(1,*) of domain elements): Samples
        cachedErrors_ % (struct of double(1,*)): Error information
    end

    methods

        %% Abstract

        function [X1 err] = project(self, X)
        % Projects any ``nR x nC`` matrix in the equivariant subspace
        %
        % Performs the integration
        %
        % `` X1 = int{g in G} dg rhoR.image(g) * X * rhoC.inverseImage(g) ``
        %
        % Args:
        %   X (double(*,*)): Matrix to project
        %
        % Returns
        % -------
        % X1:
        %   double(*,*): Projected matrix
        % err:
        %   double: Estimation of the numerical error, expressed as the distance of the returned ``X1`` to
        %           the invariant subspace in Frobenius norm or ``NaN`` if no such estimation can be performed
            error('Abstract');
        end

    end

    methods

        function self = Equivariant(repC, repR, special)
        % Constructor; use `+replab.makeEquivariant` in user code
        %
        % That function selects an optimized implementation depending on the use case.
            self.repR = repR;
            self.nR = repR.dimension;
            self.repC = repC;
            self.nC = repC.dimension;
            assert(isequal(repR.field, repC.field), ...
                   'Both representations must have be defined on the same field');
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
        %   context (`+replab.+equivariant.Context`): Context in which samples are cached
        %   i (double): 1-based index of the sample
        %
        % Returns
        % -------
        % X:
        %   double(*,*): A sample from this equivariant space
        % err:
        %   double: Estimation of the numerical error, expressed as the distance of the returned ``X`` to
        %           the invariant subspace in Frobenius norm
            [X err] = self.project(self.domain.sample);
        end

        function clearCache(self, context)
        % Clears the samples cached for the given context
        %
        % Args:
        %   context (`+replab.+equivariant.Context`): Context to clear
            self.cachedSamples_ = rmfield(self.cachedSamples_, context.id);
            self.cachedErrors_ = rmfield(self.cachedErrors_, context.id);
        end

        function [X err] = sampleInContext(self, context, ind)
        % Returns an approximate sample from this equivariant space along with estimated numerical error
        %
        % The samples are cached in a context.
        %
        % Args:
        %   context (`+replab.+equivariant.Context`): Context in which samples are cached
        %   ind (double): 1-based index of the sample
        %
        % Returns
        % -------
        % X:
        %   double(*,*): A sample from this equivariant space
        % err:
        %   double: Estimation of the numerical error, expressed as the distance of the returned ``X`` to
        %           the invariant subspace in Frobenius norm
            assert(~context.closed);
            id = context.id;
            if ~isfield(self.cachedSamples_, id)
                context.register(self);
                self.cachedSamples_.(id) = cell(1, 0);
                self.cachedErrors_.(id) = zeros(1, 0);
            end
            n = length(self.cachedSamples_.(id));
            samples = self.cachedSamples_.(id);
            errors = self.cachedErrors_.(id);
            if ind > n
                for i = n+1:ind
                    [X err] = self.sampleWithError;
                    samples{1, i} = X;
                    errors{1, i} = err;
                end
                self.cachedSamples_.(id) = samples;
                self.cachedErrors_.(id) = errors;
            end
            X = samples{ind};
            err = errors(ind);
        end

        function E1 = subEquivariant(self, subC, subR, special)
        % Constructs a invariant subspace of an equivariant space
        %
        % Args:
        %   subC (`+replab.SubRep`): A subrepresentation of ``self.repC``
        %   subR (`+replab.SubRep`): A subrepresentation of ``self.repR``
        %   special (charstring): Whether the equivariant subspace has special structure
            assert(isa(subC, 'replab.SubRep'));
            assert(isa(subR, 'replab.SubRep'));
            assert(subC.parent == self.repC);
            assert(subR.parent == self.repR);
            E1 = replab.equivariant.ForSubReps(subC, subR, special, self);
        end

        %% Str methods

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

        %% Domain methods

        function b = eqv(self, X, Y)
            b = self.domain.eqv(X, Y);
        end

        function X = sample(self)
            X = self.sampleWithError;
        end

    end

end
