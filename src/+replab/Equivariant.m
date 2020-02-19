classdef Equivariant < replab.Domain
% Describes a vector space of group-equivariant matrices
%
% Let ``repR`` and ``repC`` be two representations of the same group ``G``.
%
% This describes the set of matrices ``X`` such that ``repR.image(g) * X = X * repC.image(g)``
%
% See Proposition 4 of
% J.-P. Serre, Linear Representations of Finite Groups (Springer, 1977).
%
% There are two special cases of equivariant spaces.
%
% - ``commutant`` equivariant spaces have ``repR == repC``,
% - ``hermitian`` equivariant spaces have ``repR == repC.conjugate.dual``.
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
        parent % (`+replab.Domain`): Parent domain, real or complex matrices
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
        %           the invariant subspace in Frobenius norm
            error('Abstract');
        end

    end

    methods

        function [X err] = sampleWithError(self)
        % Returns an approximate sample from this equivariant space along with estimated numerical error
        %
        % Returns
        % -------
        % X:
        %   double(*,*): A sample from this equivariant space
        % err:
        %   double: Estimation of the numerical error, expressed as the distance of the returned ``X`` to
        %           the invariant subspace in Frobenius norm
            [X err] = self.project(self.parent.sample);
        end

        function self = Equivariant(repR, repC, special)
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
            self.parent = replab.domain.Matrices(self.field, self.nR, self.nC);
            self.special = special;
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
            b = self.parent.eqv(X, Y);
        end

        function X = sample(self)
            X = self.sampleWithError;
        end

    end

end
