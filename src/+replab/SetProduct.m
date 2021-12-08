classdef SetProduct < replab.Domain
% Describes a multiset of elements obtained by reduction of the elements of a cartesian product
%
% Let ``T1, ..., Tn`` be finite subsets of a monoid ``G``. The set product of such ``T`` is the multiset
% of elements of the form ``x = t1 t2 ... tn`` iterating over each ``ti`` in each ``Ti``. The binary operation
% is given by the composition of ``G``.
%
% Thus, the multiset has ``|T1| * |T2| * ... * |Tn|`` elements, possibly with repetition.
%
% RepLAB uses `.SetProduct` in two areas:
%
% - to describe a unique decomposition of a finite group, in which case the multiset is actually a set whose
%   cardinality is the group order; we additionally require that the first element of each ``T{i}`` is the
%   identity;
%
% - to describe the left cosets of the connected component in a compact Lie group.

    properties (SetAccess = protected)
        monoid % (`+replab.Monoid`): Monoid providing the binary operation
        identityFirst % (logical): True if each element `.T` contains the identity as its first element
        sets % (cell(1,\*) of cell(1,\*) of elements): Stores the sets Ti as ``{T1 T2 ... Tn}``
    end

    methods

        function self = SetProduct(monoid, sets, identityFirst)
            self.monoid = monoid;
            self.sets = sets;
            self.identityFirst = identityFirst;
        end

        function S = imap(self, isomorphism)
        % Maps this multiset under an isomorphism
        %
        % Args:
        %   isomorphism (`.Isomorphism`): Isomorphism with its containing all elements of this `.SetProduct`
        %
        % Returns:
        %   `.SetProduct`: Multiset with elements in ``isomorphism.target``
            sets1 = cellfun(@(s) cellfun(@(el) isomorphism.imageElement(el), s, 'uniform', 0), self.sets, 'uniform', 0);
            S = replab.SetProduct(isomorphism.target, sets1, self.identityFirst);
        end

    end

    methods % Implementations

        % Obj

        function l = laws(self, finiteSet)
            if nargin == 1
                l = replab.laws.SetProductLaws(self);
            else
                l = replab.laws.SetProductLaws(self, finiteSet);
            end
        end

        % Domain

        function l = eqv(self, x, y)
            l = self.monoid.eqv(x, y);
        end

        function s = sample(self)
            s = self.monoid.identity;
            for i = 1:length(self.sets)
                S = self.sets{i};
                s = self.monoid.compose(s, S{randi(length(S))});
            end
        end

    end

    methods (Static)

        function P = identity(monoid)
        % Constructs a SetProduct containing only the identity of a monoid
        %
        % Args:
        %   monoid (`.Monoid`): Monoid to take the identity from
        %
        % Returns:
        %   `.SetProduct`: A singleton SetProduct
            P = replab.SetProduct(monoid, {{monoid.identity}}, true);
        end

    end

end
