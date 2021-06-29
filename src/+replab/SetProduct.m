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

        function P = fromFiniteSet(finiteSet, identityFirst)
        % Constructs a SetProduct decomposition from a single set
        %
        % No efficiency gains are expected.
        %
        % Args:
        %   finiteSet (`.FiniteSet`): Finite set to describe using a (trivial) `.SetProduct`
        %   identityFirst (logical, optional): Whether to move the identity in the first place, if present, default: false
        %
        % Returns:
        %   `.SetProduct`: The set product
            if nargin < 2 || isempty(identityFirst)
                identityFirst = false;
            end
            c = double(finiteSet.nElements);
            S = finiteSet.elements.toCell;
            group = finiteSet.type;
            if identityFirst
                ind = double(finiteSet.elements.find(group.identity));
                if ind > 0 && ind ~= 1
                    S([1 ind]) = S([ind 1]);
                end
            end
            P = replab.SetProduct(group, {S}, identityFirst);
        end

    end

end
