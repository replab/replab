classdef FiniteGroupType < replab.FiniteGroupType
% Finite group type that delegates computations through an isomorphism
%
% Every object created in this group type has an isomorphic counterpart, its "nice"
% counterpart, stored in `+replab.+gen.FiniteSet.nice`.
%
% The classes defined in the `+replab.+gen` package solve two challenges:
%
% * Construct an isomorphism for any generating set of type elements.
%   Different strategies are used depending on the type of group considered.
%
% * Implement the different finite objects in RepLAB through this isomorphism scheme.
%   The implementations in `+replab.+gen.FiniteSet`, `+replab.+gen.LeftCoset`, ...
%   are straightforward.
%
% This class has an important subclass, `+replab.+gen.StaticFiniteGroupType`, used for
% group types where the same isomorphism can be reused for all generating sets of elements.
% For example, all signed permutation groups acting on the same domain can be mapped to
% the permutation group (which acts on a domain of doubled size).
%
% Other group types, such as the one for matrix groups with cyclotomic coefficients
% (see `+replab.+matrix.FiniteGroupType` ), construct tailor-made isomorphisms for
% each finite matrix group.

    methods (Access = protected)

        function G = makeGenericGroup(self, generators, nice, niceIsomorphism)
        % Creates a generic group of this type
        %
        % Args:
        %   generators (cell(1,\*) of elements of this type): Group generators
        %   nice (`+replab.FiniteGroup`): Subgroup of the isomorphism target whose generators are in 1-to-1 correspondance with ``generators``
        %   niceIsomorphism (`+replab.+gen.NiceIsomorphism`): Isomorphism whose source contains all the generators
            G = replab.gen.FiniteGroup(self, generators, 'nice', nice, 'niceIsomorphism', niceIsomorphism);
        end

    end

    methods

        function iso = constructNiceIsomorphism(self, elements)
        % Constructs a generic isomorphism for the given elements
        %
        % It is not guaranteed that the type of the isomorphism target is always the same.
        %
        % Args:
        %   elements (cell(1,\*) of group type elements): Elements
        %
        % Returns:
        %   `.FiniteIsomorphism`: An order preserving isomorphism whose source contains all given elements
            error('Abstract');
        end

    end

    methods % Implementations

        function G = groupWithGenerators(self, generators, varargin)
            assert(all(ismember(varargin(1:2:end), {'generatorNames', 'order', 'relators'})));
            niceIso = self.constructNiceIsomorphism(generators);
            target = niceIso.target;
            targetGenerators = cellfun(@(g) niceIso.imageElement(g), generators, 'uniform', 0);
            nice = target.subgroup(targetGenerators, varargin{:});
            G = self.makeGenericGroup(generators, nice, niceIso);
        end

    end

end
