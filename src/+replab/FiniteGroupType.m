classdef FiniteGroupType < replab.Group & replab.TotalOrder
% Describes a type of finite groups
%
% It is a possibly infinite group, and this object is able to construct finite subgroups of itself.
%
% It also defines a (somewhat arbitrarily chosen) total order used, for example, to define
% minimal canonical representatives of conjugacy classes and cosets.
%
% Examples of group types include:
%
% - permutations of a given domain ``1..n`` (the set of such permutations is finite),
% - invertible matrices of a given size with coefficients in the cyclotomic field
%   (the set of such matrices is not finite, but we assume that any explicitly constructed group
%    is finite),
% - direct products of finite group types.
%
% Note:
%   For implementations purposes, there are several families of group types in RepLAB.
%
%   * The first is `+replab.PermutationGroupType` where computational group theory algorithms
%     are implemented and used directly.
%
%   * The second is `+replab.+gen.FiniteGroupType`, in which the constructed objects delegate all
%     computations to objects where the algorithms are actually implemented, through a
%     "nice isomorphism".
%
%   For now, those two families are the only ones implemented. Natural extensions would be:
%
%   * Implementing algorithms for solvable groups through the polycyclic presentation as in
%     e.g. GAP System.
%
%   * Implementing black box group algorithms, either through a naive enumeration approach
%     (to validate existing algorithms), or through recent algorithms proposed for black box groups.
%
%   * Implementing matrix groups through a BSGS approach using the action of matrices on vectors;
%     possibly also through recent algorithms based on composition trees.

    methods

        function G = groupWithGenerators(self, generators, varargin)
        % Constructs a group with the given generators
        %
        % Args:
        %   generators (cell(1,\*) of group elements): Group generators
        %
        % Keyword Args:
        %   generatorNames (cell(1,\*) of charstring, optional): Names of the generators
        %   order (vpi or integer, optional): Group order
        %   relators (cell(1,\*) of charstring or integer(1,\*), optional): Relators
        %
        % Returns:
        %   `.FiniteGroup`: Constructed group
            error('Abstract');
        end

        function l = isSameTypeAs(self, otherType)
        % Tests whether this type is the same as another type
        %
        % Args:
        %   otherType (`.FiniteGroupType`): Other type to test with
        %
        % Returns:
        %   logical: True if the types are the same
            error('Abstract');
        end

        function G = subgroupWithGenerators(self, group, generators, varargin)
        % Constructs a subgroup of a given group with the given parameters
        %
        % Args:
        %   group (`.FiniteGroup`): Finite group whose subgroup we are constructing
        %   generators (cell(1,\*) of group elements): Subgroup generators
        %
        % Keyword Args:
        %   generatorNames (cell(1,\*) of charstring): Names of the generators
        %   order (vpi or integer): Group order
        %   relators (cell(1,\*) of charstring or integer(1,\*)): Relators
        %
        % Returns:
        %   `.FiniteGroup`: Constructed group
            G = self.groupWithGenerators(generators, varargin{:});
        end

    end

    methods % Implementations

        % Domain

        function l = laws(self)
            l = replab.laws.FiniteGroupTypeLaws(self);
        end

    end
    
end
