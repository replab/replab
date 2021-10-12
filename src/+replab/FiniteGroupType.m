classdef FiniteGroupType < replab.Group & replab.TotalOrder
% Describes a type of finite groups
%
% Examples of types include:
%
% - permutations of a given domain ``1..n`` (the set of such permutations is finite),
% - invertible matrices of a given size with coefficients in the cyclotomic field
%   (the set of such matrices is not finite, but we assume that any explicitly constructed group
%    is finite),
% - direct products of finite group types.
%
% On top of providing the group structure, this class provides methods to construct finite groups and subgroups.
% It also defines a more or less arbitrary total order of elements, which is used to define representatives of
% various `.FiniteSet` instances, such as conjugacy classes and cosets.

    methods

        function G = groupWithGenerators(self, generators, varargin)
        % Constructs a group with the given generators
        %
        % Args:
        %   generators (cell(1,\*) of group elements): Group generators
        %
        % Keyword Args:
        %   generatorNames (cell(1,\*) of charstring): Names of the generators
        %   order (vpi or integer): Group order
        %   relators (cell(1,\*) of charstring): Relators
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
        %   group (`.GenericFiniteGroup`): Finite group whose subgroup we are constructing
        %   generators (cell(1,\*) of group elements): Subgroup generators
        %
        % Keyword Args:
        %   generatorNames (cell(1,\*) of charstring): Names of the generators
        %   order (vpi or integer): Group order
        %   relators (cell(1,\*) of charstring): Relators
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
