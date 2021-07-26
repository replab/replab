classdef FiniteGroupType < replab.Group
% Describes a type of finite groups

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

end
