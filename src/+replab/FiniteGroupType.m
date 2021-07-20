classdef FiniteGroupType < replab.Group

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

    end

end
