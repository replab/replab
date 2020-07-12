classdef FreeGroup < replab.Group
% Describes a free group
%
% Example:
%   >>> [F, a, x] = replab.FreeGroup.of('a', 'x');
%   >>> a * inv(a) == F.identity
%       1

    properties (Access = protected)
        groupId % (integer): Unique group id
    end

    properties (SetAccess = protected)
        generatorNames % (cell(1,\*) of charstring): Names of the generators
        generators % (cell(1,\*) of `.FreeGroupWord`): Generators
    end

    methods % Constructor

        function self = FreeGroup(generatorNames)
        % Creates a free group with the given generator names
        %
        % Example:
        %   >>> replab.FreeGroup({'x', 'a'})
        %     Free group < a, x >
        %     generatorNames: {'a', 'x'}
        %           identity: 1
        %
        % Args:
        %   generatorNames (cell(1,\*) of charstring): Generator names
        %
        % Returns:
        %   `.FreeGroup`: The constructed free group
            self.identity = replab.FreeGroupWord.empty(self);
            self.groupId = replab.globals.nextUniqueId;
            self.generatorNames = generatorNames;
            self.generators = arrayfun(@(i) replab.FreeGroupWord.make(self, i), 1:length(generatorNames), 'uniform', 0);
        end

    end

    methods (Static) % FreeGroup construction

        function [F, varargout] = of(varargin)
        % Creates a free group with the given generator names
        %
        % Example:
        %   >>> [F, x, a] = replab.FreeGroup.of({'x', 'a'});
        %   >>> x * x
        %       x^2
        %
        % Args:
        %   varargin: Generator names, given as charstrings
        %
        % Returns:
        %   `.FreeGroup`: The constructed free group
            F = replab.FreeGroup(varargin);
            for i = 1:F.rank
                varargout{i} = F.generator(i);
            end
        end

    end

    methods % Implementations

        function res = eq(self, rhs)
            res = (self.groupId == rhs.groupId);
        end

        function res = ne(self, rhs)
            res = self.groupId ~= rhs.groupId;
        end

        function s = headerStr(self)
            s = ['Free group < ' strjoin(self.generatorNames, ', ') ' >'];
        end

        function x = sample(self)
            l = 10;
            letters = randi([-self.rank self.rank], 1, l);
            letters = letters(letters ~= 0);
            x = replab.FreeGroupWord.make(self, letters);
        end

        function res = eqv(self, x, y)
            res = (x == y);
        end

        function z = compose(self, x, y)
            z = x*y;
        end

        function z = inverse(self, x)
            z = inv(x);
        end

    end

    methods % Group construction

        function sub = mrdivide(self, relators)
        % Constructs the quotient of this group by relators
        %
        % RepLAB assumes that the group thus described is finite and has small order.
        %
        % Args:
        %   relators (cell(1,\*) of `.FreeGroupWord`): List of relators
        %
        % Returns:
        %   `+replab.FiniteFPGroup`: The constructed finite group
            sub = replab.FiniteFPGroup(self, relators);
        end

    end

    methods % Word operations

        function w = parse(self, str)
        % Constructs a word from a string
        %
        % Example:
        %   >>> [F x y] = replab.FreeGroup.of('x', 'y');
        %   >>> F.parse('x (x y)^2 / x')
        %       x^2 y x y x^-1
        %
        % Args:
        %   str (charstring): String describing a word
        %
        % Returns:
        %   `.FreeGroupWord`: The parsed word
        %
        % Raises:
        %   An error if the string is malformed
            w = replab.FreeGroupWord.parse(self, str);
        end

        function w = word(self, letters)
        % Constructs a word from a sequence of letters
        %
        % The constructed word is reduced.
        %
        % Example:
        %   >>> [F x y] = replab.FreeGroup.of('x', 'y');
        %   >>> F.word([1 2 -2 1])
        %       x^2
        %
        % Args:
        %   letters (integer(1,\*)): Sequence of letters with indices in ``{-r,...,-1,1,...,r}`` where ``r`` is `.rank`
        %
        % Returns:
        %   `.FreeGroupWord`: The reduced word
            w = replab.FreeGroupWord.make(self, letters);
        end

        function r = rank(self)
        % Returns the rank of this free group
        %
        % Returns:
        %   integer: Number of generators of this free group
            r = length(self.generatorNames);
        end

        function p = generator(self, i)
        % Returns the i-th group generator
        %
        % Args:
        %   i (integer): Generator index
        %
        % Returns:
        %   `+replab.FreeGroupWord`: i-th group generator

            p = self.word([i]);
        end

        function p = generatorInverse(self, i)
        % Returns the inverse of the i-th group generator
        %
        % Args:
        %   i (integer): Generator index
        %
        % Returns:
        %   element: Inverse of the i-th group generator
            p = self.word([-i]);
        end

    end

end
