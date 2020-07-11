classdef FreeGroup < replab.Group
% Describes a free group

    properties (SetAccess = protected)
        id % (integer): Unique group id
        names % (cell(1,\*) of charstring): Names of the generators
    end

    methods

        function self = FreeGroup(names)
            self.identity = replab.FreeGroupWord(
            self.id = replab.globals.nextUniqueId;
            self.names = names;
        end

    end

    methods (Static) % FreeGroup construction

        function [F, varargout] = of(varargin)
            F = replab.FreeGroup(varargin);
            for i = 1:F.nGenerators
                varargout{i} = F.generator(i);
            end
        end

    end

    methods % Implementations

        function res = eq(self, rhs)
            res = (self.id == rhs.id);
        end

        function res = ne(self, rhs)
            res = self.id ~= rhs.id;
        end

        function s = headerStr(self)
            s = ['Free group < ' strjoin(self.names, ', ') ' >'];
        end

        function x = sample(self)
            l = 10;
            x = randi([-self.n self.n], 1, l);
            x = self.reduce(x(x ~= 0));
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

    methods

        function w = word(self, arg)
        % Constructs a word either from a string or an integer vector of letters
        %
        % Example:
        %   >>> [F x y] = replab.FreeGroup.of('x', 'y');
        %   >>> F.word('x (x y)^2 / x')
        %       x^2 y x y x^-1
            if ischar(arg)
                w = replab.FreeGroupWord.parse(self, arg);
            elseif isa(arg, 'double')
                w = replab.FreeGroupWord.make(self, arg);
            end
        end

        function r = rank(self)
        % Returns the rank of this free group
        %
        % Returns:
        %   integer: Number of generators of this free group
            r = length(self.names);
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
