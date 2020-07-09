classdef FPGroup < replab.GroupWithGenerators
% Describes a finitely presented group
%
% The group is defined on ``n`` generators ``s_1`` to ``s_n``. By convention, we identify the integers ``1`` to ``n``
% with the generators, while we identify the integers ``-1`` to ``-n`` to the inverses of the generators.
%
% A relator is a product of generators and their inverses which is equal to the identity. We write
% it as a row integer vector using the numbering scheme above.

    properties
        groupId % (integer): Unique group id
        n % (integer): Number of generators
        names % (cell(1,\*) of charstring): Names of the generators
        relators % (cell(1,\*) of integer(1,\*)): Relators
    end

    methods (Access = protected)

        function self = FPGroup(names, relatorLetters, id)
            self.identity = replab.Word(self, []);
            n = length(names);
            self.groupId = id;
            self.n = n;
            self.generators = arrayfun(@(i) replab.Word(self, i), 1:n, 'uniform', 0);
            self.relators = cellfun(@(l) replab.Word(self, l), relatorLetters, 'uniform', 0);
            if nargin < 3 || isequal(names, [])
                self.names = arrayfun(@(i) ['s' num2str(i)], 1:n, 'uniform', 0);
            else
                self.names = names;
            end
        end

    end


    methods (Static) % Builders

        function [group, varargout] = parse(str)
        % Creates a finitely presented group from a description string
        %
        % Example:
        %   >>> [F a x] = replab.FPGroup.parse('<a, x>'); % free group on two elements
        %   >>> isa(F, 'replab.FreeGroup')
        %       1
        %   >>> [G x] = replab.FiniteFPGroup.parse('< x | x^3 = 1 >');
        %   >>> isa(G, 'replab.FiniteFPGroup')
        %       1
        %
        % Args:
        %   str (charstring): Single-line description string
        %
        % Returns:
        %   `+replab.FPGroup`: The parsed finitely presented group
            if nargin < 2
                permutations = [];
            end
            P = replab.fp.Parser;
            [ok, names, relators] = P.parsePresentation(str);
            assert(ok, 'Error in given presentation string');
            id = replab.FPGroup.newId;
            if isempty(relators)
                group = replab.FreeGroup(names, id);
            else
                group = replab.FiniteFPGroup(names, relators, id);
            end
            if nargout > 1
                for i = 1:length(names)
                    varargout{i} = group.generator(i);
                end
            end
        end

    end

    methods % Superclass implementations

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
                P = replab.fp.Parser;
                [ok, letters] = P.parseWord(arg, self.names);
                assert(ok, 'This string does not define a word');
                w = replab.Word(self, letters);
            elseif isa(arg, 'double')
                assert(all(arg ~= 0), 'Words have letters between -n ... -1 and 1 ... n');
                assert(arg == round(arg));
                assert(all(arg >= -self.n) && all(arg <= self.n), 'Generator index too big');
                w = replab.Word(self, letters);
            end
        end


        function nR = nRelators(self)
        % Returns the number of relators describing this group
        %
        % Returns:
        %   integer: Number of relators
            nR = length(self.relators);
        end

        function r = relator(self, i)
        % Returns the ``i``-th relator in the presentation of this group
        %
        % Args:
        %   i (integer): Index of relator
        % Returns:
        %   `+replab.Word`: ``i``-th relator
            r = self.relators{i};
        end

        function s = stringPresentation(self)
        % Returns a readable text presentation of this group
        %
        % It should be able to be parsed using `+replab.FPGroup.parse`.
        %
        % Returns:
        %   charstring: Presentation as string
            gens = strjoin(cellfun(@(w) w.word2str, self.generators, 'uniform', 0), ', ');
            if self.nRelators == 0
                s = ['< ' gens ' >'];
            else
                rels = strjoin(cellfun(@(r) r.word2str, self.relators, 'uniform', 0), ' = ');
                s = ['< ' gens ' | ' rels ' = 1 >'];
            end
        end

        function s = headerStr(self)
            s = self.stringPresentation;
        end

        function l = imagesDefineMorphism(self, target, generatorImages)
        % Checks whether the given images satisfy the relations of the presentation of this group
        %
        % If it returns true, it means those images describe a valid homomorphism from this `.FPGroup`
        % to the given target group.
        %
        % Args:
        %   target (`+replab.Group`): Target group
        %   generatorImages (cell(1,\*) of elements of ``target``): Images of the generators of this group
        %
        % Returns:
        %   logical: True if the images verify the presentation
            nR = length(self.relators);
            for i = 1:nR
                r = self.relators{i};
                if ~target.isIdentity(r.computeImage(target, generatorImages))
                    l = false;
                    return
                end
            end
            l = true;
        end

        function m = morphismByImages(self, target, generatorImages)
            assert(self.imagesDefineMorphism(target, generatorImages));
            m = replab.Morphism.lambda(self, target, @(w) w.computeImage(target, generatorImages));
        end

    end

    methods (Static)

        function id = newId
        % Static function that returns unique integers
        %
        % This is used to verify to which group words belong.
            persistent lastId
            if isempty(lastId)
                lastId = 0;
            end
            id = lastId + 1;
            lastId = id;
        end

    end
end
