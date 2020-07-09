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

    methods (Static)

        function id = newId
            persistent lastId
            if isempty(lastId)
                lastId = 0;
            end
            id = lastId + 1;
            lastId = id;
        end

        function [group, varargout] = parse(str)
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

    methods % superclass implementations

        function z = compose(self, x, y)
            z = x*y;
        end

        function z = inverse(self, x)
            z = inv(x);
        end

    end

    methods

        function w = word(self, arg)
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

        function sub = mrdivide(self, rel)
            if ~iscell(rel)
                rel = {rel};
            end
            assert(all(cellfun(@(r) r.group.groupId == self.groupId, rel)));
            relators = horzcat(self.relators, rel);
            relatorLetters = cellfun(@(r) r.letters, relators, 'uniform', 0);
            sub = replab.FiniteFPGroup(self.names, relatorLetters, replab.FPGroup.newId);
        end

        function nR = nRelators(self)
            nR = length(self.relators);
        end

        function r = relator(self, i)
            r = self.relators{i};
        end

        function s = stringPresentation(self)
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

end
