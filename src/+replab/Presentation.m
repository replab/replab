classdef Presentation < replab.Str
% Defines a presentation of a group
%
% The group is defined on ``n`` generators ``s_1`` to ``s_n``. By convention, we identify the integers ``1`` to ``n``
% with the generators, while we identify the integers ``-1`` to ``-n`` to the inverses of the generators.
%
% A relator is a product of generator and their inverses with is equal to the identity. We write
% it as a row integer vector using the numbering scheme above.

    properties
        n % (integer): Number of generators
        relators % (cell(1,\*) of integer(1,\*)): Relators
        names % (cell(1,\*) of char): Names of the generators
    end

    methods

        function self = Presentation(n, relators, names)
            self.n = n;
            self.relators = relators;
            if nargin < 3 || isequal(names, [])
                self.names = arrayfun(@(i) ['s' num2str(i)], 1:n, 'uniform', 0);
            else
                self.names = names;
            end
        end

        function s = word2str(self, word)
            s = '';
            sep = '';
            i = 1;
            while i <= length(word)
                e = sign(word(i)); % exponent
                l = abs(word(i));
                assert(e ~= 0);
                i = i + 1;
                while i <= length(word) && abs(word(i)) == l
                    e = e + sign(word(i));
                    i = i + 1;
                end
                if e ~= 0
                    s = [s sep self.names{l}];
                    sep = ' ';
                    if e ~= 1
                        s = [s '^' num2str(e)];
                    end
                end
            end
        end

        function show(self)
            fprintf('Presentation with generators: %s\n', strjoin(self.names, ','));
            rels = strjoin(cellfun(@(r) self.word2str(r), self.relators, 'uniform', 0), ' = ');
            rels = [rels ' = id'];
            fprintf('%s\n', rels);
        end

        function g = computeWord(self, group, word)
            g = group.identity;
            for i = 1:length(word)
                w = word(i);
                if w > 0
                    g = group.compose(g, group.generator(w));
                else
                    g = group.composeWithInverse(g, group.generator(-w));
                end
            end
        end


        function res = isVerifiedBy(self, group)
            nR = length(self.relators);
            for i = 1:nR
                val = self.computeWord(group, self.relators{i});
                if ~group.isIdentity(val)
                    res = false;
                    return
                end
            end
            res = true;
        end

    end

end
