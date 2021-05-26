classdef AtlasEntry < replab.Obj
% Identifies a user-defined group as a standard group present in an atlas

    properties (SetAccess = protected)
        name % (charstring): Group name
        group % (`.AbstractGroup`): Group
    end

    methods

        function self = AtlasEntry(name, group, characterTable)
            self.name = name;
            self.group = group;
        end

        function res = canMatch(self, group)
        % Returns whether this atlas entry can match the given group
            res = false;
            gAE = cellfun(@(d) d.abelianInvariants, group.derivedSeries, 'uniform', 0);
            myAE = cellfun(@(d) d.abelianInvariants, self.group.permutationGroup.derivedSeries, 'uniform', 0);
            if ~isequal(gAE, myAE)
                return
            end
            res = true;
            % TODO: tests based on conjugacy classes
        end

        function m = isomorphism(self, group)
            m = self.group.findIsomorphism(group);
        end

        function r = match(self, group)
            m = self.isomorphism(group);
            if isempty(m)
                r = [];
                return
            end
            r = replab.AtlasResult(group, self, m);
        end

        function laws = laws(self)
            laws = replab.laws.AtlasEntryLaws(self);
        end

    end

    methods (Static) % Entries for common finite groups

        function A = trivial
        % Constructs the atlas entry corresponding to the trivial group
            name = 'Trivial group';
            generators = cell(1, 0);
            relators = cell(1, 0);
            % Presentation is empty
            ag = replab.AbstractGroup(generators, relators, 'permutationGenerators', {[]});
            classes = ag.conjugacyClasses;
            characters = replab.cyclotomic.eye(1);
            classNames = {'id'};
            irrepNames = {'id'};
            irreps = {ag.trivialRep('C', 1)};
            ct = replab.ComplexCharacterTable(ag, 'C', classes, characters, 'irreps', irreps, 'classNames', classNames, 'irrepNames', irrepNames);
            A = replab.AtlasEntry(name, ag, ct);
        end

        function A = dihedral(n)
        % Constructs the atlas entry corresponding to the dihedral group of order 2*n
            assert(n > 2);
            name = sprintf('Dihedral group of order %d', 2*n);
            % Permutation realization
            X = [n:-1:1];
            A = [2:n 1];
            % Presentation from the groupprops wiki
            % < x, a | a^n = x^2 = 1, x a x^-1 = a^-1 >
            relators = {['a^' num2str(n)] 'x^2' 'x a x^-1 a'};
            ag = replab.AbstractGroup({'x', 'a'}, relators, 'permutationGenerators', {X, A}, 'order', vpi(2*n));
            ct = replab.ct.DihedralCharacterTable(n);
            assert(ct.group == ag.permutationGroup);
            A = replab.AtlasEntry(name, ag, ct.imap(ag.niceMorphism.inverse));
        end

        function A = cyclic(n)
        % Constructs the cyclic group of order n
            assert(n >= 2);
            name = sprintf('Cyclic group C(%d) of order %d', n, n);
            % Permutation realization
            X = [2:n 1];
            % standard presentation
            % < x | x^n = 1 >
            ag = replab.AbstractGroup({'x'}, {['x^' num2str(n)]}, 'permutationGenerators', {X}, 'order', vpi(n));
            ct = replab.ct.CyclicCharacterTable(n);
            assert(ct.group == ag.permutationGroup);
            A = replab.AtlasEntry(name, ag, ct.imap(ag.niceMorphism.inverse));
        end

        function A = kleinFourGroup
        % Constructs the atlas entry corresponding to the klein four-group
            name = sprintf('Klein four-group of order %d', 4);
            % Permutation realization
            X = [2,1,4,3];
            A = [3,4,1,2];
            % Presentation from the groupprops wiki
            % < x, a | a^2 = x^2 = 1, x a x^-1 = a^-1 >
            relators = {'a^2' 'x^2' 'x a x^-1 a'};
            ag = replab.AbstractGroup({'x' 'a'}, relators, 'permutationGenerators', {X, A}, 'order', vpi(4));
            classes = replab.ConjugacyClasses(ag, cellfun(@(g) ag.conjugacyClass(g), {'1' 'x' 'a' 'x a'}, 'uniform', 0));
            chars = replab.cyclotomic([1 1 1 1; 1 1 -1 -1; 1 -1 1 -1; 1 -1 -1 1]);
            classNames = {'1' 'x' 'a' 'xa'};
            irrepNames = {'trivial' 'ker x' 'ker a' 'ker xa'};
            irreps = {ag.repByImages('C', 1, 'images', {replab.cyclotomic(1) replab.cyclotomic(1)}) ...
                      ag.repByImages('C', 1, 'images', {replab.cyclotomic(1) replab.cyclotomic(-1)}) ...
                      ag.repByImages('C', 1, 'images', {replab.cyclotomic(-1) replab.cyclotomic(1)}) ...
                      ag.repByImages('C', 1, 'images', {replab.cyclotomic(-1) replab.cyclotomic(-1)})};
            ct = replab.ComplexCharacterTable(ag, 'C', classes, chars, 'classNames', classNames, 'irrepNames', irrepNames, 'irreps', irreps);
            A = replab.AtlasEntry(name, ag, ct);
        end

        function A = symmetric(n)
        % Constructs the atlas entry corresponding to the symmetric group of degree n
            assert(n > 2);
            name = sprintf('Symmetric group S(%d) of degree %d', n, n);
            % Permutation realization
            S = [2:n 1];
            T = [2 1 3:n];
            % this is the presentation from page 2100 of
            % https://www.ams.org/journals/tran/2003-355-05/S0002-9947-03-03040-X/S0002-9947-03-03040-X.pdf
            relators = {['s^' num2str(n)], 't^2', ['(s*t)^' num2str(n-1)]};
            for j = 2:floor(n/2)
                relators{1,end+1} = sprintf('(t^-1 s^-%d t s^%d)^2', j, j);
            end
            ag = replab.AbstractGroup({'s' 't'}, relators, 'permutationGenerators', {S, T}, 'order', replab.util.factorial(n));
            ct = replab.sym.SymmetricGroupCharacterTable(n);
            A = replab.AtlasEntry(name, ag, ct.imap(ag.niceMorphism.inverse));
        end

        function A = alternating(n)
        % Constructs the alternating group of degree n
            assert(n >= 4);
            name = sprintf('Alternating group A(%d) of degree %d', n, n);
            isEven = mod(n, 2) == 0;
            % Permutation realization
            T = [2 3 1 4:n];
            if isEven
                S = [2 1 4:n 3];
            else
                S = [1 2 4:n 3];
            end
            % this is the presentation from page 2100 of
            % https://www.ams.org/journals/tran/2003-355-05/S0002-9947-03-03040-X/S0002-9947-03-03040-X.pdf
            if isEven
                relators = {['s^' num2str(n-2)], 't^3', ['(s t)^' num2str(n-1)], '(t^-1 s^-1 t s)^2'};
            else
                relators = {['s^' num2str(n-2)], 't^3', ['(s t)^' num2str(n)]};
                for k = 1:floor((n-3)/2)
                    relators{1,end+1} = sprintf('(t s^-%d t s^%d)^2', k, k);
                end
            end
            ag = replab.AbstractGroup({'s' 't'}, relators, 'permutationGenerators', {S, T}, 'order', replab.util.factorial(n)/2);
            A = replab.AtlasEntry(name, ag, []);
        end

    end

    methods (Static) % Generic construction

% $$$         function A = parse(s)
% $$$         % Parses a JSON string corresponding to an AtlasEntry
% $$$         %
% $$$         % Args:
% $$$         %   s (charstring): JSON string
% $$$         %
% $$$         % Returns:
% $$$         %   `.AtlasEntry`: The parsed entry
% $$$             J = replab.util.parseJSON(s);
% $$$             if isfield(J, 'name')
% $$$                 name = J.name;
% $$$             else
% $$$                 name = [];
% $$$             end
% $$$             G = replab.atlas.parseGroup(J.group);
% $$$             C = replab.AtlasEntry.parseCharacterTable(J, G);
% $$$             A = replab.AtlasEntry(name, G, C);
% $$$         end

    end

end
