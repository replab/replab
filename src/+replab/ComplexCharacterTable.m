classdef ComplexCharacterTable < replab.CharacterTable
% Describes the standard character table of a group
%
% Example:
% %   >>> G = replab.PermutationGroup.dihedral(3); % doctest: +cyclotomic
% %   >>> G.characterTable
% %       Class  1a   3a   2a
% %        Size   1    2    3
% %
% %         X.1  1    1    1
% %         X.2  1    1   -1
% %         X.3  2   -1    0
%
% The character values are stored as elements of the cyclotomic field, using the `.cyclotomic` class which requires
% external libraries and a Java Virtual Machine available.
%
% Instances of `.ComplexCharacterTable` are immutable.
%
% See `.CharacterTable` for additional details.

    methods

        function self = ComplexCharacterTable(group, classes, values, varargin)
        % Constructs a character table
        %
        % Args:
        %   group (`.FiniteGroup`): Group represented by character table
        %   classes (`.ConjugacyClasses`): Conjugacy classes of `.group`
        %   values (`.cyclotomic` (nClasses, nClasses)): Character values
        %
        % Keyword Args:
        %   irreps (cell(1,\*) of ``[]`` or `+replab.RepByImages`, optional): Explicit matrix representations (can contain empty values)
        %   classNames (cell(1,\*) of charstring, optional): Names of conjugacy classes
        %   irrepNames (cell(1,\*) of charstring, optional): Names of irreducible representations
        %   kronecker (integer(\*,\*,\*), optional): Kronecker coefficients
            self@replab.CharacterTable(group, 'C', classes, values, varargin{:});
            nIrreps = size(values, 1);
            nClasses = size(values, 2);
            assert(nIrreps == nClasses, 'Complex character tables should have the same number of characters and conjugacy classes.');
        end

        function ind = conjugateCharacterIndices(self)
        % For each character, returns the index of the character conjugate to it
        %
        % Returns:
        %   integer(1,\*): Index of the conjugate character

        end

% $$$         function [R, C, H] = types(self)
% $$$         % Returns the indices of the real/complex/quaternion-type irreducible representations
% $$$         %
% $$$         % Returns
% $$$         % -------
% $$$         %   R: integer(1,\*)
% $$$         %     Indices of the real-type representations
% $$$         %   C: integer(1/2,\*)
% $$$         %     Indices of the complex-type representations, of conjugate pairs if `.overC` is true
% $$$         %   H: integer(1,\*)
% $$$         %     Indices of the quaternion-type representations
% $$$             F = cellfun(@(C) C.frobeniusSchurIndicator, self.characters);
% $$$             R = find(F == 1);
% $$$             R = R(:)';
% $$$             if self.overC
% $$$                 H = find(F == -1);
% $$$                 H = H(:)';
% $$$             else
% $$$                 H = find(F == -2);
% $$$                 H = H(:)';
% $$$             end
% $$$             inds = find(F == 0);
% $$$             if self.overR
% $$$                 C = inds;
% $$$             else
% $$$                 C = zeros(2, 0);
% $$$                 while ~isempty(inds)
% $$$                     ind = inds(1);
% $$$                     inds = inds(2:end);
% $$$                     v = self.values(ind,:);
% $$$                     indC = [];
% $$$                     for i = inds
% $$$                         if self.values(i,:) == conj(v)
% $$$                             indC = i;
% $$$                             break
% $$$                         end
% $$$                     end
% $$$                     if isempty(indC)
% $$$                         error('Inconsistent character table');
% $$$                     end
% $$$                     C = [C [ind; indC]];
% $$$                     inds = setdiff(inds, indC);
% $$$                 end
% $$$             end
% $$$         end

    end

    methods % Implementations

        % Obj

        function l = laws(self)
            l = replab.laws.ComplexCharacterTableLaws(self);
        end

    end

    methods (Static)

        function C = fromRealCharacterTable(R)
        % Creates a complex character table from a real character table
        %
        % This works only if its real irreducible representations are absolutely irreducible;
        % though that fact is not checked by this method.
        %
        % Args:
        %   R (`.RealCharacterTable`): Real character table
        %
        % Returns:
        %   `.ComplexCharacterTable`: The corresponding complex character table
            irreps = cell(1, R.nIrreps);
            for i = 1:R.nIrreps
                if R.hasIrrep(i)
                    images = cell(1, R.group.nGenerators);
                    for j = 1:R.group.nGenerators
                        images{j} = R.irrep(i).image(R.group.generator(j), 'exact');
                    end
                    irreps{i} = R.group.repByImages('C', R.irrep(i).dimension, 'preimages', R.group.generators, 'images', images);
                else
                    irreps{i} = [];
                end
            end
            C = replab.ComplexCharacterTable(R.group, R.classes, R.values, 'irreps', irreps, 'classNames', R.classNames, 'irrepNames', R.irrepNames);
        end

        function C = fromIrreps(group, irreps)
        % Creates a complex character table from a complete list of irreducible exact representations
        %
        % Args:
        %   group (`.FiniteGroup`): Finite group
        %   irreps (cell(1,\*) of `.Rep`): Irreducible representations
        %
        % Returns:
        %   `.ComplexCharacterTable`: The corresponding complex character table
            assert(all(cellfun(@(r) isa(r, 'replab.Rep') && r.overC, irreps)));
            nIrreps = length(irreps);
            cc = group.conjugacyClasses;
            nClasses = cc.nClasses;
            values = replab.cyclotomic.zeros(nIrreps, nClasses);
            for i = 1:nIrreps
                for j = 1:nClasses
                    r = cc.classes{j}.representative;
                    values(i,j) = trace(irreps{i}.image(r, 'exact'));
                end
            end
            C = replab.ComplexCharacterTable(group, cc, values, 'irreps', irreps);
        end

    end

    methods

        function res = imap(self, f)
        % Maps the character table under an isomorphism
        %
        % Example:
%         %   >>> D6a = replab.PermutationGroup.of([3 2 1], [2 3 1]); % doctest: +cyclotomic
%         %   >>> D6b = replab.PermutationGroup.of([1 4 3 2], [1 3 4 2]);
%         %   >>> f = D6a.isomorphismByImages(D6b, 'preimages', D6a.generators, 'images', D6b.generators);
%         %   >>> Ca = D6a.characterTable;
%         %   >>> Cb = Ca.imap(f);
%         %   >>> Cb.laws.checkSilent;
        %
        % Args:
        %   f (`.FiniteIsomorphism`): Isomorphism with ``self.group.isSubgroupOf(f.source)``
        %
        % Returns:
        %   `.ComplexCharacterTable`: The character table of the subgroup in the image of the isomorphism
            if self.group.order < f.source.order
                f = f.restrictedSource(self.group);
            end
            classes1 = self.classes.imap(f);
            group1 = f.target;
            values1 = self.values;
            irreps1 = cell(1, self.nIrreps);
            for i = 1:self.nIrreps
                if ~isempty(self.irreps{i})
                    irreps1{i} = self.irreps{i}.imap(f);
                end
            end
            res = replab.ComplexCharacterTable(group1, classes1, values1, 'irreps', irreps1, 'classNames', self.classNames, 'irrepNames', self.irrepNames);
            if self.inCache('kronecker')
                res.cache('kronecker', self.kronecker, 'error');
            end
        end

    end

end
