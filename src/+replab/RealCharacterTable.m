classdef RealCharacterTable < replab.CharacterTable
% Describes the real character table of a group
%
% The character values are stored as elements of the cyclotomic field, using the `.cyclotomic` class which requires
% external libraries and a Java Virtual Machine available.
%
% Instances of `.RealCharacterTable` are immutable.
%
% For additional help, see `.CharacterTable`.

    methods

        function self = RealCharacterTable(group, classes, values, varargin)
        % Constructs a character table
        %
        % Args:
        %   group (`.FiniteGroup`): Group represented by character table
        %   classes (`.ConjugacyClasses`): Conjugacy classes of `.group`
        %   values (`.cyclotomic` (nIrreps, nClasses)): Character values
        %
        % Keyword Args:
        %   irreps (cell(1,\*) of ``[]`` or `+replab.RepByImages`, optional): Explicit matrix representations (can contain empty values)
        %   classNames (cell(1,\*) of charstring, optional): Names of conjugacy classes
        %   irrepNames (cell(1,\*) of charstring, optional): Names of irreducible representations
        %   kronecker (integer(\*,\*,\*), optional): Kronecker coefficients
            self@replab.CharacterTable(group, 'R', classes, values, varargin{:});
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
            l = replab.laws.RealCharacterTableLaws(self);
        end

    end

    methods

        function res = imap(self, f)
            classes1 = self.classes.imap(f);
            group1 = f.target;
            values1 = self.values;
            irreps1 = cell(1, self.nIrreps);
            for i = 1:self.nIrreps
                if ~isempty(self.irreps{i})
                    irreps1{i} = self.irreps{i}.imap(f);
                end
            end
            res = replab.RealCharacterTable(group1, classes1, values1, 'irreps', irreps1, 'classNames', self.classNames, 'irrepNames', self.irrepNames);
            if self.inCache('kronecker')
                res.cache('kronecker', self.kronecker, 'error');
            end
        end

    end

end
