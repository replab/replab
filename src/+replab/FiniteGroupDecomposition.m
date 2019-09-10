classdef FiniteGroupDecomposition < replab.Str
% Describes the decomposition of a "group" into a product of coset representatives
%
% Describes the set U = { u }, where u = u_1 u_2 ... u_n and
% u_i is in transversals{i}
%
% For each i, we have transversals{i}{1} the group identity.
%
% The set U possibly has repetitions of elements, and composition
% of elements is done according to the binary operation of "group".
%
% A construction comes from the decomposition of a group into
% a chain of subgroups
% G = G_1 >= G_2 >= ... G_{n+1} = trivial group
% where transversals{i} are (left) coset representatives of G_i / G_{i+1}
    properties (SetAccess = protected)
        group;
        transversals;
    end
    methods
        function self = FiniteGroupDecomposition(group, transversals)
            self.group = group;
            self.transversals = transversals;
        end
    end
    methods (Static)
        function D = trivial(group, elements)
            if ~group.isIdentity(elements{1})
                idIndex = find(cellfun(@(g) group.isIdentity(g), elements));
                assert(length(idIndex) == 1, 'Elements should have a single copy of the identity');
                elements = elements([idIndex setdiff(1:length(elements, idIndex))]);
            end
            D = FiniteGroupDecomposition(group, {elements});
        end
    end
end
