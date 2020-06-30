classdef FiniteGroupDecomposition < replab.Obj
% Describes the decomposition of a finite group into a product of sets
%
% We assume the existence of sets T1, T2, ..., Tn such that every group elements has a unique decomposition
% g = t1 t2 ... tn, where ti is in Ti; thus, we have that length(T1)*length(T2)*...*length(Tn) = group order
%
% We require additionally that in each set Ti, Ti{1} is the group identity.
%
    properties (SetAccess = protected)
        group % (`+replab.FiniteGroup): Group decomposed
        T % (cell(1,\*) of cell(1,\*) of group elements): Stores the sets Ti as {T1 T2 ... Tn}
    end

    methods

        function self = FiniteGroupDecomposition(group, T)
            self.group = group;
            self.T = T;
        end

    end

    methods (Static)

        function D = trivial(group)
        % Constructs a group decomposition using a single set
        %
        % No efficiency gains are then expected.
        %
        % Args:
        %   group (replab.FiniteGroup): Finite group to decompose
            O = double(group.order);
            T = group.elements.toCell;
            idIndex = double(group.elements.find(group.identity));
            T = T([idIndex setdiff(1:O, idIndex)]);
            D = replab.FiniteGroupDecomposition(group, {T});
        end

    end

end
