classdef RightCosets < replab.RightCosets

    methods

        function self = RightCosets(group, subgroup)
            self.group = group;
            self.subgroup = subgroup;
        end

    end

    methods % Implementations

        function t = cosetRepresentative(self, g)
            t = replab.bsgs.Cosets.rightRepresentative(self.subgroup.lexChain, g);
        end

        function T = transversal(self)
        % Returns all the canonical representatives of cosets
        %
        % Returns:
        %   cell(1, \*) of `.group` elements: Transversal
            M = replab.bsgs.Cosets.rightTransversalMatrix(self.group.lexChain, self.subgroup.lexChain);
            T = arrayfun(@(i) M(:,i)', 1:double(self.nElements), 'uniform', 0);
        end

    end

end
