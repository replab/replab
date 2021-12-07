classdef LeftCosets < replab.LeftCosets

    methods

        function self = LeftCosets(group, subgroup)
            self.group = group;
            self.subgroup = subgroup;
        end

    end

    methods % Implementations

        function t = cosetRepresentative(self, g)
            t = replab.bsgs.Cosets.leftRepresentative(self.subgroup.lexChain, g);
        end

        function mu = leftAction(self)
            nG = self.group.nGenerators;
            ds = self.group.domainSize;
            T = self.transversalAsMatrix;
            S = replab.perm.Set(ds);
            S.insert(T);
            n = size(T, 2);
            images = cell(1, nG);
            for i = 1:nG
                g = self.group.generator(i);
                img = zeros(1, n);
                for j = 1:n
                    gt = replab.bsgs.Cosets.leftRepresentative(self.subgroup.lexChain, g(T(:,j)'));
                    loc = S.find(gt');
                    assert(length(loc) == 1);
                    img(j) = loc;
                end
                images{i} = img;
            end
            Sn = replab.S(n);
            mu = self.group.morphismByImages(Sn, 'images', images);
        end

        function T = transversal(self)
            M = self.transversalAsMatrix;
            T = arrayfun(@(i) M(:,i)', 1:double(self.nElements), 'uniform', 0);
        end

    end

    methods % Permutation group specific methods

        function M = transversalAsMatrix(self)
            M = self.cached('transversalAsMatrix', @() self.computeTransversalAsMatrix);
        end

    end

    methods (Access = protected)

        function M = computeTransversalAsMatrix(self)
            M = replab.bsgs.Cosets.leftTransversalAsMatrix(self.group.lexChain, self.subgroup.lexChain);
        end

    end

end
