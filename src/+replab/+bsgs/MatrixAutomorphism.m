classdef MatrixAutomorphism < replab.bsgs.Backtrack
% Computes the unordered partition stabilizer of a group

    properties
        matrix % (double(\*,\*)): Matrix to stabilize
    end

    methods

        function self = MatrixAutomorphism(group, matrix, knownSubgroup, debug)
            if nargin < 4 || isempty(debug)
                debug = false;
            end
            if nargin < 3 || isempty(knownSubgroup)
                knownSubgroup = [];
            end
            n = size(matrix, 1);
            assert(n == size(matrix, 2));
            assert(n == group.domainSize);
            base = 1:n;
            self@replab.bsgs.Backtrack(group, base, knownSubgroup, knownSubgroup, debug);
            self.matrix = matrix;
        end


        function ok = test(self, l, prev, ul)
        % Verifies that blocks are mapped to blocks consistently
            beta = self.base(l);
            b = prev(ul(beta));
            ok = false;
            partial = self.base(1:l);
            images = prev(ul(partial));
            if ~isequal(self.matrix(partial, partial), self.matrix(images, images))
                return
            end
            if ~isequal(sort(self.matrix(:,beta)), sort(self.matrix(:,b)))
                return
            end
            if ~isequal(sort(self.matrix(beta,:)), sort(self.matrix(b,:)))
                return
            end
            ok = true;
        end

        function ok = prop(self, g)
        % Verifies the image of every block is a block of the same partition
            ok = isequal(self.matrix, self.matrix(g, g));
        end

    end

end
