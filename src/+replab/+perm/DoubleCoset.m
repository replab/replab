classdef DoubleCoset < replab.DoubleCoset

    methods

        function self = DoubleCoset(representative, leftSubgroup, rightSubgroup, group)
            assert(group.hasSameTypeAs(leftSubgroup));
            assert(group.hasSameTypeAs(rightSubgroup));
            self.type = group.type;
            self.representative_ = representative;
            self.leftSubgroup = leftSubgroup;
            self.rightSubgroup = rightSubgroup;
            self.group = group;
        end

    end

    methods % Implementations

        function b = contains(self, el)
            if ~self.group.contains(el)
                b = false;
                return
            end
            dc = self.leftSubgroup.doubleCoset(el, self.rightSubgroup, 'group', self.group);
            b = self.type.eqv(self.representative, dc.representative);
        end

        function E = elementsSequence(self)
            S = replab.perm.Set(self.leftSubgroup.domainSize);
            rightMat = self.rightSubgroup.lexChain.allElements;
            g = self.representative';
            S.insert(g(rightMat));
            toCheck = 1;
            while ~isempty(toCheck)
                cur = S.at(toCheck(end))';
                toCheck = toCheck(1:end-1);
                for j = 1:self.leftSubgroup.nGenerators
                    h = self.leftSubgroup.generator(j);
                    newEl = h(cur);
                    if S.find(newEl') == 0
                        newEl = newEl';
                        sz = S.nElements;
                        inds = S.insert(newEl(rightMat));
                        assert(all(inds > sz));
                        toCheck = [toCheck sz+1];
                    end
                end
            end
            S.sort;
            E = replab.perm.Sequence(S.matrix);
        end

        function s = nElements(self)
            s = self.cached('nElements', @() self.computeNElements);
        end

    end

    methods (Access = protected)

        function s = computeNElements(self)
        % From Wikipedia: |left x right| = |left| |right| / |right \intersection x^-1 left x|
            leftConj = self.leftSubgroup.leftConjugateGroup(self.group.inverse(self.representative));
            inter = self.rightSubgroup.intersection(leftConj);
            s = self.leftSubgroup.order * self.rightSubgroup.order / inter.order;
        end

    end

end
