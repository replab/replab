classdef SetProductLaws < replab.laws.DomainLaws

    properties (SetAccess = protected)
        finiteSet % (`+replab.FiniteSet` or ``[]``): Represented finite set, optional
    end

    methods

        function self = SetProductLaws(setProduct, finiteSet)
            if nargin < 2
                finiteSet = [];
            end
            self@replab.laws.DomainLaws(setProduct);
            self.finiteSet = finiteSet;
        end

        function law_identityFirst_(self)
            if self.S.identityFirst
                sets = self.S.sets;
                monoid = self.S.monoid;
                for i = 1:length(sets)
                    self.assert(monoid.isIdentity(sets{i}{1}));
                end
            end
        end

        function law_check_size_(self)
            if ~isempty(self.finiteSet)
                size1 = self.finiteSet.nElements;
                size2 = replab.util.multiplyIntegers(cellfun(@length, self.S.sets));
                self.assert(size1 == size2);
            end
        end

        function law_sample_contained_S(self, s)
            if ~isempty(self.finiteSet)
                self.assert(self.finiteSet.contains(s));
            end
        end

    end

end
