classdef FiniteSetLaws < replab.laws.DomainLaws

    properties (SetAccess = protected)
        T % (`+replab.FiniteGroupType`): Finite group type
    end

    methods

        function self = FiniteSetLaws(S)
            self@replab.laws.DomainLaws(S);
            self.T = S.type;
        end

        function l = isSlow(self)
        % Returns whether some computations on this object are particularly slow
        %
        % Returns true when `.S` has more than 500 elements and is
        %
        % * a conjugacy class, as the elements of a conjugacy class cannot be enumerated in order,
        % * a coset, as some computations require computing all elements.
            l = isa(self, 'replab.DoubleCoset') || isa(self, 'replab.Coset') || isa(self, 'replab.ConjugacyClass');
            l = l && self.nElements > 500;
        end

    end

    methods % Laws

        % contains

        function law_set_contains_random_element_S(self, s)
            assertTrue(self.S.contains(s));
        end

        % elements

        function law_elements_(self)
            if self.S.nElements < 100
                E = self.S.elements;
                assertTrue(length(E) == self.S.nElements);
                for i = 1:length(E)
                    assertTrue(self.S.contains(E{i}));
                end
            end
        end

        % elementsSequence

        function laws = laws_elementsSequence(self)
            if self.isSlow
                laws = replab.Laws.empty;
            else
                laws = self.S.elementsSequence.laws;
            end
        end

        % nElements

        function law_nElements_(self)
            assert(self.S.nElements >= 0);
        end

        % representative

        function law_representative_is_minimal_S(self, s)
            assertTrue(self.S.type.compare(self.S.representative, s) <= 0);
        end

        function law_representative_is_inside_(self)
            assertTrue(self.S.contains(self.S.representative));
        end

        % setProduct

        function laws = laws_setProduct(self)
            laws = self.S.setProduct.laws(self.S);
        end

% $$$         function laws = laws_elementsSequence(self)
% $$$             if isa(self.S, 'replab.ConjugacyClass') && self.S.nElements > 500
% $$$                 return
% $$$             end
% $$$             laws = self.S.
% $$$             E = self.S.elementsSequence;
% $$$             assertTrue(E.nElements == self.S.nElements);
% $$$             % test one representative
% $$$             s = E.sample;
% $$$             assertTrue(S.contains(s));
% $$$             ind = E.find(s);
% $$$             s1 = E.at(ind);
% $$$             self.S.assertEqv(s, s1);
% $$$         end




    end

end
