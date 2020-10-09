classdef FactorizationTrivial < replab.mrp.Factorization
% Computes the factorization of a permutation group elements in its generators

    methods

        function self = FactorizationTrivial(group, useInverses)
            if nargin < 2 || isempty(useInverses)
                useInverses = false;
            end
            self.group = group;
            self.generators = cell(1, 0);
            self.useInverses = useInverses;
        end

        function letters = preimageElement(self, g)
            letters = zeros(1, 0);
        end

    end

end
