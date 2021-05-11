classdef FactorizationTrivial < replab.mrp.Factorization
% Computes the factorization of a permutation group elements in its generators

    methods

        function self = FactorizationTrivial(group, useInverses)
            assert(group.isTrivial);
            if nargin < 2 || isempty(useInverses)
                useInverses = false;
            end
            self.group = group;
            self.generators = cell(1, 0);
            self.useInverses = useInverses;
        end

    end

    methods % Implementations

        function letters = factorize(self, g)
            letters = zeros(1, 0);
        end

        function [l, r] = factorizeRepresentativeOfLeftCoset(self, leftCoset)
            l = zeros(1, 0);
            r = self.group.identity;
        end

    end

end
