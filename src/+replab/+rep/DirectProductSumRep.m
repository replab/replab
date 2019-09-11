classdef DirectProductSumRep < replab.Rep
    
    properties (SetAccess = protected)
        reps;
    end
    
    methods
        
        function self = DirectProductSumRep(group, reps)
            assert(isa(group, 'replab.DirectProductGroup'));
            n = group.nFactors;
            assert(n >= 1, 'Direct sum cannot be empty');
            assert(length(reps) == n);
            self.field = reps{1}.field;
            self.reps = reps;
            d = 0;
            for i = 1:n
                assert(reps{i}.group == group.factor(i));
                assert(isequal(reps{i}.field, self.field));
                d = d + reps{i}.dimension;
            end
            self.dimension = d;
            self.group = group;
        end
        
        function rho = image(self, g)
            n = length(g);
            rhos = arrayfun(@(i) self.reps{i}.image(g{i}), 1:n, 'uniform', 0);
            rho = blkdiag(rhos{:});
        end
        
    end
    
end
