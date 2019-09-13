classdef SumRep < replab.Rep
    
    properties (SetAccess = protected)
        reps % row cell array of replab.Rep: Representations of factors
    end
    
    methods
        
        function self = SumRep(group, reps)
            assert(isa(group, 'replab.directproduct.OfCompactGroups'));
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
            self.isUnitary = replab.domain.Trilean.and(cellfun(@(x) x.isUnitary, reps, 'uniform', 0));
            self.group = group;
        end
        
        function rho = image(self, g)
            n = length(g);
            rhos = arrayfun(@(i) self.reps{i}.image(g{i}), 1:n, 'uniform', 0);
            rho = blkdiag(rhos{:});
        end
        
    end
    
end
