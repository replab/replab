classdef DirectProductTensorRep < replab.Rep
    
    properties (SetAccess = protected)
        reps;
    end
    
    methods
        
        function self = DirectProductTensorRep(group, reps)
            assert(isa(group, 'replab.DirectProductGroup'));
            n = group.nFactors;
            assert(n >= 1, 'Direct product cannot be empty');
            assert(length(reps) == n);
            self.field = reps{1}.field;
            self.reps = reps;
            d = 1;
            for i = 1:n
                assert(isequal(reps{i}.group, group.factor(i)));
                assert(isequal(reps{i}.field, self.field));
                d = d * reps{i}.dimension;
            end
            self.dimension = d;
            self.group = group;
        end
        
        function rho = image(self, g)
            rho = self.reps{1}.image(g{1});
            for i = 2:self.group.nFactors
                rho = kron(rho, self.reps{i}.image(g{i}));
            end
        end
        
    end
    
end
