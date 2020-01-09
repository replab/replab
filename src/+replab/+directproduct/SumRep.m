classdef SumRep < replab.Rep
% Representation of a direct product using a direct sum of factor representations
    
    properties (SetAccess = protected)
        factorReps % row cell array of replab.Rep: Representations of factors
    end
    
    methods
        
        function self = SumRep(group, factorReps)
        % Constructs a representation of a direct product
        %
        % See `+replab.+directproduct.OfCompactGroups.directSumRep`
            assert(isa(group, 'replab.directproduct.OfCompactGroups'));
            n = group.nFactors;
            assert(n >= 1, 'Direct sum cannot be empty');
            assert(length(factorReps) == n);
            self.field = factorReps{1}.field;
            self.factorReps = factorReps;
            d = 0;
            for i = 1:n
                assert(factorReps{i}.group == group.factor(i));
                assert(isequal(factorReps{i}.field, self.field));
                d = d + factorReps{i}.dimension;
            end
            self.dimension = d;
            factorRepsAreUnitary = cellfun(@(x) x.isUnitary, factorReps, 'uniform', 0);
            self.isUnitary = replab.trileanAnd(factorRepsAreUnitary{:});
            self.group = group;
        end
        
        function rho = image(self, g)
            n = length(g);
            rhos = arrayfun(@(i) self.factorReps{i}.image(g{i}), 1:n, 'uniform', 0);
            rho = blkdiag(rhos{:});
        end
        
    end
    
end
