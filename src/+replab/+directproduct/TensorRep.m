classdef TensorRep < replab.Rep
% Representation of a direct product using a tensor product of factor representations
    
    properties (SetAccess = protected)
        factorReps % row cell array of replab.Rep: Representations of factors
    end
    
    methods
        
        function self = TensorRep(group, factorReps, varargin)
        % Constructs a representation of a direct product
        %
        % See `+replab.+directproduct.OfCompactGroups.tensorRep`
            
            assert(isa(group, 'replab.directproduct.OfCompactGroups'));
            n = group.nFactors;
            assert(n >= 1, 'Direct product cannot be empty');
            assert(length(factorReps) == n);
            self.field = factorReps{1}.field;
            self.factorReps = factorReps;
            d = 1;
            for i = 1:n
                assert(factorReps{i}.group == group.factor(i));
                assert(isequal(factorReps{i}.field, self.field));
                d = d * factorReps{i}.dimension;
            end
            self.dimension = d;
            factorRepsAreUnitary = cellfun(@(x) x.isUnitary, factorReps, 'uniform', 0);            
            self.isUnitary = replab.trileanAnd(factorRepsAreUnitary{:});
            self.group = group;
        end
        
        function rho = image(self, g)
            rho = self.factorReps{1}.image(g{1});
            for i = 2:self.group.nFactors
                rho = kron(rho, self.factorReps{i}.image(g{i}));
            end
        end
        
    end
    
end
