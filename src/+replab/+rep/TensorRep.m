classdef TensorRep < replab.Rep
% A tensor product of representations
%
% All factor representations must be defined on the same group
    
    properties
        factors % row cell array of replab.Rep: Factor representations
    end
    
    methods
        
        function self = TensorRep(factors)
        % Constructs a tensor representation from a cell array of representations
        %
        % All the subrepresentations should be defined on the same group, and on the same field.
        %
        % Args:
        %   blocks (row cell array of replab.Rep): Factor representations
            assert(length(factors) >= 1);
            d = 1;
            for i = 1:length(factors)
                assert(isa(factors{i}, 'replab.Rep'));
                d = d * factors{i}.dimension;
            end
            self.dimension = d;
            factorsAreUnitary = cellfun(@(x) x.isUnitary, factors, 'uniform', 0);
            self.isUnitary = replab.trileanAnd(factorsAreUnitary{:});
            for i = 2:length(factors)
                assert(factors{1}.group == factors{i}.group);
                assert(isequal(factors{1}.field, factors{i}.field));
            end
            self.factors = factors;
            self.group = factors{1}.group;
            self.field = factors{1}.field;
        end
        
        function n = nFactors(self)
            n = length(self.factors);
        end
        
        function factor = factor(self, i)
            factor = self.factors{i};
        end
        
        % Str
        
        function names = hiddenFields(self)
            names = hiddenFields@replab.Rep(self);
            names{1, end+1} = 'factors';
        end
        
        function [names values] = additionalFields(self)
            [names values] = additionalFields@replab.Rep(self);
            for i = 1:self.nFactors
                names{1, end+1} = sprintf('factor(%d)', i);
                values{1, end+1} = self.factor(i);
            end
        end
        
        % Rep
        
        function rho = image(self, g)
            rho = self.factors{1}.image(g);
            for i = 2:self.nFactors
                rho = kron(rho, self.factors{i}.image(g));
            end
        end
        
    end
end
