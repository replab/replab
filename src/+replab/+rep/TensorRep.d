classdef TensorRep < replab.Rep
% A tensor product of representations
%
% All factor representations must be defined on the same group

    properties
        factors % (cell(1,\*) of `+replab.Rep`): Factor representations
    end

    methods

        function self = TensorRep(group, field, factors)
        % Constructs a tensor representation from a cell array of representations
        %
        % All the subrepresentations should be defined on the same group, and on the same field.
        %
        % Args:
        %   group (`+replab.CompactGroup`): Common group
        %   field ({'R', 'C'}): Real or complex field
        %   blocks (cell(1,\*) of `+replab.Rep`): Factor representations
            replab.rep.assertCompatibleFactors(group, field, factors);
            % own properties
            self.factors = factors;
            % replab.Rep immutable
            self.group = group;
            self.field = field;
            self.dimension = prod(cellfun(@(f) f.dimension, factors));
            % replab.Rep mutable
            factorsAreUnitary = cellfun(@(x) x.isUnitary, factors, 'uniform', 0);
            self.isUnitary = replab.trileanAnd(factorsAreUnitary{:});
        end

        function n = nFactors(self)
        % Returns the number of factors in the tensor product
        %
        % Returns:
        %   integer: Number of factors
            n = length(self.factors);
        end

        function f = factor(self, i)
        % Returns a factor in the tensor product
        %
        % Args:
        %   i (integer): Index of a factor
        %
        % Returns:
        %   `+replab.Rep`: Representation corresponding to the i-th factor
            f = self.factors{i};
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

        function s = headerStr(self)
            s = headerStr@replab.Rep(self); % logic in parent class
        end

        % Rep

        function rho = image_internal(self, g)
            if isempty(self.factors)
                rho = 1;
            else
                rho = self.factors{1}.image_internal(g);
                for i = 2:self.nFactors
                    rho = kron(rho, self.factors{i}.image_internal(g));
                end
            end
        end

        function rho = imageInverse_internal(self, g)
            if isempty(self.factors)
                rho = 1;
            else
                rho = self.factors{1}.imageInverse_internal(g);
                for i = 2:self.nFactors
                    rho = kron(rho, self.factors{i}.imageInverse_internal(g));
                end
            end
        end

    end
end
