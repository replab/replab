classdef DirectSumRep < replab.Rep
% A direct sum of representations, such that images are diagonal by blocks

    properties
        factors % (cell(1,*) of `+replab.Rep`): Contained subrepresentations
    end

    methods

        function self = DirectSumRep(group, field, factors)
        % Constructs a direct sum from a cell array of representations
        %
        % All the subrepresentations should be defined on the same group, and on the same field.
        %
        % Args:
        %   group (`+replab.CompactGroup`): Common group
        %   field ({'R', 'C'}): Real or complex field
        %   blocks (cell(1,*) of `+replab.Rep`): Subrepresentations
            replab.rep.assertCompatibleFactors(group, field, factors);
            % own properties
            self.factors = factors;
            % replab.Rep immutable
            self.group = group;
            self.field = field;
            self.dimension = sum(cellfun(@(f) f.dimension, factors));
            % replab.Rep mutable
            factorsAreUnitary = cellfun(@(x) x.isUnitary, factors, 'uniform', 0);
            self.isUnitary = replab.trileanAnd(factorsAreUnitary{:});
        end

        function n = nFactors(self)
        % Returns the number of factors in the direct sum
        %
        % Returns:
        %   integer: Number of subrepresentations composing this representation
            n = length(self.factors);
        end

        function f = factor(self, i)
        % Returns a factor in the direct sum
        %
        % Args:
        %   i (integer): Index of a factor/block
        %
        % Returns:
        %   `+replab.Rep`: Representation corresponding to the i-th factor
            f = self.factors{i};
        end

        %% Str methods

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

        %% Rep methods

        function rho = image_internal(self, g)
            rhos = cellfun(@(rep) rep.image_internal(g), self.factors, 'uniform', 0);
            rho = blkdiag(rhos{:});
        end

        function rho = inverseImage_internal(self, g)
            rhos = cellfun(@(rep) rep.inverseImage_internal(g), self.factors, 'uniform', 0);
            rho = blkdiag(rhos{:});
        end

    end

end
