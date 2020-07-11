classdef OfCompactGroups < replab.CompactGroup
% Describes an external direct product of compact groups

    properties (SetAccess = protected)
        factors % (factors{1,:} of subclass of `+replab.CompactGroup`): Factor groups
    end

    methods

        function self = OfCompactGroups(factors)
        % Constructs a direct product of groups
        %
        % Args:
        %   factors (row cell array): Factor groups
            n = length(factors);
            for i = 1:n
                assert(isa(factors{i}, self.requiredType), ['All factors must be instances of' self.requiredType]);
            end
            self.factors = factors;
            self.identity = cellfun(@(f) f.identity, factors, 'uniform', 0);
        end

        function t = requiredType(self)
            t = 'replab.CompactGroup';
        end

        function n = nFactors(self)
            n = length(self.factors);
        end

        function f = factor(self, i)
            f = self.factors{i};
        end

        function rep = directSumRep(self, factorReps)
        % Constructs a direct sum representation
        %
        % Args:
        %   factorReps (row cell array): Representations for each of the factor groups
        %                                i.e. factorReps{i} is a representation of factor(i)
        %
        % Returns:
        %   `+replab.Rep`: A direct sum representation
            rep = replab.directproduct.SumRep(self, factorReps);
        end

        function rep = directSumRepFun(self, fun)
        % Constructs a direct sum representation from a function that maps factors to representations
        %
        % Enables constructions such as ``directProduct.sumRepFun(@(x) x.definingRep)``
        %
        % Args:
        %   fun (function_handle): A function valid for each factor group that maps the group to one of
        %                          its representations
        %
        % Returns:
        %   replab.Rep: A direct sum representation
            reps = arrayfun(@(i) fun(self.factor(i)), 1:self.nFactors, 'uniform', 0);
            rep = self.directSumRep(reps);
        end

        function rep = directSumRepFunWithIndex(self, fun)
        % Constructs a direct sum representation from a function that maps factors to representations
        %
        % Enables constructions such as ``directProduct.sumRepFunWithIndex(@(x, i) x.definingRep)``
        %
        % Args:
        %   fun (function_handle): A function valid for each factor group that maps the group and its index
        %                          to one of its representations
        %
        % Returns:
        %   replab.Rep: A direct sum representation
            reps = arrayfun(@(i) fun(self.factor(i), i), 1:self.nFactors, 'uniform', 0);
            rep = self.directSumRep(reps);
        end

        function rep = tensorRep(self, factorReps)
        % Constructs a tensor product representation
        %
        % Args:
        %   factorReps (row cell array): Representations for each of the factor groups
        %                                i.e. factorReps{i} is a representation of factor(i)
        %
        % Returns:
        %   replab.Rep: A tensor representation
            rep = replab.directproduct.TensorRep(self, factorReps);
        end

        function rep = tensorRepFun(self, fun)
        % Constructs a direct sum representation from a function that maps factors to representations
        %
        % Enables constructions such as ``directProduct.tensorRepFun(@(x) x.definingRep)``
        %
        % Args:
        %   fun (function_handle): A function valid for each factor group that maps the group to one of
        %                          its representations
        %
        % Returns:
        %   replab.Rep: A tensor representation
            reps = arrayfun(@(i) fun(self.factor(i)), 1:self.nFactors, 'uniform', 0);
            rep = self.tensorRep(reps);
        end

        function rep = tensorRepFunWithIndex(self, fun)
        % Constructs a direct sum representation from a function that maps factors to representations
        %
        % Enables constructions such as ``directProduct.tensorRepFunWithIndex(@(x, i) x.definingRep)``
        %
        % Args:
        %   fun (function_handle): A function valid for each factor group that maps the group and its index
        %                          to one of its representations
        %
        % Returns:
        %   replab.Rep: A tensor representation
            reps = arrayfun(@(i) fun(self.factor(i), i), 1:self.nFactors, 'uniform', 0);
            rep = self.tensorRep(reps);
        end


        %% Str methods

        function names = hiddenFields(self)
            names = hiddenFields@replab.CompactGroup(self);
            names{1, end+1} = 'factors';
        end

        function [names values] = additionalFields(self)
            [names values] = additionalFields@replab.CompactGroup(self);
            for i = 1:self.nFactors
                names{1, end+1} = sprintf('factor(%d)', i);
                values{1, end+1} = self.factor(i);
            end
        end

        %% Domain methods

        function b = eqv(self, x, y)
            b = true;
            for i = 1:self.nFactors
                if ~self.factor(i).eqv(x{i}, y{i})
                    b = false;
                    return
                end
            end
        end

        function h = hash(self, x)
            hf = zeros(1, self.nFactors);
            for i = 1:self.nFactors
                hf(i) = self.factor(i).hash(x{i});
            end
            h = replab.Domain.hashVector(hf);
        end

        %% Monoid methods

        function z = compose(self, x, y)
            z = cell(1, self.nFactors);
            for i = 1:self.nFactors
                z{i} = self.factor(i).compose(x{i}, y{i});
            end
        end

        %% Group methods

        function xInv = inverse(self, x)
            xInv = cell(1, self.nFactors);
            for i = 1:self.nFactors
                xInv{i} = self.factor(i).inverse(x{i});
            end
        end

        %% CompactGroup methods

        function g = sample(self)
            g = cell(1, self.nFactors);
            for i = 1:self.nFactors
                g{i} = self.factor(i).sample;
            end
        end

    end

end
