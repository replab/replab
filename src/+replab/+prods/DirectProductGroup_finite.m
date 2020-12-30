classdef DirectProductGroup_finite <  replab.DirectProductGroup & replab.NiceFiniteGroup
% External direct product of finite groups
%
% In particular, the permutation image of an element of a direct product group
% is simply the concatenation of the permutation images of the factors (which
% are nice finite groups themselves).
%
% We overload a bunch of methods to make sure we use the more efficient variants, that do not require the BSGS chain construction.

    methods

        function self = DirectProductGroup_finite(factors)
            assert(all(cellfun(@(x) isa(x, 'replab.FiniteGroup'), factors)));
            self.factors = factors;
            self.identity = cellfun(@(f) f.identity, factors, 'uniform', 0);
            % the generators of a direct product of finite groups is
            % the union of the generators of the factors, lifted into the
            % proper tuples
            generators = {};
            for i = 1:length(factors)
                factor = factors{i};
                for j = 1:factor.nGenerators
                    generator = self.identity;
                    generator{i} = factor.generator(j);
                    generators{1, end+1} = generator;
                end
            end
            self.generators = generators;
            if all(cellfun(@(f) f.type == f, factors))
                self.type = self;
            else
                self.type = replab.prods.DirectProductGroup_finite(cellfun(@(f) f.type, factors, 'uniform', 0));
            end
        end

        function res = hasSameTypeAs(self, rhs)
            lhs = self.type;
            rhs = rhs.type;
            if isa(rhs, 'replab.DirectProductGroup') && isa(rhs, 'replab.FiniteGroup') && lhs.nFactors == rhs.nFactors
                res = all(arrayfun(@(i) lhs.factor(i) == rhs.factor(i), 1:self.nFactors));
            else
                res = false;
            end
        end

    end

    methods (Access = protected) % Implementations

        function g = atFun(self, ind)
        % See comments in self.elements
            g = self.identity;
            ind = ind - 1;
            for i = self.nFactors:-1:1
                f = self.factor(i);
                this = mod(ind, f.order);
                ind = (ind - this)/f.order;
                g{i} = f.elements.at(this + 1);
            end
        end

        function ind = findFun(self, g)
        % See comments in self.elements
            ind = vpi(0);
            for i = 1:self.nFactors
                f = self.factor(i);
                ind = ind * f.order;
                ind = ind + f.elements.find(g{i}) - 1;
            end
            ind = ind + 1;
        end

        function o = computeOrder(self)
            o = vpi(1);
            % The order of a direct product is the product of the
            % order of the factors
            for i = 1:self.nFactors
                o = o * self.factor(i).order;
            end
        end

        function e = computeElements(self)
            e = replab.IndexedFamily.lambda(self.order, ...
                                            @(ind) self.atFun(ind), ...
                                            @(g) self.findFun(g));
            % The elements of a direct product of finite groups can be
            % enumerated by considering the direct product as a cartesian
            % product of sets, and decomposing the index a la ind2sub/sub2ind
            % which is the role of the `atFun` and `findFun` functions

            % TODO: verify ordering
        end

        function gd = computeDecomposition(self)
            T = {};
            % The decomposition of a direct product into sets
            % is simply the concatenation of the sequence of sets
            % corresponding to each factor
            for i = 1:self.nFactors
                D = self.factor(i).decomposition.T;
                Ti = cell(1, length(D));
                for j = 1:length(D)
                    Dj = D{j};
                    Tij = cell(1, length(Dj));
                    for k = 1:length(Dj)
                       Djk = Dj{k};
                       Tijk = self.identity;
                       Tijk{i} = Djk;
                       Tij{k} = Tijk;
                    end
                    Ti{j} = Tij;
                end
                T = horzcat(T, Ti);
            end
            gd = replab.FiniteGroupDecomposition(self, T);
        end

    end

    methods % Implementations

        % Str

        function names = hiddenFields(self)
            names = replab.str.uniqueNames( ...
                hiddenFields@replab.DirectProductGroup(self), ...
                hiddenFields@replab.NiceFiniteGroup(self) ...
                );
        end

        function [names values] = additionalFields(self)
            [names1 values1] = additionalFields@replab.DirectProductGroup(self);
            [names2 values2] = additionalFields@replab.NiceFiniteGroup(self);
            names = replab.str.horzcatForce(names1, names2);
            values = replab.str.horzcatForce(values1, values2);
        end

        % Domain

        function g = sample(self)
            g = sample@replab.DirectProductGroup(self);
        end

        % FiniteGroup

        function l = contains(self, g)
            for i = 1:self.nFactors
                if ~self.factor(i).contains(g{i})
                    l = false;
                    return
                end
            end
            l = true;
        end

        % NiceFiniteGroup

        function p = niceImage(self, g)
            shift = 0;
            p = [];
            % concatenates the permutation images of the factors
            for i = 1:self.nFactors
                pf = self.factor(i).niceMorphism.imageElement(g{i});
                p = [p pf+shift];
                shift = shift + length(pf);
            end
        end

    end

end
