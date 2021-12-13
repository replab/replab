classdef SemidirectProductGroup_finite < replab.SemidirectProductGroup & replab.gen.FiniteGroup

    methods

        function self = SemidirectProductGroup_finite(type, generators, nice, niceIsomorphism)
            self@replab.gen.FiniteGroup(type, generators, 'nice', nice, 'niceIsomorphism', niceIsomorphism);
            self.H = type.H;
            self.N = type.N;
            self.phi = type.phi;
        end

    end

% $$$     methods (Access = protected)
% $$$
% $$$         function g = atFun(self, ind)
% $$$             ind = ind - 1;
% $$$             indN = mod(ind, self.N.order);
% $$$             indH = (ind - indN)/self.N.order;
% $$$             g = {self.H.elementsSequence.at(indH + 1) self.N.elementsSequence.at(indN + 1)};
% $$$         end
% $$$ % $$$
% $$$         function ind = findFun(self, g)
% $$$             indH = self.H.elementsSequence.find(g{1});
% $$$             indN = self.N.elementsSequence.find(g{2});
% $$$             ind = (indH - 1)*self.N.order + indN;
% $$$         end
% $$$
% $$$     end

% $$$     methods (Access = protected) % Implementations
% $$$
% $$$         % FiniteGroup
% $$$
% $$$         function o = computeOrder(self)
% $$$             o = self.H.order * self.N.order;
% $$$         end
% $$$
% $$$         function e = computeElementsSequence(self)
% $$$             e = replab.Sequence.lambda(self.order, ...
% $$$                                        @(ind) self.atFun(ind), ...
% $$$                                        @(g) self.findFun(g));
% $$$         end
% $$$     end

    methods % Implementations

        % Str

        function s = headerStr(self)
            s = headerStr@replab.gen.FiniteGroup(self);
        end

        function [names, values] = additionalFields(self)
            [names, values] = additionalFields@replab.FiniteGroup(self);
        end

        function names = hiddenFields(self)
            names = hiddenFields@replab.gen.FiniteGroup(self);
        end

        % Obj

        function l = laws(self)
            l = laws@replab.FiniteGroup(self);
        end

        % Group

        function m = innerAutomorphism(self, by)
            m = innerAutomorphism@replab.FiniteGroup(self, by);
        end

        function m = isomorphismByFunctions(self, target, preimageElementFun, imageElementFun)
            m = isomorphismByFunctions@replab.FiniteGroup(self, target, preimageElementFun, imageElementFun);
        end

        function m = morphismByFunction(self, target, imageElementFun, torusMap)
            m = morphismByFunction@replab.FiniteGroup(self, target, imageElementFun, torusMap);
        end

        % Domain

        function b = eqv(self, x, y)
            b = eqv@replab.SemidirectProductGroup(self, x, y);
        end

        function s = sample(self)
            s = sample@replab.SemidirectProductGroup(self);
        end

        % Monoid

        function z = compose(self, x, y)
            z = compose@replab.SemidirectProductGroup(self, x, y);
        end

        % Group

        function xInv = inverse(self, x)
            xInv = inverse@replab.SemidirectProductGroup(self, x);
        end

        % FiniteSet

        function s = setProduct(self)
            TH = self.H.setProduct.sets;
            TN = self.N.setProduct.sets;
            idN = self.N.identity;
            idH = self.H.identity;
            TH1 = cellfun(@(t) cellfun(@(h) {h idN}, t, 'uniform', 0), TH, 'uniform', 0);
            TN1 = cellfun(@(t) cellfun(@(n) {idH n}, t, 'uniform', 0), TN, 'uniform', 0);
            s = replab.SetProduct(self, horzcat(TH1, TN1), true);
        end

        function G = withGeneratorNames(self, newNames)
            nice1 = self.nice.withGeneratorNames(newNames);
            G = replab.prods.SemidirectProduct_finite(self.type, self.generators, nice1, self.niceIsomorphism);
        end

    end

    methods % Bugfix for Octave method selection

        function b = isequal(lhs, rhs)
            b = isequal@replab.FiniteGroup(lhs, rhs);
        end

    end

    methods

        % Workaround for Octave bug

        function res = eq(self, rhs)
            res = eq@replab.gen.FiniteGroup(self, rhs);
        end

        function res = ne(self, rhs)
            res = ne@replab.gen.FiniteGroup(self, rhs);
        end

    end

end
