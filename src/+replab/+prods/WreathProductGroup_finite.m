classdef WreathProductGroup_finite < replab.WreathProductGroup & replab.gen.FiniteGroup
% Wreath product of a permutation group acting on a finite group

    methods

        function self = WreathProductGroup_finite(H, A, type, generators, nice, niceIsomorphism)
            assert(isa(H, 'replab.PermutationGroup'));
            assert(isa(A, 'replab.FiniteGroup'));
            self@replab.gen.FiniteGroup(type, generators, 'nice', nice, 'niceIsomorphism', niceIsomorphism);
            n = H.domainSize;
            base = A.directPower(n);
            phi = replab.perm.PermutationCellAction(H, base);
            self.H = H;
            self.N = base;
            self.phi = phi;
            self.A = A;
            self.n = n;
        end

    end

    methods % Implementation

        % Str

        function s = headerStr(self)
            s = headerStr@replab.gen.FiniteGroup(self);
        end

        function [names, values] = additionalFields(self)
            [names, values] = additionalFields@replab.FiniteGroup(self);
        end

        function names = hiddenFields(self)
            names = hiddenFields@replab.gen.FiniteGroup(self)
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

        function g = sample(self)
            g = {self.H.sample, arrayfun(@(i) self.A.sample, 1:self.n, 'uniform', 0)};
        end

        function b = eqv(self, x, y)
            b = self.type.eqv(x, y);
        end

        % Monoid

        function z = compose(self, x, y)
            z = self.type.compose(x, y);
        end

        % Group

        function xInv = inverse(self, x)
            xInv = self.type.inverse(x);
        end

        % FiniteGroup

        function G = withGeneratorNames(self, newNames)
            nice1 = self.nice.withGeneratorNames(newNames);
            G = replab.prods.WreathProductGroup_finite(self.H, self.A, self.generators, nice1, niceIsomorphism);
        end

    end

    methods (Static)

        function group = lambda(header, eqvFun, sampleFun, composeFun, identity, inverseFun)
            group = lambda@replab.FiniteGroup(header, eqvFun, sampleFun, composeFun, identity, inverseFun);
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
