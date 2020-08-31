classdef Character < replab.Obj
% Describes a character
%
% Example:
%   >>> D = replab.PermutationGroup.dihedral(6);
%   >>> rep = D.naturalRep.tensorPower(2);
%   >>> dec = rep.decomposition;
%   >>> irr = dec.irrep(4);
%   >>> c = replab.Character.fromApproximateRep(irr);
%   >>> g = [2 3 4 5 6 1];
%   >>> c.value(g)
%       -1
%   >>> c1 = c + c;
%   >>> c1.value(g)
%       -2
%   >>> c2 = c * c;
%   >>> c2.value(g)
%       1

    properties (SetAccess = protected)
        group % (`.FiniteGroup`): Group on which this class function is defined
        conjugacyClasses % (`.ConjugacyClasses`): List of conjugacy classes for which the order of `.values` is defined
        values % (`.cyclotomic`(1,\*)): Values of the character over the conjugacy classes
    end

    methods (Static)

        function c = fromApproximateRep(rep)
        % Computes an exact character from an approximate representation using Dixon's trick
        %
        % See J. D. Dixon, "Computing irreducible representations of groups", Math. Comp., vol. 24, no. 111, pp. 707–712, 1970
        % `<https://www.ams.org/mcom/1970-24-111/S0025-5718-1970-0280611-6/>`_
        %
        % Note that the function has complexity in $O(N e^2)$ where $N$ is the number of conjugacy classes and $e$ is the group
        % exponent.
        %
        % Args:
        %   rep (`.Rep`): Approximate representation
        %
        % Returns:
        %   `.Character`: Exact character
            group = rep.group;
            classes = group.conjugacyClasses;
            N = classes.nClasses;
            e = rep.group.exponent;
            m = zeros(N, e);
            zeta = exp(2i*pi/e);
            for i = 1:N
                r = classes.classes{i}.representative;
                for l = 1:e
                    v = 0;
                    for n = 0:e-1
                        v = v + trace(rep.image(group.composeN(r, n)))*zeta^(-l*n);
                    end
                    m(i, l) = v/e;
                end
            end
            m = round(m);
            zeta = replab.cyclotomic.E(e);
            values = replab.cyclotomic.zeros(1, N);
            for i = 1:N
                v = replab.cyclotomic.zeros(1, 1);
                for l = 1:e
                    v = v + m(i,l)*zeta^l;
                end
                values(i) = v;
            end
            c = replab.Character(classes, values);
        end

    end

    methods

        function self = Character(conjugacyClasses, values)
            self.group = conjugacyClasses.group;
            self.conjugacyClasses = conjugacyClasses;
            self.values = values;
        end

        function res = forClasses(self, newConjugacyClasses)
        % Returns the character with conjugacy classes reordered
            ind = self.conjugacyClasses.indicesOfClasses(newConjugacyClasses);
            newValues = self.values(ind);
            res = replab.Character(newConjugacyClasses, newValues);
        end

        function v = value(self, arg)
        % Returns the value of the class function over a group element or conjugacy class
        %
        % Args:
        %   element (element of `.group` or `.ConjugacyClass`): Element to compute the value of
        %
        % Returns:
        %   cyclotomic: Class function value
            if isa(arg, 'replab.ConjugacyClass')
                ind = self.conjugacyClasses.classIndexRepresentative(arg.representative);
            else
                ind = self.conjugacyClasses.classIndex(arg);
            end
            ind;
            v = self.values(ind);
        end

        function K = kernel(self)
        % Returns the kernel of this character
        %
        % The kernel of this character $\chi$ is defined as all $g \in G$ such that $\chi(g) = \chi(1)$ where $1$ is the identity.
        %
        % Returns:
        %   `.FiniteGroup`: A subgroup of `.group`
            idval = self.value(self.group.identity);
            mask = arrayfun(@(i) self.value(self.conjugacyClasses.classes{i}.representative) == idval, 1:self.conjugacyClasses.nClasses);
            K = self.conjugacyClasses.normalSubgroupClasses(find(mask));
        end

        function v = dotRep(self, rhs)
        % Computes the approximate inner product of this character with the character of an approximate representation
        %
        % Args:
        %   rhs (`.Rep`): Approximate representation
        %
        % Returns:
        %   double: Value of the dot product
            v = 0;
            n = self.conjugacyClasses.nClasses;
            sizes = self.conjugacyClasses.classSizes;
            chi = zeros(1, n);
            for i = 1:n
                cl = self.conjugacyClasses.classes{i};
                chi(i) = trace(rhs.image(cl.representative));
            end
            v = sum(conj(double(self.values)) .* chi .* double([sizes{:}]));
            v = v / double(self.group.order);
        end

        function v = dot(self, rhs)
        % Computes the inner product of this character with another character
        %
        % Args:
        %   rhs (`.Character`): Character
        %
        % Returns:
        %   `.cyclotomic`: Value of the dot product
            if self.conjugacyClasses.id ~= rhs.conjugacyClasses.id
                rhs = rhs.forClasses(self.conjugacyClasses);
            end
            sizes = self.conjugacyClasses.classSizes;
            sizes = [sizes{:}];
            v = sum(conj(self.values) .* rhs.values .* replab.cyclotomic.fromVPIs(sizes));
            v = v / replab.cyclotomic.fromVPIs(self.group.order);
        end

        function res = plus(lhs, rhs)
            if lhs.conjugacyClasses.id ~= rhs.conjugacyClasses.id
                rhs = rhs.forClasses(lhs.conjugacyClasses);
            end
            res = replab.Character(lhs.conjugacyClasses, lhs.values + rhs.values);
        end

        function res = mtimes(lhs, rhs)
            if lhs.conjugacyClasses.id ~= rhs.conjugacyClasses.id
                rhs = rhs.forClasses(lhs.conjugacyClasses);
            end
            res = replab.Character(lhs.conjugacyClasses, lhs.values .* rhs.values);
        end

        function res = eq(lhs, rhs)
            if lhs.conjugacyClasses.id ~= rhs.conjugacyClasses.id
                rhs = rhs.forClasses(lhs.conjugacyClasses);
            end
            res = all(lhs.values == rhs.values);
        end

        function res = ne(lhs, rhs)
            res = ~(lhs == rhs);
        end

    end

end