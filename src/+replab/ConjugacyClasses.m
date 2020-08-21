classdef ConjugacyClasses < replab.Obj
% Stores information about the conjugacy classes of a group
%
% In particular, it can compute power maps.

    properties (SetAccess = protected)
        group % (`.FiniteGroup`): Group under study
        classes % (cell(1,nC) of `.ConjugacyClass`): Conjugacy classes
        primes % (integer(1,nP)): List of primes that divide the group order (computed as needed)
        primePowerMap % (integer(nP,nC)): Power maps (computed as needed)
    end

    methods (Access = protected)

        function m = getPrimePowerMap(self, p)
        % Returns or computes&stores the power map corresponding to the given exponent
        %
        % Args:
        %   p (integer): Exponent
        %
        % Returns:
        %   integer(1,\*): Index of the power of each conjugacy class
            i = find(self.primes == p, 1);
            if isempty(i)
                % compute the power map
                m = zeros(1, self.nClasses);
                for j = 1:self.nClasses
                    m(j) = self.classIndex(self.group.composeN(self.classes{j}.representative, p));
                end
                self.primePowerMap = [self.primePowerMap; m];
                self.primes = [self.primes p];
                [~, I] = sort(self.primes);
                self.primes = self.primes(I);
                self.primePowerMap = self.primePowerMap(I, :);
            else
                m = self.primePowerMap(i,:);
            end
        end

    end

    methods

        function self = ConjugacyClasses(group, classes)
            self.group = group;
            self.classes = classes;
            self.primes = zeros(1, 0);
            self.primePowerMap = zeros(0, length(classes));
        end

        function n = nClasses(self)
        % Returns the number of conjugacy classes in the group
        %
        % Returns:
        %   integer: Number of classes
            n = length(self.classes);
        end

        function s = centralizerSizes(self)
        % Returns the sizes of the centralizers
        %
        % Returns:
        %   integer(1,\*): Size of the centralizer for each conjugacy class
            s = cellfun(@(c) c.representativeCentralizer.order, self.classes, 'uniform', 0);
        end

        function s = classOrders(self)
        % Returns the element order for each of the conjugacy classes
        %
        % Returns:
        %   integer(1,\*): Elements orders
            s = cellfun(@(c) self.group.elementOrder(c.representative), self.classes, 'uniform', 0);
        end

        function s = classSizes(self)
        % Returns the sizes of the conjugacy classes
        %
        % Returns:
        %   integer(1,\*): Size of each conjugacy class
            o = self.group.order;
            s = cellfun(@(c) o/c.representativeCentralizer.order, self.classes, 'uniform', 0);
        end

        function ind = classIndex(self, g)
        % Finds the conjugacy class where a given element is located
        %
        % Args:
        %   g (element of `.group`): Group element
        %
        % Returns:
        %   integer: Index of the class containing ``g``
            r = replab.ConjugacyClass.representative(self.group, g);
            for k = 1:length(self.classes)
                if self.group.eqv(self.classes{k}.representative, r)
                    ind = k;
                    return
                end
            end
            assert(self.group.contains(g));
            error('Error in the list of conjugacy classes');
        end

        function m = powerMap(self, n)
        % Returns the power map corresponding to the given exponent
        %
        % We have ``m(i) = self.classIndex(self.group.composeN(self.classes{i}.representative, n))``.
        %
        % Args:
        %   p (integer): Exponent
        %
        % Returns:
        %   integer(1,\*): Index of the power of each conjugacy class
            m = 1:self.nClasses;
            for p = double(factor(n))
                mp = self.getPrimePowerMap(p);
                m = mp(m);
            end
        end

        function m = powerMaps(self, exponents)
        % Returns or computes&stores the power map corresponding to the given exponent
        %
        % Args:
        %   exponents (integer): Exponents
        %
        % Returns:
        %   integer(length(exponents),\*): Index of the power of each conjugacy class
            m = zeros(length(exponents), self.nClasses);
            for i = 1:length(exponents)
                m(i,:) = self.powerMap(exponents(i));
            end
        end

        function a = powerMapMatrix(self)
        % Returns an adjacency-like matrix that describes (incomplete) relationships between conjugacy classes
            primes = unique(double(factor(self.group.order)));
            m = self.powerMaps(primes);
            a = zeros(self.nClasses, self.nClasses);
            for r = 1:self.nClasses
                powers = m(:,r);
                for c = 1:self.nClasses
                    a(r,c) = prod(primes(powers == c));
                end
            end
        end

        function c1 = imap(self, f, imageGroup, preserveLexOrder)
        % Maps the conjugacy classes under an isomorphism
        %
        % Args:
        %   f (`.FiniteIsomorphism`): Isomorphism with ``self.group.isSubgroupOf(f.source)``
        %   imageGroup (`.FiniteGroup`, optional): Image of `.group` under ``f``, default ``[]`` (recompute)
        %   preserveLexOrder (logical, optional): Whether the isomorphism preserves the lexicographic order of group elements, default false
        %
        % Returns:
        %   `.ConjugacyClasses`: The conjugacy classes mapped under ``f``, expressed as a subset of ``f.image``
            if nargin < 3 || isempty(imageGroup)
                imageGroup = f.imageGroup(self.group);
            end
            if nargin < 4 || isempty(preserveLexOrder)
                preserveLexOrder = false;
            end
            classes1 = cellfun(@(c) c.imap(f, imageGroup, preserveLexOrder), self.classes, 'uniform', 0);
            c1 = replab.ConjugacyClasses(imageGroup, classes1);
        end

    end

end
