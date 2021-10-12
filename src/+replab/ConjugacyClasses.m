classdef ConjugacyClasses < replab.Obj
% Stores information about the conjugacy classes of a group
%
% In particular, it can compute power maps.

    properties (Access = protected)
        primes % (integer(1,nP)): List of primes that divide the group order (computed as needed)
        primePowerMap % (integer(nP,nC)): Power maps (computed as needed)
    end

    properties (SetAccess = protected)
        group % (`.FiniteGroup`): Group whose conjugacy classes we describe
        classes % (cell(1,nC) of `.ConjugacyClass`): Conjugacy classes
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
                    m(j) = self.classIndexOf(self.group.composeN(self.classes{j}.representative, p));
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

    methods % Implementations

        % Str

        function [names values] = additionalFields(self)
            [names values] = additionalFields@replab.Obj(self);
            for i = 1:self.nClasses
                names{1, end+1} = sprintf('classes{%d}', i);
                values{1, end+1} = self.classes{i};
            end
        end

        function names = hiddenFields(self)
            names = hiddenFields@replab.Obj(self);
            names{1, end+1} = 'classes';
        end

        % Obj

        function l = laws(self)
            l = replab.laws.ConjugacyClassesLaws(self);
        end

    end

% $$$     methods (Static)
% $$$
% $$$         function C = sorted(group, classes)
% $$$             reps = cellfun(@(c) group.niceMorphism.imageElement(c.representative), classes, 'uniform', 0);
% $$$             repMat = cellfun(@transpose, reps, 'uniform', 0);
% $$$             repMat = [repMat{:}]';
% $$$             [~, I] = sortrows(repMat);
% $$$             C = replab.ConjugacyClasses(group, classes(I));
% $$$         end
% $$$
% $$$     end

    methods

        function self = ConjugacyClasses(group, classes)
            self.group = group;
            self.classes = classes;
            self.primes = zeros(1, 0);
            self.primePowerMap = zeros(0, length(classes));
        end

    end

    methods % Properties

        function s = centralizerSizes(self)
        % Returns the sizes of the centralizers
        %
        % Returns:
        %   cell(1,\*) of vpi: Size of the centralizer for each conjugacy class
            s = cellfun(@(c) c.representativeCentralizer.order, self.classes, 'uniform', 0);
        end

        function s = classElementOrders(self)
        % Returns the element order for each of the conjugacy classes
        %
        % Returns:
        %   integer(1,\*): Elements orders
            s = cellfun(@(c) c.elementOrder, self.classes);
        end

        function r = classRepresentatives(self)
        % Returns the canonical representatives of the conjugacy classes
        %
        % Returns:
        %   cell(1,\*) of `.ConjugacyClass`: Representatives
            r = cellfun(@(c) c.representative, self.classes, 'uniform', 0);
        end

        function s = classSizes(self)
        % Returns the sizes of the conjugacy classes
        %
        % Returns:
        %   cell(1,\*) of vpi: Size of each conjugacy class
            o = self.group.order;
            s = cellfun(@(c) c.nElements, self.classes, 'uniform', 0);
        end

        function n = nClasses(self)
        % Returns the number of conjugacy classes in the group
        %
        % Returns:
        %   integer: Number of classes
            n = length(self.classes);
        end

    end

    methods

        function ind = classIndexOf(self, g, varargin)
        % Finds the conjugacy class where a given element is located
        %
        % Args:
        %   g (element of `.group`): Group element
        %
        % Keyword Args:
        %   isCanonical (logical, optional): Whether the given ``g`` is the class canonical representative, default: false
        %
        % Returns:
        %   integer: Index of the class containing ``g``
            args = replab.util.populateStruct(struct('isCanonical', false), varargin);
            if args.isCanonical
                rep = g;
            else
                cc = self.group.conjugacyClass(g);
                rep = cc.representative;
            end
            for k = 1:length(self.classes)
                if self.group.eqv(self.classes{k}.representative, rep)
                    ind = k;
                    return
                end
            end
        end

        function ind = indicesOfClasses(self, classes)
        % Returns the indices of the given conjugacy classes in this list
        %
        % Args:
        %   classes (`.ConjugacyClasses`): Conjugacy classes of `.group`
        %
        % Returns:
        %   integer(1,\*): Indices of the given classes in these classes
            n = classes.nClasses;
            assert(self.nClasses == n);
            ind = arrayfun(@(i) self.classIndexOf(classes.classes{i}.representative, 'isCanonical', true), 1:n);
        end

        function G = normalSubgroupClasses(self, indices)
        % Returns the normal subgroup consisting of the conjugacy classes whose positions are given
        %
        % Args:
        %   indices (integer(1,\*)): Indices of conjugacy classes
        %
        % Returns:
        %   `.FiniteGroup`: Normal subgroup representing that union
            reps = arrayfun(@(ind) self.classes{ind}.representative, indices, 'uniform', 0);
            G = self.group.normalClosure(self.group.subgroup(reps));
        end

    end

    methods % Power maps

        function m = powerMap(self, n)
        % Returns the power map corresponding to the given exponent
        %
        % We have ``m(i) = self.classIndexOf(self.group.composeN(self.classes{i}.representative, n))``.
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

        function pp = powerMapDefaultPrimes(self)
        % Returns the list of primes that are necessary to reconstruct any power map
        %
        % Returns:
        %   integer(1,\*): Prime numbers that are less or equal to the largest class element order
            pp = primes(max(self.classElementOrders));
        end

        function a = powerMapMatrix(self)
        % Returns an adjacency-like matrix that describes (incomplete) relationships between conjugacy classes
            pp = self.powerMapDefaultPrimes;
            m = self.powerMaps(pp);
            a = zeros(self.nClasses, self.nClasses);
            for r = 1:self.nClasses
                powers = m(:,r);
                for c = 1:self.nClasses
                    a(r,c) = prod(pp(powers == c));
                end
            end
        end

        function m = powerMaps(self, exponents)
        % Returns or computes&stores the power map corresponding to the given exponent
        %
        % Args:
        %   exponents (integer(1,\*)): Exponents
        %
        % Returns:
        %   integer(length(exponents),\*): Index of the power of each conjugacy class
            m = zeros(length(exponents), self.nClasses);
            for i = 1:length(exponents)
                m(i,:) = self.powerMap(exponents(i));
            end
        end

    end

% $$$     methods
% $$$
% $$$         function c1 = imap(self, f)
% $$$         % Maps the conjugacy classes under an isomorphism
% $$$         %
% $$$         % Args:
% $$$         %   f (`.FiniteIsomorphism`): Isomorphism with ``self.group.isSubgroupOf(f.source)``
% $$$         %
% $$$         % Returns:
% $$$         %   `.ConjugacyClasses`: The conjugacy classes mapped under ``f``, expressed as a subset of ``f.image``
% $$$             if self.group.order < f.source.order
% $$$                 f = f.restrictedSource(self.group);
% $$$             end
% $$$             classes1 = cellfun(@(c) c.imap(f), self.classes, 'uniform', 0);
% $$$             c1 = replab.ConjugacyClasses(f.target, classes1);
% $$$         end
% $$$
% $$$     end

end
