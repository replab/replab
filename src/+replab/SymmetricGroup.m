classdef SymmetricGroup < replab.PermutationGroup
% The symmetric group ``S(n)``, i.e. all permutations over n = "domainSize" elements
%
% Example:
%   >>> S5 = replab.S(5);
%   >>> S5.order
%      ans =
%      120

    methods (Static)

        function G = make(n)
        % Constructs the symmetric over a given domain size
        %
        % This static method keeps the constructed copies of ``S(n)`` in cache.
        %
        % Args:
        %   domainSize (integer): Domain size, must be > 0
        %
        % Returns:
        %   `.SymmetricGroup`: The constructed or cached symmetric group
            persistent cache
            if isempty(cache)
                cache = cell(1, 0);
            end
            if n+1 > length(cache) || isempty(cache{n+1})
                cache{1,n+1} = replab.SymmetricGroup(n);
            end
            G = cache{n+1};
        end

    end

    methods (Access = protected) % Constructor

        function self = SymmetricGroup(domainSize)
        % Constructs the symmetric over a given domain size
        %
        % Instead of the constructor, use `.make`, which caches the constructed group.
        %
        % Args:
        %   domainSize (integer): Domain size, must be >= 0
            if domainSize < 2
                generators = cell(1, 0);
            elseif domainSize == 2
                generators = {[2 1]};
            else
                generators = {[2:domainSize 1] [2 1 3:domainSize]};
            end
            self@replab.PermutationGroup(domainSize, generators, 'order', @() replab.util.factorial(domainSize));
        end

        function b = hasFastOrder(self)
            b = true;
        end

    end

    methods % Implementations

        % replab.Str

        function s = headerStr(self)
            s = sprintf('Symmetric group acting on %d elements', self.domainSize);
        end

        % replab.Domain

        function s = sample(self)
            s = randperm(self.domainSize); % overriden for efficiency
        end

        % replab.FiniteGroup

        function b = contains(self, g)
            assert(length(g) == self.domainSize, 'Permutation in wrong domain');
            assert(all(g > 0), 'Permutation should have positive coefficients');
            b = true;
        end

        % replab.PermutationGroup

        function c = conjugacyClasses(self)
            c = self.cached('conjugacyClasses', @() self.computeConjugacyClasses);
        end

    end

    methods (Access = protected)

        function classes = computeConjugacyClasses(self)
            I = replab.sym.IntegerPartition.all(self.domainSize);
            classes = replab.ConjugacyClasses(self, cellfun(@(ip) ip.conjugacyClass, I, 'uniform', 0));
            classes = classes.sorted;
        end

% $$$         function E = computeElementsSequence(self)
% $$$             E = replab.Sequence.lambda(self.order, ...
% $$$                                        @(ind) self.enumeratorAt(ind), ...
% $$$                                        @(el) self.enumeratorFind(el));
% $$$         end
% $$$
% $$$         function ind = enumeratorFind(self, g)
% $$$             n = self.domainSize;
% $$$             ind0 = vpi(0);
% $$$             els = 1:n;
% $$$             for i = 1:n
% $$$                 ind0 = ind0 * (n - i + 1);
% $$$                 ind0 = ind0 + (find(els == g(i)) - 1);
% $$$                 els = setdiff(els, g(i));
% $$$             end
% $$$             ind = ind0 + 1;
% $$$         end
% $$$
% $$$         function g = enumeratorAt(self, ind)
% $$$             n = self.domainSize;
% $$$             ind0 = ind - 1; % make it 0-based
% $$$             inds = zeros(1, n);
% $$$             for i = 1:n
% $$$                 r = mod(ind0, i);
% $$$                 ind0 = (ind0 - r)/i;
% $$$                 inds(i) = double(r + 1);
% $$$             end
% $$$             inds = fliplr(inds);
% $$$             els = 1:n;
% $$$             g = zeros(1, n);
% $$$             for i = 1:n
% $$$                 e = els(inds(i));
% $$$                 g(i) = e;
% $$$                 els = setdiff(els, e);
% $$$             end
% $$$         end

    end

    methods % Representations

        function [rep data] = irrep(self, partition, form)
        % Returns the irreducible representation of this symmetric group corresponding the given Young Diagram
        %
        % Example:
        %   >>> S5 = replab.S(5); % doctest: +cyclotomic
        %   >>> rep = S5.irrep([3 2], 'specht');
        %   >>> rep.dimension
        %         5
        %   >>> rep = S5.irrep([3 2], 'seminormal');
        %   >>> rep.dimension
        %         5
        %   >>> rep = S5.irrep([3 2], 'orthogonal');
        %   >>> rep.dimension
        %         5
        %
        % Args:
        %   partition (integer(1,\*)): The partition corresponding the Young Diagram, with elements listed in
        %                              decreasing order (e.g ``[3 3 1]`` represents the partition of 7 elements: ``7 = 3+3+1``
        %   form ('specht', 'seminormal', or 'orthogonal'): The form the irrep takes. Default is 'specht'.
        %
        % Note: 'specht' is slower to construct than 'seminormal' or 'orthogonal' but, unlike them, has integer entries.
        %
        % Returns:
        %   `+replab.Rep`: The corresponding irreducible representation
            assert(sum(partition) == self.domainSize);
            if nargin == 2
                form = 'specht';
            end
            partition = sort(partition,'descend');
            switch form
              case 'specht'
                data = replab.sym.SymmetricSpechtIrrep(self, partition);
              case 'seminormal'
                data = replab.sym.SymmetricYoungIrrep(self, partition, 'seminormal');
              case 'orthogonal'
                data = replab.sym.SymmetricYoungIrrep(self, partition, 'orthogonal');
              otherwise
                error('That is not a valid irreducible representation form')
            end
            rep = data.rep;
        end

    end

    methods (Static) % Representations

        function dim = irrepDimension(partition)
        % Returns the dimension of the irreducible representation of a symmetric group corresponding a Young Diagram
        %
        % Example:
        %   >>> S5 = replab.S(5);
        %   >>> rep = S5.irrep([3 2]);
        %   >>> rep.dimension
        %         5
        %   >>> S5.irrepDimension([3 2])
        %         5
        %
        % Args:
        %   partition (integer(1,\*)): The partition corresponding the Young Diagram, with elements listed
        %                              in decreasing order (e.g ``[5 3 1]`` represents the partition of 9 elements ``9 = 5+3+1``
        % Returns:
        %   integer: The dimension of the corresponding irreducible representation.
            dim = replab.sym.findPartitions.dimension(partition);
        end

    end

end
