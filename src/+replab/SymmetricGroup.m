classdef SymmetricGroup < replab.PermutationGroup
% Describes permutations over n = "domainSize" elements, i.e. the symmetric group Sn
%
% Example:
%   >>> S5 = replab.SymmetricGroup(5);
%   >>> S5.order
%      ans =
%      120

    methods

        function self = SymmetricGroup(domainSize)
        % Constructs the symmetric over a given domain size
        %
        % Args:
        %   domainSize (integer): Domain size, must be > 0
            o = factorial(vpi(domainSize));
            if domainSize < 2
                generators = cell(1, 0);
            elseif domainSize == 2
                generators = {[2 1]};
            else
                generators = {[2:domainSize 1] [2 1 3:domainSize]};
            end
            self = self@replab.PermutationGroup(domainSize, generators, o, 'self');
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

    end

    methods % Element creation methods

        function p = transposition(self, i, j)
        % Returns the transposition permuting ``i`` and ``j``.
        %
        % Args:
        %   i (integer): First domain element to be transposed.
        %   j (integer): Second domain element to be transposed.
        %
        % Returns:
        %   permutation: The constructed transposition
            n = self.domainSize;
            assert(1 <= i);
            assert(i <= n);
            assert(1 <= j);
            assert(j <= n);
            assert(i ~= j);
            p = 1:n;
            p([i j]) = [j i];
        end

        function p = shift(self, i)
        % Returns the cyclic permutation that shifts the domain indices by ``i``.
        %
        % Args:
        %   i: Shift so that ``j`` is sent to ``j + i`` (wrapping around).
        %
        % Returns:
        %   permutation: The constructed cyclic shift
            p = mod((0:n-1)+i, n)+1;
        end

        function p = fromCycles(self, varargin)
        % Constructs a permutation from a product of cycles.
        %
        % Each cycle is given as a row vector, and the sequence of cycles is given as variable arguments.
        %
        % Args:
        %   varargin (cell(1,\*) of integer(1,\*)): Sequence of cycles as row vectors of indices
        %
        % Returns:
        %   permutation: The permutation corresponding to the product of cycles.
            n = self.domainSize;
            p = 1:n;
            for i = length(varargin):-1:1
                cycle = varargin{i};
                % cycle 2 3 1 means that 2 -> 3, 3 -> 1, 1 -> 2
                cycleImage = [cycle(2:end) cycle(1)];
                newEl = 1:n;
                newEl(cycle) = cycleImage;
                p = newEl(p); % compose(newEl, p);
            end
        end

    end

    methods % Property computation

        function o = computeOrder(self)
            o = factorial(vpi(self.domainSize));
        end

        function E = computeElements(self)
            E = replab.IndexedFamily.lambda(self.order, ...
                                            @(ind) self.enumeratorAt(ind), ...
                                            @(el) self.enumeratorFind(el));
        end

        function d = computeDecomposition(self)
            G = self.subgroup(self.generators, self.order);
            d = G.decomposition;
        end

    end

    methods (Access = protected)


        function ind = enumeratorFind(self, g)
            n = self.domainSize;
            ind0 = vpi(0);
            els = 1:n;
            for i = 1:n
                ind0 = ind0 * (n - i + 1);
                ind0 = ind0 + (find(els == g(i)) - 1);
                els = setdiff(els, g(i));
            end
            ind = ind0 + 1;
        end

        function g = enumeratorAt(self, ind)
            n = self.domainSize;
            ind0 = ind - 1; % make it 0-based
            inds = zeros(1, n);
            for i = 1:n
                r = mod(ind0, i);
                ind0 = (ind0 - r)/i;
                inds(i) = double(r + 1);
            end
            inds = fliplr(inds);
            els = 1:n;
            g = zeros(1, n);
            for i = 1:n
                e = els(inds(i));
                g(i) = e;
                els = setdiff(els, e);
            end
        end

    end

    methods % Representations

        function rep = irrep(self, partition)
        % Returns the irreducible representation of this symmetric group corresponding the given Young Diagram
        %
        % Args:
        %   partition (integer(1,\*)): The partition corresponding the Young Diagram, with elements listed in
        %                              decreasing order (e.g ``[3 3 1]`` represents the partition of 7 elements: ``7 = 3+3+1``
        %
        % Returns:
        %   `+replab.Rep`: The corresponding irreducible representation
            assert(sum(partition) == self.domainSize);
            rep = replab.sym.SymmetricGroupIrrep(self, partition);
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
            dim = replab.sym.SymmetricGroupIrrep.dimension(partition);
        end

    end

end
