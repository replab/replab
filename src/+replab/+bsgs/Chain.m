classdef Chain < replab.Str
% A BSGS chain data structure for a permutation group
%
% The represented group acts on $\{1, ..., n\}$. We follow loosely the
% notation in Handbook of Computational Group Theory, Derek Holt.
%
% The BSGS chain is stored as follows.
%
% We write $B = [\beta_1, \beta_2 ... \beta_k]$ where $k$ is the length of the stabilizer chain.
% We write $G^i = \{ g \in G : g(\beta_1) = \beta_1, ..., g(\beta_{i-1}) = \beta_{i-1} \}$, so that
% $G^1 = G$.
% The (currently known) strong generators are stored as column vectors in a matrix S.
%
% We write the set $S^i = S \cap G^i$, and the subgroups $H^i = <S^i>$, with $H^{k+1} = <1>$.
%
% The order of strong generators as columns in the matrix is such that
% $S^i =$ ``S(:, Sind(i):end)`` where ``Sind`` is a row vector of starting indices in ``S``.
%
% The properties ``Delta`` ($\Delta$), ``U``, ``Uinv`` are cell arrays indexed by ``i``, while ``iDelta``
% is stored as an integer matrix.
%
% For each level $i$ in the BSGS chain, we have the following.
%
% $\Delta^i = H^i(\beta_i) = \{ h(\beta_i) : h \in H^i \}$ is the orbit of $\beta_i$ under the current
% strong generators. $\Delta^i$ is stored in ``Delta{i}`` in the order of discovery with $\beta_i$ is sorted
% in ``Delta{i}(1)``.
%
% To check quickly if a point $b$ is a member of an orbit $\Delta^i$, the matrix ``iDelta`` stores its index
% in ``Delta{i}`` as follows: ``Delta{i}(iDelta(i,b)) == b``, and ``iDelta(i,b) == 0`` when $b$ is not part of
% $\Delta^i$.
%
% We define the transversal elements $u^i_b$ such that $u^i_b(\beta_i) = b$, and the inverse transversal
% elements $uInv^i_(\beta_i) = b$.
%
% The transversal elements for a particular level $i$ in the chain are stored as column vectors in the matrix ``U{i}``
% with an order corresponding to the orbit $\Delta^i$.
%
% Accordingly, transversal element inverses are stored as column vectors in a matrix in the matrix ``Uinv{i}``.
%
% The following invariants are maintained by the code below.
%
% When the chain is immutable: the strong generating set must be a strong generating set for the $G^i$,
% not only for the $H^i$, i.e. the BSGS construction is complete.
%
% When the chain is mutable:
%
% * The chain describes the orbits/transversals of $H^i$.
%
% In both cases:
%
% * The base ``B``, the data structures for orbits and transversals ``Delta``, ``U``, ``Uinv`` are consistent.
%
% * The strong generating set is ordered, so that the starting indices ``Sind`` describe the sequence $S^i$

    properties (SetAccess = protected)
        isMutable % (logical): Whether the chain can be modified
        n % (integer): Domain size
        B % (integer(1,\*)): Row vector of base points (all between 1..n without duplicates)

        S % (integer(n, nS)): Matrix of strong generators stored as columns vectors (nS = # of strong generators)
        Sind % (integer(1, k+1)): Starting index of strong generators for each stabilizer subgroup, where k = length(B)
             %
             %                    We have ``Sind(k+1) == nS + 1``.

        Delta % (cell(1,k) of integer(1,\*)): Each orbit is a row vector containing orbit elements
        iDelta % (integer(n,k)): For each orbit, maps a domain element to its position in Delta{i} or 0 if not present

        U % (cell(1,k) of integer(n,\*)): Transversal elements
        Uinv % (cell(1,k) of integer(n,\*)): Inverse transversal elements
    end

    methods

        function self = Chain(n, B, S, Sind, Delta, iDelta, U, Uinv)
        % Constructs an empty mutable chain for a group of permutations acting on ``n`` elements
        %
        % Args:
        %   n (integer): Domain size
        %                               The default value is `+replab.+bsgs.TrivialGroup`
        %
        % Returns:
        %   `+replab.bsgs.Chain`: A constructed BSGS chain
            self.isMutable = true;
            self.n = n;
            if nargin == 8
                self.B = B;
                self.S = S;
                self.Sind = Sind;
                self.Delta = Delta;
                self.iDelta = iDelta;
                self.U = U;
                self.Uinv = Uinv;
            else
                self.B = [];
                self.S = [];
                self.Sind = [1];
                self.Delta = {};
                self.iDelta = zeros(n,0);
                self.U = {};
                self.Uinv = {};
            end
        end

        function b = base(self)
            b = self.B;
        end

        function s = strongGeneratorsForLevel(self, l)
        % Returns the strong generators that are present in H^l but not H^l+1
            s = self.S(:,self.Sind(l):self.Sind(l+1)-1);
        end

        function replaceStrongGeneratorsForLevel(self, l, snew)
            sold = self.strongGeneratorsForLevel(l);
            Sind = self.Sind;
            self.S = [self.S(:,1:self.Sind(l)-1) snew self.S(:,self.Sind(l+1):end)];
            Sind(l+1:end) = Sind(l+1:end) + size(snew, 2) - size(sold, 2);
            self.Sind = Sind;
        end

        function c = baseChange(self, newBase)
            c = replab.bsgs.Chain.make(self.n, self.S, newBase, self.order);
        end

        function show(self, i)
            if nargin < 2
                table = cell(2, self.length+1);
                table{1,1} = 'base: ';
                table{2,1} = 'orbit: ';
                table{3,1} = 'sgens: ';
                spec = ['r' repmat('l', 1, self.length)];
                for i = 1:self.length
                    table{1,i+1} = sprintf('%d', self.B(i));
                    table{2,i+1} = strrep(replab.shortStr(sort(self.Delta{i})), ' ', '');
                    for j = self.Sind(i):self.Sind(i+1)-1
                        table{3+j-self.Sind(i),i+1} = strrep(replab.shortStr(self.S(:,j)'), ' ', '');
                    end
                end
                disp(strjoin(replab.str.alignspace(table, spec), '\n'));
            else
                fprintf('Level %d Base point %d\n', i, self.B(i));
                table = cell(length(self.Delta{i})+1, 3);
                table{1,1} = 'Orbit pt ';
                table{1,2} = 'Transversal ';
                table{1,3} = 'Trans. inv. ';
                for j = 1:length(self.Delta{i})
                    table{j+1,1} = sprintf('%d', self.Delta{i}(j));
                    table{j+1,2} = strrep(replab.shortStr(self.U{i}(:,j)'), ' ', '');
                    table{j+1,3} = strrep(replab.shortStr(self.Uinv{i}(:,j)'), ' ', '');
                end
                disp(strjoin(replab.str.alignspace(table, 'cll'), '\n'));
            end
        end

        function c = stabilizer(self, b)
            if self.B(1) == b
                newB = self.B(2:end);
                newS = self.S(:, self.Sind(2):end);
                newSind = self.Sind(2:end) - self.Sind(2) + 1;
                newDelta = self.Delta(2:end);
                newiDelta = self.iDelta(:, 2:end);
                newU = self.U(2:end);
                newUinv = self.Uinv(2:end);
                c = replab.bsgs.Chain(self.n, newB, newS, newSind, newDelta, newiDelta, newU, newUinv);
            else
                if ismember(b, self.B)
                    newBase = [b setdiff(self.B, b)];
                else
                    newBase = [b self.B];
                end
                c = self.baseChange(newBase).stabilizer(b);
            end
        end

        function checkLevel(self, i)
        % Checks the ``i``-th level of this BSGS chain
            n = self.n;
            betai = self.B(i);
            % check strong generators
            for j = self.Sind(i):self.Sind(i+1)-1
                % strong generators particular to this step move beta_i
                assert(self.S(betai, j) ~= betai);
            end
            for j = self.Sind(i+1):self.Sind(end)-1
                % strong generators of stabilizer subgroups are stabilized
                assert(self.S(betai, j) == betai);
            end
            iD = self.iDelta(:,i)';
            D = self.Delta{i};
            assert(isequal(sort(find(iD ~= 0)), sort(D)));
            assert(isequal(iD(D), 1:length(D)));
            for j = 1:length(D)
                b = D(j);
                Ui = self.U{i};
                ub1 = Ui(:,j)';
                ub2 = self.u(i, b);
                assert(isequal(ub1, ub2), 'inconsistent transversal retrieval');
                assert(ub1(betai) == b, 'inconsistent transversal element');
                Uinvi = self.Uinv{i};
                ubinv1 = Uinvi(:,j)';
                ubinv2 = self.uinv(i, b);
                assert(isequal(ubinv1, ubinv2), 'inconsistent transversal retrieval');
                assert(ubinv1(b) == betai, 'inconsistent transversal element');
                for l = self.Sind(i):self.Sind(end)-1
                    imgD = self.S(D, l);
                    assert(all(ismember(imgD, D)) > 0, 'All images under strong generators should be present in the orbit');
                end
            end
        end


        function check(self)
        % Checks that the BSGS chain is well-formed
            for i = 1:self.length
                self.checkLevel(i);
            end
        end

        %% Immutable functions

        function s = orbitSize(self, l)
        % Returns the size of the l-orbit orbit
            s = length(self.Delta{l});
        end

        function s = orbitSizes(self)
        % Returns the size of orbits for each level
            s = cellfun(@(x) length(x), self.Delta);
        end

        function k = length(self)
        % Returns the length of this BSGS chain
        %
        % Returns:
        %   integer: Chain length
            k = length(self.B);
        end

        function o = order(self)
        % Returns the order of this BSGS chain
        %
        % Returns:
        %   vpi: Size of the group stored in the chain
            if self.length == 0
                o = vpi(1);
                return
            end
            if exist('java.math.BigInteger')
                o = java.math.BigInteger.valueOf(self.orbitSize(1));
                for i = 2:self.length
                    o = o.multiply(java.math.BigInteger.valueOf(self.orbitSize(i)));
                end
                o = vpi(char(o.toString));
            else
                o = self.orbitSize(1);
                i = 2;
                while i <= self.length && log2(o) + log2(self.orbitSize(i)) < 53
                    o = o * length(self.orbitSize(i));
                    i = i + 1;
                end
                o = vpi(o);
                while i <= self.length
                    o = o * vpi(self.orbitSize(i));
                end
            end
        end

        function g = sampleUniformly(self)
        % Samples an element uniformly from the group
        %
        % Returns:
        %  permutation: Random group element
            g = 1:self.n;
            for i = 1:self.length
                gi = self.randomTransversal(i);
                g = g(gi); % compose(g, gi)
            end
        end

        function nS = nStrongGenerators(self)
        % Returns the number of strong generators in this BSGS chain
            nS = self.Sind(end) - 1;
        end

        function p = strongGenerator(self, i)
        % Returns the i-th strong generator
        %
        % Args:
        %   i (integer): Strong generator index
        %
        % Returns:
        %   A permutation given as a row vector
            p = self.S(:, i)';
        end

        function u = randomTransversal(self, i)
        % Returns a random transversal element from a specific transversal set
        %
        % Args:
        %  i (integer): Index of the transversal set
        %
        % Returns:
        %   permutation: Random transversal element
            j = randi(length(self.Delta{i}));
            Ui = self.U{i};
            u = Ui(:,j)';
        end

        function j = orbitIndex(self, i, b)
        % Looks up an orbit element
        %
        % Args:
        %   i (integer): Index of orbit
        %   b (integer): Orbit element to lookup
        %
        % Returns:
        %   integer: When ``b`` is part of the i-th orbit, returns an integer ``j`` such that ``self.Delta{i}(j) = b``,
        %            otherwise returns 0.
            j = self.iDelta(b, i);
        end

        function g = u(self, i, b)
        % Looks up a transversal element that maps beta_i to b
        %
        % Args:
        %   i (integer): Index of transversal
        %   b (integer): Orbit element
        %
        % Returns:
        %   The corresponding transversal element, or [] if b is not part of the orbit
            j = self.iDelta(b, i);
            if j == 0
                g = [];
            else
                Ui = self.U{i};
                g = Ui(:,j)';
            end
        end

        function g = uinv(self, i, b)
        % Looks up the inverse transversal element that maps b to beta_i
        %
        % Args:
        %   i (integer): Index of transversal
        %   b (integer): Orbit element
        %
        % Returns:
        %   The corresponding inverse transversal element or [] if b is not part of the orbit
            j = self.iDelta(b, i);
            if j == 0
                g = [];
                return
            end
            Uinvi = self.Uinv{i};
            g = Uinvi(:,j)';
        end

        function b = contains(self, g)
        % Tests whether the BSGS chain contains an element
            k = self.length;
            n = self.n;
            [h, i] = self.strip(g);
            b = (i > k) && all(h == 1:n); % i > k and h is the identity
        end


        function [h i] = strip(self, g)
        % Strips a permutation through the chain
        %
        % Args:
        %   g (row permutation vector): A permutation group element
        %
        % Returns
        % -------
        %   h: row permutation vector
        %     The part of the group element that could not be sifted
        %   i: integer
        %     Index ``i`` such that ``h(beta_i)`` was not part of the orbit ``Delta^i``
            k = self.length;
            h = g;
            for i = 1:k
                b = h(self.B(i));
                j = self.iDelta(b, i);
                if j == 0
                    return
                end
                Uinvi = self.Uinv{i};
                uinv = Uinvi(:, j)';
                % note order is reversed compared to Holt, as
                % we use a left action
                h = uinv(h); % compose(uinv, h)
            end
            i = k + 1; % marker that we striped through the chain
        end

        %% Element indexing

        function g = elementFromIndices(self, indices)
        % Computes the group element from transversal indices
        %
        % Args:
        %   indices (integer(1, \*)): Transversal indices
        %
        % Returns:
        %   permutation: Chain element
            g = 1:self.n;
            for i = 1:self.length
                Ui = self.U{i};
                gi = Ui(:,indices(i));
                g = g(gi); % compose(g, gi)
            end
        end

        function indices = indicesFromElement(self, g)
        % Computes the transversal indices decomposition for a group element
        %
        % The indices are such that
        % ``g = self.u(1, indices(1)) * ... * self.u(k, indices(k))``
        %
        % Args:
        %   g (permutation): A permutation group element
        %
        % Returns:
        %   integer(1,\*): Transversal indices
            k = self.length;
            h = g;
            indices = zeros(1, k);
            for i = 1:k
                b = h(self.B(i));
                j = self.iDelta(b, i);
                indices(i) = j;
                if j == 0
                    indices = [];
                    return
                end
                Uinvi = self.Uinv{i};
                uinv = Uinvi(:, j)';
                % note order is reversed compared to Holt, as
                % we use a left action
                h = uinv(h); % compose(uinv, h)
            end
        end

        %% Mutable methods

        function insertInOrbit(self, i, b, u)
        % Inserts a new orbit element in an orbit
        %
        % Modifies the chain in place.
        %
        % Args:
        %   i (integer): Level of the orbit, 1 <= i <= self.length
        %   b (integer): New orbit point
        %   u (permutation): New transversal element
            assert(self.isMutable);
            assert(self.iDelta(b, i) == 0, 'Orbit point must be new');
            self.Delta{i} = [self.Delta{i} b];
            idx = length(self.Delta{i});
            self.iDelta(b, i) = idx;
            self.U{i} = [self.U{i} u(:)];
            n = self.n;
            uinv = zeros(1, n);
            uinv(u) = 1:n;
            self.Uinv{i} = [self.Uinv{i} uinv(:)];
        end

        function completeOrbit(self, i)
        % Completes the i-th orbit
        %
        % Iterates over current orbit elements and strong generators, and add new orbit points if necessary
        %
        % Args:
        %   i (integer): Level of the orbit, ``1 <= i <= self.length``
            assert(self.isMutable);
            n = self.n;
            toTest = self.Delta{i}; % we need to test all elements in the currently known orbit
            Srange = self.Sind(i):self.Sind(end)-1;
            while ~isempty(toTest)
                imgs = self.S(toTest, Srange); % images of the tested point for all the strong generators
                iDelta = self.iDelta(:, i);
                mask = iDelta(imgs) == 0;
                mask = reshape(mask, size(imgs)); % in case imgs is a vector, Matlab may change btw row/column vectors
                [b_ind S_ind] = find(mask);
                for j = 1:length(b_ind)
                    s = self.S(:, Srange(S_ind(j))); % strong generator
                    b = toTest(b_ind(j)); % orbit element
                    newb = s(b);
                    if self.iDelta(newb, i) == 0
                        % new orbit point discovered, add it
                        newu = s(self.u(i, b)); % new transversal compose(s, self.u(i, b))
                        self.insertInOrbit(i, newb, newu);
                    end
                end
                toTest = unique(imgs(mask));
            end
        end

        function insertEndBasePoint(self, newBeta)
        % Adds a new basis point at the end of the BSGS chain
        %
        % Args:
        %   newBeta (integer): New basis point not part of the current basis
            assert(self.isMutable, 'Chain needs to be mutable');
            assert(~ismember(newBeta, self.B), 'Base point already exists');
            n = self.n;
            k = self.length; % previous chain length
            self.B = [self.B newBeta];
            % We add the data structures for the new orbit/transversal
            self.Delta = horzcat(self.Delta, newBeta);
            self.iDelta = [self.iDelta zeros(n, 1)];
            self.iDelta(newBeta, k+1) = 1;
            self.U = horzcat(self.U, {[1:n]'});
            self.Uinv = horzcat(self.Uinv, {[1:n]'});
            self.Sind = [self.Sind self.Sind(end)];
        end

        function addStrongGenerator(self, i, newS)
        % Adds a strong generator at a particular place in the BSGS chain
        %
        % Args:
        %   i (integer): Smallest i such that the strong generator ``newS`` is part of S^(i)
        %   newS (permutation): New strong generator
            self.S = [self.S(:, 1:(self.Sind(i+1)-1)) newS(:) self.S(:, self.Sind(i+1):end)];
            self.Sind((i+1):end) = self.Sind((i+1):end) + 1;
            for j = 1:i
                self.completeOrbit(j);
            end
        end

        function newSG = stripAndAddStrongGenerator(self, g)
        % Strips and, when relevant, adds a new strong generator to this BSGS chain
        %
        % Performs part of the RANDOMSCHREIER procedure in Holt et al., Handbook of CGT, page 98
        %
        % Args:
        %   g (row permutation vector): Element to strip
        %
        % Returns:
        %   logical: True if a new strong generator has been found
            n = self.n;
            [h j] = self.strip(g);
            newSG = false;
            if j <= self.length
                % New strong generator h at level j
                newSG = true;
            else
                % Check if h is the identity
                gamma = find(h ~= 1:n, 1);
                if length(gamma) > 0
                    % New strong generator h fixes all base points
                    % We have a new base point gamma
                    % and insert it at the end of the chain
                    self.insertEndBasePoint(gamma);
                    newSG = true;
                end
            end
            if newSG
                % if we have a new strong generator, add it to the chain
                self.addStrongGenerator(j, h);
            end
        end

        function makeImmutable(self)
            assert(self.isMutable);
            self.isMutable = false;
        end

        function randomizedSchreierSims(self, order)
        % Runs the randomized Schreier-Sims algorithm
        %
        % Failure probability can be tuned using replab.Parameters.randomizedSchreierSimsTries
            nTries = replab.Parameters.randomizedSchreierSimsTries;
            R = replab.bsgs.RandomBag(self.n, self.S, [], []);
            c = 0;
            if isempty(order)
                while c <= nTries
                    g = R.sample;
                    if self.stripAndAddStrongGenerator(R.sample)
                        c = 0;
                    else
                        c = c + 1;
                    end
                end
            else
                while self.order < order
                    self.stripAndAddStrongGenerator(R.sample);
                end
            end
        end

    end

    methods (Static)

        function C = make(n, generators, base, order)
            if nargin < 3
                base = [];
            end
            if nargin < 4
                order = [];
            end
            C = replab.bsgs.Chain(n);
            for i = 1:length(base)
                C.insertEndBasePoint(base(i));
            end
            for i = 1:length(generators)
                C.stripAndAddStrongGenerator(generators{i});
            end
            C.randomizedSchreierSims(order);
            C.makeImmutable;
        end

    end

end
