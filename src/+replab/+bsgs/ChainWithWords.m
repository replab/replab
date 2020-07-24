classdef ChainWithWords < replab.Str
% A BSGS chain data structure for a permutation group
%
% This is a variant of `.Chain` which stores

    properties (SetAccess = protected)
        isMutable % (logical): Whether the chain can be modified
        n % (integer): Domain size
        B % (integer(1,nB)): Row vector of base points (all between 1..n without duplicates)

        S % (integer(n,nS)): Matrix of strong generators stored as columns vectors (nS = # of strong generators)
        T % (cell(1,nS) of integer(1,\*)): Strong generators as words
        Sind % (integer(1,nB+1)): Starting index of strong generators for each stabilizer subgroup, where k = length(B)
             %
             %                    We have ``Sind(k+1) == nS + 1``.

        Delta % (cell(1,nB) of integer(1,\*)): Each orbit is a row vector containing orbit elements
        iDelta % (integer(n,nB)): For each orbit, maps a domain element to its position in Delta{i} or 0 if not present

        U % (cell(1,nB) of cell(1,\*) of integer(n,\*)): Transversal elements
        V % (cell(1,nB) of integer(1,\*)): Transversal elements as words
        Uinv % (cell(1,k) of integer(n,\*)): Inverse transversal elements
    end

    methods

        function self = ChainWithWords(n, B, S, T, Sind, Delta, iDelta, U, V, Uinv)
        % Constructs an empty mutable chain for a group of permutations acting on ``n`` elements
        %
        % Args:
        %   n (integer): Domain size
        %
        % Returns:
        %   `+replab.bsgs.Chain`: A constructed BSGS chain
            self.isMutable = true;
            self.n = n;
            if nargin == 8
                self.B = B;
                self.S = S;
                self.T = T;
                self.Sind = Sind;
                self.Delta = Delta;
                self.iDelta = iDelta;
                self.U = U;
                self.V = V;
                self.Uinv = Uinv;
            else
                self.B = zeros(1, 0);
                self.S = zeros(n, 0);
                self.T = cell(1, 0);
                self.Sind = [1];
                self.Delta = cell(1, 0);
                self.iDelta = zeros(n,0);
                self.U = cell(1, 0);
                self.V = cell(1, 0);
                self.Uinv = cell(1, 0);
            end
        end


        %% Immutable functions

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
            else
                o = replab.util.multiplyIntegers(self.orbitSizes);
            end
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

        function img = v(self, i, b)
        % Looks up the word corresponding to the transversal element that maps beta_i to b
        %
        % The element needs to exist.
        %
        % Args:
        %   i (integer): Index of transversal
        %   b (integer): Orbit element
        %
        % Returns:
        %   integer(1,\*): The corresponding transversal element image
        %
        % Raises:
        %   An error if ``b`` is not part of the orbit Delta^i
            j = self.iDelta(b, i);
            assert(j ~= 0, 'Element not part of orbit');
            Vi = self.V{i};
            img = Vi{j};
        end

        function img = vinv(self, i, b)
        % Looks up the word corresponding to the inverse transversal element that maps b to beta_i
        %
        % The element needs to exist.
        %
        % Args:
        %   i: Index of transversal
        %   b: Orbit element
        %
        % Returns:
        %   integer(1,\*): The corresponding image of the inverse transversal element
        %
        % Raises:
        %   An error if ``b`` is not part of the orbit Delta^i
            j = self.iDelta(b, i);
            assert(j ~= 0, 'Element not part of orbit');
            Vi = self.V{i};
            img = -fliplr(Vi{j});
        end

        function img = image(self, g)
        % Returns the image of a chain element
        %
        % Args:
        %   g (permutation row vector): Permutation part of this chain
        %
        % Returns:
        %   The image of the given element ``g``
        %
        % Raises:
        %   An error if the element is not part of the chain
            h = g;
            img = [];
            for i = 1:self.length
                beta_i = self.B(i);
                b = h(beta_i);
                j = self.iDelta(b, i);
                assert(j ~= 0, 'Element is not member of the chain');
                Uinvi = self.Uinv{i};
                uinv = Uinvi(:, j)';
                Vi = self.V{i};
                v = Vi{j};
                h = uinv(h); % compose(uinv, h)
                img = replab.fp.composeLetters(img, v);
            end
        end

        function [h i w] = strip(self, g, v)
        % Strips a permutation through the chain
        %
        % Args:
        %   g (permutation): Permutation
        %   v (integer(1,\*)): Word
        %
        % Returns
        % -------
        %   h: permutation
        %     The part of the group element that could not be sifted
        %   i: integer
        %     Index ``i`` such that ``h(beta_i)`` was not part of the orbit ``Delta^i``
        %   w: integer(1,\*)
        %     The part of the word that could not be sifted
            k = self.length;
            h = g;
            w = v;
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
                Vi = self.V{i};
                vinv = -fliplr(Vi{j});
                w = replab.fp.composeLetters(vinv, w);
            end
            i = k + 1; % marker that we striped through the chain
        end

        %% Mutable methods

        function insertInOrbit(self, i, b, u, v)
        % Inserts a new orbit element in an orbit
        %
        % Modifies the chain in place.
        %
        % Args:
        %   i (integer): Level of the orbit, 1 <= i <= self.length
        %   b (integer): New orbit point
        %   u (permutation): New transversal element
        %   v (integer(1,\*)): Word corresponding to the new transversal element
            assert(self.isMutable);
            assert(self.iDelta(b, i) == 0, 'Orbit point must be new');
            self.Delta{i} = [self.Delta{i} b];
            idx = length(self.Delta{i});
            self.iDelta(b, i) = idx;
            n = self.n;
            uinv = zeros(1, n);
            uinv(u) = 1:n;
            self.U{i} = [self.U{i} u(:)];
            self.Uinv{i} = [self.Uinv{i} uinv(:)];
            self.V{i} = {self.V{i}{:} v};
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
                [b_ind, S_ind] = find(mask);
                for j = 1:length(b_ind)
                    s = self.S(:, Srange(S_ind(j))); % strong generator
                    b = toTest(b_ind(j)); % orbit element
                    newb = s(b);
                    if self.iDelta(newb, i) == 0
                        % new orbit point discovered, add it
                        Ui = self.U{i};
                        ind = self.iDelta(b, i);
                        newu = s(Ui(:, ind)'); % new transversal compose(s, self.u(i, b))
                        newv = replab.fp.composeLetters(self.T{Srange(S_ind(j))}, self.v(i, b));
                        self.insertInOrbit(i, newb, newu, newv);
                    end
                end
                toTest = imgs(mask);
                toTest = unique(toTest(:)');
            end
        end

        function insertEndBasePoint(self, newBeta)
        % Adds a new basis point at the end of the BSGS chain
        %
        % Args:
        %   newBeta (integer): New basis point not part of the current basis
            assert(self.isMutable, 'Chain needs to be mutable');
            assert(all(self.B ~= newBeta), 'Base point already exists');
            n = self.n;
            k = self.length; % previous chain length
            self.B = [self.B newBeta];
            % We add the data structures for the new orbit/transversal
            self.Delta = horzcat(self.Delta, newBeta);
            self.iDelta = [self.iDelta zeros(n, 1)];
            self.iDelta(newBeta, k+1) = 1;
            self.U = {self.U{:} [1:n]'};
            self.Uinv = {self.Uinv{:} [1:n]'};
            self.V = {self.V{:} {[]}};
            self.Sind = [self.Sind self.Sind(end)];
        end

        function addStrongGenerator(self, i, newS, newT)
        % Adds a strong generator at a particular place in the BSGS chain
        %
        % Args:
        %   i (integer): Smallest i such that the strong generator ``newS`` is part of S^(i)
        %   newS (permutation): New strong generator
        %   newT (integer(1,\*)): Word corresponding to the strong generator
            self.S = [self.S(:, 1:(self.Sind(i+1)-1)) newS(:) self.S(:, self.Sind(i+1):end)];
            self.T = {self.T{1:(self.Sind(i+1)-1)} newT self.T{self.Sind(i+1):end}};
            self.Sind((i+1):end) = self.Sind((i+1):end) + 1;
        end

        function addStrongGeneratorAndRecomputeOrbits(self, i, newS, newT)
        % Adds a strong generator at a particular place in the BSGS chain and recomputes orbits
        %
        % Args:
        %   i (integer): Smallest i such that the strong generator ``newS`` is part of S^(i)
        %   newS (permutation): New strong generator
        %   newT (integer(1,\*)): Word corresponding to the strong generator
            self.addStrongGenerator(i, newS, newT);
            for j = 1:i
                self.completeOrbit(j);
            end
        end

        function newSG = stripAndAddStrongGenerator(self, g, v)
        % Strips and, when relevant, adds a new strong generator to this BSGS chain
        %
        % Performs part of the RANDOMSCHREIER procedure in Holt et al., Handbook of CGT, page 98
        %
        % Args:
        %   g (permutation): Element to strip
        %   v (integer(1,\*)): Word
        %
        % Returns:
        %   logical: True if a new strong generator has been found
            n = self.n;
            [h j w] = self.strip(g, v);
            newSG = false;
            if j <= self.length
                % New strong generator h at level j
                newSG = true;
            else
                % Check if h is the identity
                gamma = find(h ~= 1:n, 1);
                if ~isempty(gamma)
                    % New strong generator h fixes all base points
                    % We have a new base point gamma
                    % and insert it at the end of the chain
                    self.insertEndBasePoint(gamma);
                    newSG = true;
                end
            end
            if newSG
                % if we have a new strong generator, add it to the chain
                self.addStrongGeneratorAndRecomputeOrbits(j, h, w);
            end
        end

        function makeImmutable(self)
            assert(self.isMutable);
            self.isMutable = false;
        end

        function inew = schreierSimsTest(self, i)
        % Tests a level of the stabilizer chain for completeness of strong generators
        %
        % Adds the strong generators found to the chain
        %
        % Args:
        %   i (integer): Level to test
        %
        % Returns:
        %   integer: Either $i-1$ if the level is complete, or the new level to test
            orbit = self.Delta{i};
            iOrbit = self.iDelta;
            U = self.U{i};
            Uinv = self.Uinv{i};
            V = self.V{i};
            Srange = self.Sind(i):self.Sind(end)-1;
            n = self.n;
            for o = 1:length(orbit)
                betai = self.B(i);
                b = orbit(o);
                ub = U(:,o)';
                vb = V{o};
                for k = Srange
                    x = self.S(:,k)';
                    y = self.T{k};
                    uxb_inv = Uinv(:,iOrbit(x(b),i))';
                    vxb_inv = -fliplr(V{iOrbit(x(b),i)});
                    toStrip = ub(x(uxb_inv));
                    toStripW = replab.fp.composeLetters(vb, replab.fp.composeLetters(y, vxb_inv));
                    if any(toStrip ~= 1:n)
                        flag = true;
                        [h j w] = self.strip(toStrip, toStripW);
                        if j <= self.length
                            % new strong generator h at level j
                            flag = false;
                        else % j = self.length + 1
                            gamma = find(h ~= 1:n, 1);
                            if ~isempty(gamma)
                                % h fixes
                                flag = false;
                                self.insertEndBasePoint(gamma);
                            end
                        end
                        if ~flag
                            self.addStrongGeneratorAndRecomputeOrbits(j, h, w);
                            inew = j;
                            return
                        end
                    end
                end
            end
            inew = i - 1;
        end

        function deterministicSchreierSims(self)
        % Runs the deterministic Schreier-Sims algorithm, with a bound on the maximum order
        %
        % If ``self.order > maxOrder`` after this function returns, the chain may be incomplete.
        %
        % Args:
        %   maxOrder (integer or vpi or ``inf``): Order cutoff
            i = self.length;
            while i >= 1
                i = self.schreierSimsTest(i);
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
            C = replab.bsgs.ChainWithWords(n);
            for i = 1:length(base)
                C.insertEndBasePoint(base(i));
            end
            for i = 1:length(generators)
                C.stripAndAddStrongGenerator(generators{i}, i);
            end
            C.deterministicSchreierSims;
            C.makeImmutable;
        end

    end

end
