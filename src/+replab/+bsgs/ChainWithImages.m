classdef ChainWithImages < replab.Str
% A BSGS chain data structure for a permutation group
%
% This class uses essentially the same data structures as `+replab.+bsgs.Chain`, but it also stores
% the images of a group homomorphism, with target group ``J``.
%
% For that, images of the strong generators are stored in ``T``, and the images of transversal
% elements are stored in row cell arrays (containing row cell arrays of group elements)
% ``V`` and ``Vinv``, with conventions similar as ``U`` and ``Uinv``.

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

        J % image group
        T % row cell array of images of strong generators
        V % row cell array of row cell arrays of transversal images
        Vinv % row cell array of row cell arrays of inverse transversal images
    end

    methods

        function self = ChainWithImages(n, J, B, S, T, Sind, Delta, iDelta, U, Uinv, V, Vinv)
        % Constructs an empty mutable chain for a group of permutations acting on ``n`` elements
        %
        % Args:
        %   n (integer): Domain size
        %
        % Returns:
        %   `+replab.bsgs.ChainWithImages`: A constructed BSGS chain
            self.isMutable = true;
            self.n = n;
            self.J = J;
            if nargin > 2
                self.B = B;
                self.S = S;
                self.T = T;
                self.Sind = Sind;
                self.Delta = Delta;
                self.iDelta = iDelta;
                self.U = U;
                self.Uinv = Uinv;
                self.V = V;
                self.Vinv = Vinv;
            else
                self.B = [];
                self.S = [];
                self.Sind = [1];
                self.T = {};
                self.Delta = {};
                self.iDelta = zeros(n,0);
                self.U = {};
                self.Uinv = {};
                self.J = J;
                self.T = {};
                self.V = {};
                self.Vinv = {};
            end
        end

        function c = toChain(self)
        % Strips the image information and returns the BSGS chain
        %
        % Returns:
        %   `+replab.+bsgs.Chain`: BSGS chain with images stripped
            c = replab.bsgs.Chain(self.n, self.B, self.S, self.Sind, self.Delta, self.iDelta, self.U, self.Uinv)
        end

        function show(self, i)
            if nargin < 2
                self.toChain.show;
            else
                fprintf('Level %d Base point %d\n', i, self.B(i));
                table = cell(length(self.Delta{i})+1, 5);
                table{1,1} = 'Orbit pt ';
                table{1,2} = 'Transversal ';
                table{1,3} = 'Trs. inv. ';
                table{1,4} = 'Trs image';
                table{1,5} = 'Trs inv image';
                for j = 1:length(self.Delta{i})
                    table{j+1,1} = sprintf('%d', self.Delta{i}(j));
                    table{j+1,2} = strrep(replab.shortStr(self.U{i}(:,j)'), ' ', '');
                    table{j+1,3} = strrep(replab.shortStr(self.Uinv{i}(:,j)'), ' ', '');
                    table{j+1,4} = replab.shortStr(self.T{i}{j});
                    table{j+1,5} = replab.shortStr(self.Tinv{i}{j});
                end
                disp(strjoin(replab.str.alignspace(table, 'cllll'), '\n'));
            end
        end


        %% Immutable functions

        function s = orbitSizes(self)
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
            o = vpi(1);
            for i = 1:self.length
                o = o * vpi(length(self.Delta{i}));
            end
        end

        function [g v] = sample(self)
        % Samples an element uniformly from the group along with its image
        %
        % Returns
        % -------
        %  g: permutation
        %    Random group element
        %  v: element of `J`
        %    Image of ``g``
            g = 1:self.n;
            v = self.J.identity;
            for i = 1:self.length
                [gi vi] = self.randomTransversal(i);
                g = g(gi); % compose(g, gi)
                v = self.J.compose(v, vi);
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

        function img = strongGeneratorImage(self, i)
        % Returns the image of the i-th strong generator under the stored homomorphism
        %
        % Args:
        %  i (integer): Strong generator index
        %
        % Returns:
        %   A group element of `J`
            p = self.T{i};
        end

        function [u v]  = randomTransversal(self, i)
        % Returns a random transversal element from a specific transversal set along with its image
        %
        % Args:
        %  i (integer): Index of the transversal set
        %
        % Returns
        % -------
        %  u: permutation
        %    Random transversal element
        %  v: element of `J`
        %    Image of ``g``
            j = randi(length(self.Delta{i}));
            Ui = self.U{i};
            u = Ui(:,j)';
            Vi = self.V{i};
            v = Vi{j};
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

        function img = v(self, i, b)
        % Looks up the image of the transversal element that maps beta_i to b
        %
        % The element needs to exist.
        %
        % Args:
        %   i (integer): Index of transversal
        %   b (integer): Orbit element
        %
        % Returns:
        %   The corresponding transversal element image
        %
        % Raises:
        %   An error if ``b`` is not part of the orbit Delta^i
            j = self.iDelta(b, i);
            assert(j ~= 0, 'Element not part of orbit');
            Vi = self.V{i};
            img = Vi{j};
        end

        function img = vinv(self, i, b)
        % Looks up the image of the inverse transversal element that maps b to beta_i
        %
        % The element needs to exist.
        %
        % Args:
        %   i: Index of transversal
        %   b: Orbit element
        %
        % Returns:
        %   The corresponding image of the inverse transversal element
        %
        % Raises:
        %   An error if ``b`` is not part of the orbit Delta^i
            j = self.iDelta(b, i);
            assert(j ~= 0, 'Element not part of orbit');
            Vinvi = self.Vinv{i};
            img = Vinvi{j};
        end

        function img = image(self, g)
        % Returns the image of a chain element
        %
        % Args:
        %   g (permutation): Permutation to compute the image of
        %
        % Returns:
        %   The image of the given element ``g``
        %
        % Raises:
        %   An error if the element is not part of the chain
            h = g;
            img = self.J.identity;
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
                img = self.J.compose(img, v);
            end
        end

        function img = inverseImage(self, g)
        % Returns the inverse image of a chain element
        %
        % Satisfies
        %
        % self.J.compose(self.image(g), self.inverseImage(g)) == self.J.identity
        %
        % Args:
        %   g (permutation): Permutation part of this chain
        %
        % Returns:
        %   The inverse image of the given element ``g``
        %
        % Raises:
        %   An error if the element is not part of the chain
            h = g;
            img = self.J.identity;
            for i = 1:self.length
                beta_i = self.B(i);
                b = h(beta_i);
                j = self.iDelta(b, i);
                assert(j ~= 0, 'Element is not member of the chain');
                Uinvi = self.Uinv{i};
                uinv = Uinvi(:, j)';
                Vinv = self.Vinv{i};
                vinv = Vinv{j};
                h = uinv(h); % compose(uinv, h)
                img = self.J.compose(vinv, img);
            end
        end

        function b = contains(self, g)
        % Tests whether the BSGS chain contains an element
            k = self.length;
            n = self.n;
            [h, i] = self.strip(g);
            b = (i > k) && all(h == 1:n); % i > k and h is the identity
        end


        function T = imagesDecomposition(self)
        % Returns the homomorphism group images into a product of transversals
        %
        % Guarantees that the first element of each transversal is the identity
            k = self.length;
            T = self.V;
        end

        function [h i w] = strip(self, g, v)
        % Strips a permutation through the chain along with its image
        %
        % Args:
        %   g (row permutation vector): A permutation group element
        %   v (element of `J`): The image of ``g`` under the encoded homomorphism
        %
        % Returns
        % -------
        %   h: row permutation vector
        %     The part of the group element that could not be sifted
        %   i: integer
        %     Index i-th such that h(beta_i) was not part of the orbit Delta^i
        %   w: element of `J`
        %     The part of the image that could not be sifted
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
                Vinvi = self.Vinv{i};
                vinv = Vinvi{j};
                w = self.J.compose(vinv, w);
            end
            i = k + 1; % marker that we striped through the chain
        end

        function [res, errorBound] = double(self, knownUnitary)
        % Returns a new BSGS chain with the images approximated
        %
        % Args:
        %   knownUnitary (logical): Whether the representation described by the chain is known to be unitary
        %
        % Returns
        % -------
        %   res: `+replab.+bsgs.ChainWithImages`
        %     Chain with approximate images
        %   errorBound: double
        %     Bound on the computed images
            if nargin < 2
                knownUnitary = false;
            end
            k = self.length;
            newT = cellfun(@(t) double(t), self.T, 'uniform', 0);
            newV = cell(1, k);
            newVinv = cell(1, k);
            maxCondNums = zeros(1, k);
            maxErrors = zeros(1, k);
            for i = 1:k
                Vi = self.V{i};
                Vinvi = self.Vinv{i};
                l = length(Vi);
                newVi = cell(1, l);
                newVinvi = cell(1, l);
                maxError = 0;
                maxCondNum = 1;
                for j = 1:length(Vi)
                    if isa(Vi{j}, 'replab.cyclotomic')
                        [approx, err] = Vi{j}.doubleApproximation;
                    else
                        approx = Vi{j};
                        err = 0;
                    end
                    newVi{j} = replab.numerical.bestStorage(approx);
                    maxError = max(maxError, norm(err, 'fro'));
                    if knownUnitary
                        newVinv{j} = newVi{j}';
                    else
                        if isa(Vinvi{j}, 'replab.cyclotomic')
                            [approx, err] = Vinvi{j}.doubleApproximation;
                        else
                            approx = Vinvi{j};
                            err = 0;
                        end
                        newVinvi{j} = replab.numerical.bestStorage(approx);
                        maxCondNum = max(maxCondNum, replab.numerical.condUpperBound(newVi{j}, newVinvi{j}));
                        maxError = max(maxError, norm(err, 'fro'));
                    end
                end
                maxCondNums(i) = maxCondNum;
                maxErrors(i) = maxError;
                newV{i} = newVi;
                newVinv{i} = newVinvi;
            end
            errorBound = 0;
            for i = 1:k
                errorBound = errorBound + maxErrors(i)*prod(maxCondNums(1:i-1))*prod(maxCondNums(i+1:end));
            end
            res = replab.bsgs.ChainWithImages(self.n, self.J, self.B, self.S, newT, self.Sind, self.Delta, self.iDelta, ...
                                              self.U, self.Uinv, newV, newVinv);
            if ~self.isMutable
                res.makeImmutable;
            end
        end

        function res = mapImages(self, mu)
        % Returns a new BSGS chain with the images mapped through a function
        %
        % The returned chain has the same mutability (`.isMutable`) as this chain.
        %
        % Args:
        %   mu (`+replab.Morphism`): Group homomorphism from the current `J` to ``newJ``
        %
        % Returns:
        %   `+replab.+bsgs.ChainWithImages`: A copy of this BSGS chain with updated images
            k = self.length;
            f = @(v) mu.imageElement(v);
            newT = cellfun(f, self.T, 'uniform', 0);
            newV = cell(1, k);
            newVinv = cell(1, k);
            for i = 1:k
                newV{i} = cellfun(f, self.V{i}, 'uniform', 0);
                newVinv{i} = cellfun(f, self.Vinv{i}, 'uniform', 0);
            end
            newJ = mu.target;
            res = replab.bsgs.ChainWithImages(self.n, newJ, self.B, self.S, newT, self.Sind, self.Delta, self.iDelta, ...
                                              self.U, self.Uinv, newV, newVinv);
            if ~self.isMutable
                res.makeImmutable;
            end
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
        %   v (element of `J`): Image of the new transversal element
            assert(self.isMutable);
            assert(self.iDelta(b, i) == 0, 'Orbit point must be new');
            self.Delta{i} = [self.Delta{i} b];
            idx = length(self.Delta{i});
            self.iDelta(b, i) = idx;
            n = self.n;
            uinv = zeros(1, n);
            uinv(u) = 1:n;
            vinv = self.J.inverse(v);
            self.U{i} = [self.U{i} u(:)];
            self.Uinv{i} = [self.Uinv{i} uinv(:)];
            self.V{i} = {self.V{i}{:} v};
            self.Vinv{i} = {self.Vinv{i}{:} vinv};
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
                        newv = self.J.compose(self.T{Srange(S_ind(j))}, self.v(i, b));
                        self.insertInOrbit(i, newb, newu, newv);
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
            self.V = {self.V{:} {self.J.identity}};
            self.Vinv = {self.Vinv{:} {self.J.identity}};
            self.Sind = [self.Sind self.Sind(end)];
        end

        function addStrongGenerator(self, i, newS, newT)
        % Adds a strong generator at a particular place in the BSGS chain
        %
        % Args:
        %   i (integer): Smallest i such that the strong generator ``newS`` is part of S^(i)
        %   newS (permutation): New strong generator
        %   newT (element of `J`): Image of the new strong generator
            self.S = [self.S(:, 1:(self.Sind(i+1)-1)) newS(:) self.S(:, self.Sind(i+1):end)];
            self.T = {self.T{1:(self.Sind(i+1)-1)} newT self.T{self.Sind(i+1):end}};
            self.Sind((i+1):end) = self.Sind((i+1):end) + 1;
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
        %   g (row permutation vector): Element to strip
        %   v (element of `J`, optional): Image of the element
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
                self.addStrongGenerator(j, h, w);
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
            R = replab.bsgs.RandomBagWithImages(self.n, self.S, [], [], self.J, self.T);
            c = 0;
            if isempty(order)
                while c <= nTries
                    [g j] = R.sample;
                    if self.stripAndAddStrongGenerator(g, j)
                        c = 0;
                    else
                        c = c + 1;
                    end
                end
            else
                while self.order < order
                    [g j] = R.sample;
                    self.stripAndAddStrongGenerator(g, j);
                end
            end
        end

    end

    methods (Static)

        function C = make(n, target, preimages, images, base, order)
            if nargin < 6
                order = [];
            end
            if nargin < 5
                base = [];
            end
            C = replab.bsgs.ChainWithImages(n, target);
            for i = 1:length(base)
                C.insertEndBasePoint(base(i));
            end
            for i = 1:length(preimages)
                C.stripAndAddStrongGenerator(preimages{i}, images{i});
            end
            C.randomizedSchreierSims(order);
            C.makeImmutable;
        end

    end

end
