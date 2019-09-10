classdef Chain < replab.Str
% A BSGS chain data structure for a permutation group
%
% The represented group acts on {1, ..., n}. We follow loosely the
% notation in Handbook of CGT, Derek Holt.
%
% The BSGS chain is stored as follows.
%
% We write B = [beta_1, beta_2 ... beta_k] where k is the length of the stabilizer chain.
% We define the base order on 1...n as follows: the smallest elements in the order are given
% by the base elements, such that beta_1 < beta_2 < ... < beta_k, followed by the domain elements
% not in the base, sorted in (the standard) increasing order.
% We write G^i = {g \in G : g(beta_1) = beta_1, ..., g(beta_(i-1)) = beta_(i-1)}, so that
% G^1 = G.
% The (currently known) strong generators are stored as column vectors in a matrix S.
%
% We write the set S^i = S intersect G^i, and the subgroups H^i = <S^i>.
%
% The order of strong generators as columns in the matrix is such that
% S^i = S(:, Sind(i):end) where Sind is a double row vector of starting indices in S.
%
% We define Delta^i = H^i(beta_i) = { h(beta_i) : h \in H^i }, and those orbits are stored
% as a cell array Delta = { Delta^1, Delta^2, ... }; each Delta^i is a double row vector
% sorted according to the base order.
%
% We define the transversal elements u^i_b such that u^i_b(beta_i) = b, and the inverse transversal
% elements uInv^i_(beta_i) = b.
%
% The transversal elements for a particular node in the chain are stored as column vectors in a matrix,
% with an order corresponding to the orbit Delta^i; those matrices are stored in the row cell array U.
%
% Transversal element inverses are stored as column vectors in a matrix, in a row cell array Uinv.
%
% Optionally, this class can store the images of a group homomorphism. 
    
% For that, images of the strong generators are stored in J, and the images of transversal
% elements are stored in row cell arrays (containing row cell arrays of group elements) V and Vinv, with
% conventions similar as U and Uinv.
%
% The following invariants are maintained by the code below.
%
% When the chain is immutable: the strong generating set must be a strong generating set for the G^i, not
% the H^i, i.e. the BSGS construction is complete.
%
% When the chain is mutable:
%
% * The chain describes the orbits/transversals of H^i.
%
% * The base B, base order bo, boInv, the data structures for orbits and transversals Delta, U, Uinv are consistent.
%
% * The strong generating set is ordered, so that the starting indices Sind describe the sequence S^i
%
% * (When withImages is true) The above is consistent with the image data structures T, V, Vinv

    properties (SetAccess = protected)
        isMutable % whether the chain can be modified
        n % domain size
        G % Group of permutation on 1..n
        B % row vector of base points (all between 1..n without duplicates)   
        bo % base order, i.e. for each domain element, gives its order
        boinv % inverse permutation of bo, permutation that starts with B and then the remaining domain elements increasing
        S % n x nS matrix of strong generators stored as columns vectors (nS = # of strong generators)
        Sind % starting index of strong generators for each stabilizer subgroup, of length k+1 if k = length(B)
             % S^(k+1) corresponds to strong generators that are stabilized by all base points
             % if no strong generators are present, use the default starting index 1
        Delta % row cell array of orbits, each orbit is a row double vector with orbit elements sorted according to the base order
        U % row cell array of transversal elements stored as n x orbitSize matrices
        Uinv % row cell array of inverse transversal elements
        withImages % whether we compute images of a morphism
        J % image group
        T % row cell array of images of strong generators
        V % row cell array of row cell arrays of transversal images
        Vinv % row cell array of row cell arrays of inverse transversal images
    end
    
    methods
        
        function self = Chain(n, J)
        % Constructs an empty mutable chain for a subgroup of the permutation group
        %
        % The subgroup acts on 1..n
        %
        % Args:
        %   n: Domain size
        %   J (replab.Group, optional): Group structure for morphism images
        %
        % Returns:
        %   A constructed empty BSGS chain
            
            if nargin < 2
                J = [];
            end
            self.isMutable = true;
            self.n = n;
            self.G = replab.Permutations(n);
            self.B = [];
            self.bo = 1:n;
            self.boinv = 1:n;
            self.S = [];
            self.Sind = [1];
            self.Delta = {};
            self.U = {};
            self.Uinv = {};
            if isequal(J, [])
                self.withImages = false;
            else
                self.withImages = true;
                self.J = J;
                self.T = {};
                self.V = {};
                self.Vinv = {};
            end
        end
        
        function check(self)
            self.checkTransversals;
            self.checkStrongGenerators;
        end

        function checkTransversals(self)
            n = self.n;
            k = self.length;
            for i = 1:k
                members = false(1, n);
                beta_i = self.B(i);
                orbit = self.Delta{i};
                members(orbit) = true;
                for j = 1:length(orbit)
                    b = orbit(j);
                    Ui = self.U{i};
                    ub1 = Ui(:,j);
                    ub2 = self.u(i, b);
                    assert(isequal(ub1, ub2), 'inconsistent transversal retrieval');
                    assert(ub1(beta_i) == b, 'inconsistent transversal element');
                    Uinvi = self.Uinv{i};
                    ubinv1 = Uinvi(:,j);
                    ubinv2 = self.uinv(i, b);
                    assert(isequal(ubinv1, ubinv2), 'inconsistent transversal retrieval');
                    assert(ubinv1(b) == beta_i, 'inconsistent transversal element');
                    for l = self.Sind(i):self.nStrongGenerators
                        imgb = self.S(b, l);
                        assert(members(imgb));
                    end
                end
            end
        end
        
        function checkStrongGenerators(self)
        % Performs checks on the strong generators of the current chain
            k = self.length;
            nS = self.nStrongGenerators;
            Sind = [self.Sind nS+1];
            for i = 1:k
                betai = self.B(i);
                for j = Sind(i):Sind(i+1)-1
                    % strong generators particular to this step move beta_i
                    assert(S(betai, j) ~= betai);
                end
                for j = Sind(i+1):nS
                    % strong generators of stabilizer subgroups are stabilized
                    assert(S(betai, j) == betai);
                end
            end
        end
        
        %% Immutable functions
        function k = length(self)
            k = length(self.B);
        end
        
        function nS = nStrongGenerators(self)
        % Returns the number of strong generators in this BSGS chain
            nS = size(self.S, 2);
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
        %   A group element of `self.J`
        end
        
        function g = u(self, i, b)
        % Looks up a transversal element that maps beta_i to b
        %
        % Args:
        %   i: Index of transversal
        %   b: Orbit element
        %
        % Returns:
        %   The corresponding transversal element, or [] if b is not part of the orbit
            [~, j] = ismember(b, self.Delta{i});
            if j == 0
                g = [];
                return
            end
            Ui = self.U{i};
            g = Ui(:,j)';
        end
        
        function img = v(self, i, b)
        % Looks up the image of the transversal element that maps beta_i to b
        %
        % Args:
        %   i: Index of transversal
        %   b: Orbit element
        %
        % Returns:
        %   The corresponding transversal element image if it exists
        %
        % Raises:
        %   An error if `self.withImages` is false
        %   An error if `b` is not part of the orbit Delta^i
            [~, j] = ismember(b, self.Delta{i});
            assert(j ~= 0);
            Vi = self.V{i};
            img = Vi{j};
        end
                    
        function g = uinv(self, i, b)
        % Looks up the inverse transversal element that maps b to beta_i
        % Args:
        %   i: Index of transversal
        %   b: Orbit element
        %
        % Returns:
        %   The corresponding inverse transversal element or [] if b is not part
        %   of the orbit
            [~, j] = ismember(b, self.Delta{i});
            if j == 0
                g = [];
                return
            end
            Uinvi = self.Uinv{i};
            g = Uinvi(:,j)';
        end
        
        function img = vinv(self, i, b)
        % Looks up the image of the inverse transversal element that maps b to beta_i
        %
        % Args:
        %   i: Index of transversal
        %   b: Orbit element
        %
        % Returns:
        %   The corresponding image of the inverse transversal element
        %
        % Raises:
        %   An error if `self.withImages` is false
        %   An error if `b` is not part of the orbit Delta^i
            [~, j] = ismember(b, self.Delta{i});
            assert(length(j) == 1, 'Element not part of orbit');
            Vinvi = self.Vinv{i};
            img = Vinvi{j};
        end
        
        function img = image(self, g)
        % Returns the image of a chain element
        %
        % Args:
        %   g (permutation row vector): Permutation part of this chain
        %
        % Returns:
        %   The image of the given element `g`
        %
        % Raises:
        %   An error if the element is not part of the chain
            assert(self.withImages);
            h = g;
            img = self.J.identity;
            for i = 1:self.length
                beta_i = self.B(i);
                b = h(beta_i);
                [~, j] = ismember(b, self.Delta{i});
                assert(j ~= 0, 'Element is not member of the chain');
                Uinvi = self.Uinv{i};
                uinv = Uinvi(:, j)';
                Vi = self.V{i};
                v = Vi{j};
                h = self.G.compose(uinv, h);
                img = self.J.compose(img, v);
            end
        end

        function [h i w] = strip(self, g, v)
        % Strips a permutation through the chain
        %
        % Args:
        %   g (row permutation vector): A permutation group element
        %   v (element of `self.J`, optional): The image of `g` under the encoded homomorphism
        %
        % Returns
        % -------
        %   h: row permutation vector
        %     The part of the group element that could not be sifted
        %   i: integer
        %     Index i-th such that h(beta_i) was not part of the orbit Delta^i
        %   w: element of `self.J`, optional
        %     The part of the image that could not be sifted
            k = self.length;
            h = g;
            if nargin > 2
                assert(self.withImages);
                img = true;
                w = v;
            else
                img = false;
            end
            for i = 1:k
                b = h(self.B(i));
                [~, j] = ismember(b, self.Delta{i});
                if j == 0
                    return
                end
                Uinvi = self.Uinv{i};
                uinv = Uinvi(:, j)';                
                % note order is reversed compared to Holt, as
                % we use a left action
                h = self.G.compose(uinv, h);
                if img
                    Vinvi = self.Vinv{i};
                    vinv = Vinvi{j};
                    w = self.J.compose(vinv, w);
                end
            end
            i = k + 1; % marker that we striped through the chain
        end
        
        %% Mutable methods
        
        function mutableMapImages(self, newJ, f)
        % Maps in place the homomorphism images
        %
        % Args:
        %   newJ: New image parent group
        %   f: Group homomorphism from the current `self.J` to `newJ`
            k = self.length;
            self.T = cellfun(f, self.T, 'uniform', 0);
            for i = 1:k
                Vi = self.V{i};
                self.V{i} = cellfun(f, Vi, 'uniform', 0);
                Vinvi = self.Vinv{i};
                self.Vinv{i} = cellfun(f, Vinvi, 'uniform', 0);
            end
            self.J = newJ;
        end
        
        function recomputeBaseOrder(self)
        % Recomputes the base order, and reorders all orbits according to the corrected base order
            assert(self.isMutable);
            remaining = setdiff(1:self.n, self.B); % domain elements not part of base
            self.boInv = [self.B remaining];
            self.bo = self.G.inverse(self.boInv);
            for i = 1:self.length
                self.reorderOrbit(i);
            end
        end

        function reorderOrbit(self, i)
        % Reorders the i-th orbit elements according to the base order
            assert(self.isMutable);
            D = self.Delta{i};
            Ui = self.U{i};
            Uinvi = self.Uinv{i};
            [~, I] = sort(self.bo(D));
            % sort index vector in I
            self.Delta{i} = D(I);
            self.U{i} = Ui(:,I);
            self.Uinv{i} = Uinvi(:,I);
            if self.withImages
                Vi = self.V{i};
                Vinvi = self.Vinv{i};
                self.V{i} = Vi(I);
                self.Vinv{i} = Vinvi(I);
            end
        end
        
        function insertInOrbit(self, i, b, iS)
        % Inserts a new orbit element in an orbit
        %
        % The new orbit point is given by g(b), where g is the iS-th strong generator.
        %
        % Modifies the chain in place.
        %
        % Args:
        %   i (double): Index of the orbit, 1 <= i <= self.length
        %   b (double): Existing orbit point
        %   iS (double): Index of strong generator
            D = self.Delta{i};
            pos = 1;
            g = self.S(:,iS)';
            newb = g(b);
            newu = self.G.compose(g, self.u(i, b));
            newuinv = self.G.inverse(newu);
            if self.withImages
                newv = self.J.compose(self.T{iS}, self.v(i, b));
                newvinv = self.J.inverse(newv);
            end
            % Look for the position `pos` where to insert the new orbit point
            while pos <= length(D) && self.bo(D(pos)) < self.bo(newb)
                pos = pos + 1;
            end
            self.Delta{i} = [D(1:pos-1) newb D(pos:end)];
            Ui = self.U{i};
            self.U{i} = [Ui(:,1:pos-1) newu' Ui(:,pos:end)];
            Uinvi = self.Uinv{i};
            self.Uinv{i} = [Uinvi(:,1:pos-1) newuinv' Uinvi(:,pos:end)];
            if self.withImages
                Vi = self.V{i};
                self.V{i} = {Vi{1:pos-1} newv Vi{pos:end}};
                Vinvi = self.Vinv{i};
                self.Vinv{i} = {Vinvi{1:pos-1} newvinv Vinvi{pos:end}};
            end
        end
        
        function completeOrbit(self, i)
        % Completes the i-th orbit
        %
        % Iterates over current orbit elements and strong generators, and add new orbit
        % points if necessary.
            assert(self.isMutable);
            n = self.n;
            touched = false(1, n); % elements that were already considered
            toTest = self.Delta{i}; % we need to test all elements in the currently known orbit
            touched(toTest) = true; % and note that they were already considered
            nS = self.nStrongGenerators;
            while length(toTest) > 0
                b = toTest(end);
                toTest = toTest(1:end-1);
                for j = self.Sind(i):nS
                    imgb = self.S(b, j);
                    if ~touched(imgb)
                        touched(imgb) = true;
                        toTest = [toTest imgb];
                        % new orbit point discovered, add it
                        self.insertInOrbit(i, b, j);
                    end
                end
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
            nS = self.nStrongGenerators;
            % We add the data structures for the new orbit/transversal
            self.Delta = horzcat(self.Delta, newBeta);
            self.U = horzcat(self.U, {[1:n]'});
            self.Uinv = horzcat(self.Uinv, {[1:n]'});
            % We check the strong generators that were stabilizing all previous base points
            % They may move the new base point
            I = self.Sind;
            unchanged = self.S(:, 1:I(k+1)-1);
            candidates = self.S(:, I(k+1):nS);
            if size(candidates, 2) ~= 0
                stabilized = [candidates(newBeta, :) == newBeta];
                % Sort the strong generators
                self.S = [unchanged candidates(:, ~stabilized) candidates(:, stabilized)];
            else
                stabilized = [];
            end
            % New starting indices
            self.Sind = [self.Sind(1:k+1) self.Sind(k+1)+sum(~stabilized)];
            if self.withImages
                % Take care of homomorphism images if necessary
                self.V = horzcat(self.V, {{self.J.identity}});
                self.Vinv = horzcat(self.Vinv, {{self.J.identity}});
                if size(candidates, 2) ~= 0
                    Tunchanged = self.T(1:I(k+1)-1);
                    Tcandidates = self.T(I(k+1):nS);
                    self.T = horzcat(Tunchanged, Tcandidates(~stabilized), Tcandidates(stabilized));
                end
            end
            self.completeOrbit(k+1);
        end
        
        function addStrongGenerator(self, i, newS, newT)
        % Adds a strong generator at a particular place in the BSGS chain
        %
        % Args:
        %   i (integer): Smallest i such that the strong generator `newS` is part of S^(i)
        %   newS (row permutation vector): New strong generator
        %   newT (element of `self.J`, optional): Image of the new strong generator 
            I = self.Sind;
            self.S = [self.S(:, 1:(I(i+1)-1)) newS' self.S(:, I(i+1):end)];
            self.Sind((i+1):end) = self.Sind((i+1):end) + 1;
            if self.withImages
                assert(nargin == 4, 'When encoding homomorphism images, image of strong generator must be provided.');
                self.T = {self.T{1:(I(i+1)-1)} newT self.T{I(i+1):end}};
            end
            for j = 1:i
                % TODO: optimization opportunity by starting the orbit completion at 2 if all strong generators
                % have already been added
                %
                % Also: need only to check completeness with respect to the newly added strong generator
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
        %   v (element of `self.J`, optional): Image of the element
        %
        % Returns:
        %   true if a new strong generator has been found
            n = self.n;
            if self.withImages
                [h j w] = self.strip(g, v);
            else
                [h j] = self.strip(g);
            end
            newSG = false; 
            if j <= self.length
                % New strong generator h at level j
                newSG = true;
            elseif ~self.G.isIdentity(h)
                % New strong generator h fixes all base points
                gamma = find(h ~= 1:n, 1);
                % We find a new base point gamma
                assert(length(gamma) == 1, 'Cannot be identity');
                % and insert it at the end of the chain
                self.insertEndBasePoint(gamma);
                newSG = true;
            end
            if newSG
                % if we have a new strong generator, add it to the chain
                if self.withImages
                    self.addStrongGenerator(j, h, w);
                else
                    self.addStrongGenerator(j, h);
                end
            end
        end
        
        function insertStrongGenerators(self, newS, newT)
        % Inserts the given strong generators
        %
        % Args:
        %   newS (integer matrix): Strong generators given in a n x nGens matrix
        %   newT (row cell vector, optional): Images of those strong generators
            n = self.n;
            for i = 1:size(newS, 2)
                h = newS(:,i)'; % candidate
                % find the place j where to insert the current strong generator
                j = 1;
                while j <= self.length && h(self.B(j)) == self.B(j)
                    j = j + 1;
                end
                if j > self.length
                    % the new strong generator fixes all base points,
                    % so we insert a new base point, see stripAndAddStrongGenerator
                    gamma = find(h ~= 1:n, 1);
                    assert(length(gamma) == 1, 'Cannot be identity');
                    self.insertEndBasePoint(gamma);
                end
                if self.withImages
                    self.addStrongGenerator(j, h, newT{i});
                else
                    self.addStrongGenerator(j, h);
                end
            end
        end
        
        function makeImmutable(self)
            assert(self.isMutable);
            self.isMutable = false;
        end
        
        function randomizedSchreierSims(self, nTries)
        % Runs the randomized Schreier-Sims algorithm
        %
        % Args:
        %   nTries: Number of attempts of successive failed attempts to generate a new strong generator
        %           before deciding the chain is complete; the probability of failure is then less than
        %           1/2^nTries
            if nargin < 2
                nTries = 1000;
            end
            if self.withImages
                R = replab.bsgs1.RandomBag(self.n, self.S, [], [], self.J, self.T);
                c = 0;
                while c <= nTries
                    [g v] = R.sample;
                    if self.stripAndAddStrongGenerator(g, v)
                        c = 0;
                    else
                        c = c + 1;
                    end
                end
            else
                R = replab.bsgs1.RandomBag(self.n, self.S);
                c = 0;
                while c <= nTries
                    g = R.sample;
                    if self.stripAndAddStrongGenerator(g)
                        c = 0;
                    else
                        c = c + 1;
                    end
                end
            end
        end

    end
    
    methods (Static)
       
        function C = makeWithImages(n, S, J, T)
            C = replab.bsgs1.Chain(n, J);
            C.insertStrongGenerators(S, T);
            C.randomizedSchreierSims;
            C.makeImmutable;
        end
        
        function C = make(n, S)
            C = replab.bsgs1.Chain(n);
            C.insertStrongGenerators(S);
            C.randomizedSchreierSims;
            C.makeImmutable;
        end
        
    end

end
