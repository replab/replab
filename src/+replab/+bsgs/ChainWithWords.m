classdef ChainWithWords < replab.Str
% A BSGS chain data structure for a permutation group
%
% This is a variant of `.Chain` which stores the words corresponding to the transversal elements.
%
% The notation below is adapted to the one in
% T. Minkwitz, "An Algorithm for Solving the Factorization Problem in Permutation Groups,", Journal of Symbolic Computation, vol. 26, no. 1, pp. 89–95, Jul. 1998
% doi: 10.1006/jsco.1998.0202


    properties (SetAccess = protected)
        group % (`+replab.PermutationGroup`): Permutation group for which this chain is computed
        chain % (`.Chain`): BSGS chain of the group being computed
        completed % (logical): Whether the chain has been completed
        n % (integer): Domain size
        k % (integer): Chain length
        B % (integer(1,k)): Row vector of base points (all between 1..n without duplicates)
        orbit % (cell(1,k) of integer(1,\*)): Each orbit is a row vector containing orbit elements
        iOrbit % (integer(n,k)): For each orbit, maps a domain element to its position in orbit{i} or 0 if not present
        newOrbit % (logical(n,k)): For each orbit, whether the element is new
        nu % (cell(1,k) of integer(n,\*)): Inverse transversal elements
        nuw % (cell(1,k) of cell(1,\*) of integer(1,\*)): Inverse transversal words
    end

    methods

        function self = ChainWithWords(group)
        % Constructs an empty mutable chain to start the Minkwitz algorithm
        %
        % Args:
        %   group (`+replab.PermutationGroup`): Group for which we compute a factorization of the elements
            n = group.domainSize;
            chain = group.chain;
            B = chain.B;
            k = length(B);
            orbit = cell(1, k);
            iOrbit = zeros(n, k);
            newOrbit = false(n, k); % no need to say that the identity is new
            nu = cell(1, k);
            nuw = cell(1, k);
            for i = 1:k
                beta = B(i);
                orbit{i} = beta;
                iOrbit(beta,i) = 1;
                nu{i} = (1:n)';
                nuw{i} = {[]};
            end
            self.group = group;
            self.chain = chain;
            self.completed = false;
            self.n = n;
            self.k = k;
            self.B = B;
            self.orbit = orbit;
            self.iOrbit = iOrbit;
            self.newOrbit = newOrbit;
            self.nu = nu;
            self.nuw = nuw;
        end

        function setCompleted(self)
        % Marks the chain as completed, and thus not modifiable
            assert(self.tableFull, 'The chain is not complete');
            self.completed = true;
        end

        function w = word(self, g)
        % Returns the word corresponding to the given permutation
        %
        % Args:
        %   g (permutation): Element of `.group`
        %
        % Returns:
        %   integer(1,\*): Letters of the word representing ``g``
            assert(self.completed);
            w = [];
            for i = 1:self.k
                beta = self.B(i);
                b = g(beta);
                ind = self.iOrbit(b,i);
                assert(ind > 0, 'Permutation not contained in the group');
                nu = self.nu{i}(:,ind)';
                nuw = self.nuw{i}{ind};
                g = nu(g);
                w = replab.fp.Letters.compose(w, -fliplr(nuw));
            end
        end

        function l = maximumWordLength(self)
        % Returns the maximal length of a word stored in this chain
        %
        % Returns:
        %   integer: Maximal length that can be returned by `.word`
            l = sum(cellfun(@(x) max(cellfun(@length, x)), self.nuw));
        end

        function s = orbitSizes(self)
        % Returns the size of orbits for each level
        %
        % Returns:
        %   integer(1,\*): List of orbit sizes
            s = cellfun(@length, self.orbit);
        end

        function o = order(self)
        % Returns the order of the group stored in this (possibly partial) chain
        %
        % Returns:
        %   vpi: Order
            o = replab.util.multiplyIntegers(self.orbitSizes);
        end

        function l = tableFull(self)
        % Returns whether the chain is complete
        %
        % See Minkwitz, Eq. (1.5), p. 90
        %
        % It is faster to check if the size of orbits matches, as we reuse the base of the standard stabilizer
        % computed for `.group`.
        %
        % Returns:
        %   logical: True if the order of this (possibly partial) chain equals the order of the underlying `.group`
            l = all(self.orbitSizes == self.chain.orbitSizes);
            % this is effectively ``l = self.order == self.group.order``
        end

        function [r rw] = step(self, i, t, tw)
        % Step procedure from Minkwitz
        %
        % Implemented as presented in Minkwitz, p. 91
        %
        % Args:
        %   i (integer): Level of the stabilizer chain at which to perform the step
        %   t (permutation): Permutation to examine
        %   tw (integer(1,\*)): Word form of the permutation ``t``
            assert(~self.completed);
            beta = self.B(i);
            img = t(beta);
            ind = self.iOrbit(img, i);
            if ind > 0
                v = self.nu{i}(:,ind)';
                vw = self.nuw{i}{ind};
                r = v(t);
                rw = replab.fp.Letters.compose(vw, tw);
                if length(tw) < length(vw)
                    t_inv = self.group.inverse(t);
                    tw_inv = -fliplr(tw);
                    self.newOrbit(img, i) = true;
                    self.nu{i}(:,ind) = t_inv';
                    self.nuw{i}{ind} = tw_inv;
                    self.step(i, t_inv, tw_inv);
                end
            else
                ind = length(self.orbit{i}) + 1;
                self.orbit{i}(ind) = img;
                self.iOrbit(img, i) = ind;
                self.newOrbit(img, i) = true;
                t_inv = self.group.inverse(t);
                tw_inv = -fliplr(tw);
                self.nu{i}(:,ind) = t_inv';
                self.nuw{i}{ind} = tw_inv;
                self.step(i, t_inv, tw_inv);
                r = 1:self.n;
                rw = [];
            end
        end

        function [t tw] = round(self, l, c, t, tw)
        % Round procedure
        %
        % Implemented as presented in Minkwitz, p. 91
        %
        % Args:
        %   l (integer): Maximal word length
        %   c (integer): Level of the chain at which to start the round
        %   t (permutation): Permutation to examine
        %   tw (integer(1,\*)): Word form of the permutation ``t``
        %
        % Returns
        % -------
        %   t:
        %     permutation: Part of the permutation that could not be sifted
        %   tw:
        %     integer(1,\*): Word form of the residual permutation
            assert(~self.completed);
            i = c;
            while ~self.group.isIdentity(t) && length(tw) < l
                beta = self.B(i);
                [r rw] = self.step(i, t, tw);
                assert(r(beta) == beta);
                t = r;
                tw = rw;
                i = i + 1;
            end
        end

        function w = next(self, w)
        % Enumerates words in the generators, by increasing word length
        %
        % Sketched in Minkwitz, p. 91
        %
        % We implement it by taking the previous word instead of a count, and simply increment the indices of the
        %  generators in the word. We also use inverses of the generators: the enumeration of an index goes
        % as ``1, ..., nG, -1, ..,. -nG``.
        %
        % When running out of digits, we restart the process with a word of length increased by one.
        %
        % We also take care of generating reduced words, so that we do not generate words of the form ``[ ... i -i ...]``.
        %
        % Args:
        %   w (integer(1,\*)): Current word in the enumeration
        %
        % Returns:
        %   integer(1,\*): Next word in the enumeration
            assert(~self.completed);
            nG = self.group.nGenerators;
            L = length(w);
            l = L;
            while l > 0
                if w(l) > -nG
                    if w(l) == nG
                        w(l) = -1;
                    elseif w(l) > 0
                        w(l) = w(l) + 1;
                    else
                        w(l) = w(l) - 1;
                    end
                    if (l == 1 || w(l-1) ~= -w(l)) && (l == L || w(l) ~= w(l+1))
                        break
                    end
                else
                    w(l) = 1;
                    l = l - 1;
                end
            end
            if l == 0
                w = ones(1, length(w) + 1);
            end
        end

        function improve(self, l)
        % Checks if the new entries computed so far can be employed to shorten the words already known
        %
        % Follows the pseudocode in Minkwitz, p. 93
        %
        % Args:
        %   l (integer): Maximal word length
            assert(~self.completed);
            while any(self.newOrbit(:))
                new = self.newOrbit;
                self.newOrbit = false(self.n, self.k);
                for j = 1:self.k
                    for ox = 1:length(self.orbit{j})
                        ix = self.orbit{j}(ox);
                        x = self.nu{j}(:,ox)';
                        xw = self.nuw{j}{ox};
                        for oy = 1:length(self.orbit{j})
                            iy = self.orbit{j}(oy);
                            y = self.nu{j}(:,oy)';
                            yw = self.nuw{j}{oy};
                            if new(ix, j) || new(iy, j)
                                t = x(y); % we do not need to inverse the product (left/right action thing)
                                tw = replab.fp.Letters.compose(xw, yw);
                                self.round(l, j, t, tw);
                            end
                        end
                    end
                end
            end
        end

        function fillOrbits(self, l)
        % Fill the current orbits
        %
        % Follows the pseudocode in Minkwitz, p. 93
        %
        % Note that there is a typo in the pseudocode, nu_i(p^{x^-1}) must be read nu_i^{-1}(p^{x^-1})
        %
        % Args:
        %   l (integer): Maximal word length
            assert(~self.completed);
            for i = 1:self.k
                Oi = self.orbit{i};
                for ii = 1:length(Oi)
                    bi = Oi(ii);
                    for j = i+1:self.k
                        Oj = self.orbit{j};
                        for jj = 1:length(Oj)
                            x = self.nu{j}(:,jj)';
                            xw = self.nuw{j}{jj};
                            p = x(bi);
                            if self.iOrbit(p, i) == 0
                                nuii = self.nu{i}(:,ii);
                                assert(nuii(bi) == self.B(i));
                                nuiiw = self.nuw{i}{ii};
                                t = self.group.composeWithInverse(x, nuii);
                                tw = replab.fp.Letters.compose(xw, -fliplr(nuiiw));
                                if length(tw) < l
                                    ind = length(self.orbit{i}) + 1;
                                    self.orbit{i}(ind) = p;
                                    self.iOrbit(p, i) = ind;
                                    self.newOrbit(p, i) = true;
                                    t_inv = self.group.inverse(t);
                                    assert(t_inv(p) == self.B(i));
                                    tw_inv = -fliplr(tw);
                                    self.nu{i}(:,ind) = t_inv';
                                    self.nuw{i}{ind} = tw_inv;
                                end
                            end
                        end
                    end
                end
            end
        end

        function sgsWord(self, l, nRounds)
        % Simplified generation procedure from Minkwitz
        %
        % Pseudocode in Minkwitz, p. 92
        %
        % Args:
        %   l (integer): Maximal word length
        %   nRounds (integer): Minimal number of rounds to perform
        %
        % Returns:
        %   integer(1,\*): The last word examined
            assert(~self.completed);
            tw = [];
            count = 0;
            while count < nRounds || ~self.tableFull
                tw = self.next(tw);
                t = self.group.imageLetters(tw);
                self.round(l, 1, t, tw);
                count = count + 1;
            end
        end

        function l = sgsWordQuick(self, nRounds, s, l)
        % Complete chain completion algorithm
        %
        % Follows closely the pseudocode in Minkwitz, pp. 92-93.
        %
        % The maximal word length is adjusted during the computation.
        %
        % Args:
        %   nRounds (integer): Minimal number of rounds to perform
        %   s (integer, optional): Perform the improvement procedure every ``s`` rounds
        %   l (integer): Maximal word length
        %
        % Returns:
        %   integer: Current word length
            assert(~self.completed);
            if nargin < 2 || isempty(nRounds)
                nRounds = 10000;
            end
            if nargin < 3 || isempty(s)
                s = self.k*self.k;
            end
            if nargin < 4 || isempty(l)
                l = 5;
            end
            tw = [];
            count = 0;
            while count < nRounds || ~self.tableFull
                tw = self.next(tw);
                t = self.group.imageLetters(tw);
                self.round(l, 1, t, tw);
                if mod(count, s) == 0
                    self.improve(l);
                    if ~self.tableFull
                        self.fillOrbits(l);
                        l = ceil(5*l/4);
                    end
                end
                count = count + 1;
            end
        end

    end

end
