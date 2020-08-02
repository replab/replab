classdef CosetTable < replab.Str
% Generates a coset table from the generators and relators of a finitely presented group and subset
%
% Based on Todd-Coxeter algorithm for coset enumeration from:
% Holt, Derek. “Coset Enumeration.” Handbook of Computational Group Theory,
% Chapman & Hall/CRC, 2004, pp. 149–198
%
% Note that the action on cosets is a right action, contrary to the conventions used everywhere else in RepLAB.

    properties
        nGenerators % (integer): Number of generators
        C % (integer(\*,\*)): Coset table, rows are cosets and columns are generators with their inverses
        n % (integer): Largest coset
        p % (integer(1,n)): Map of coincidences, ``p(i) <= i`` gives the canonical coset in case of coincidences
        M % (integer): Upper bound on the number on the cosets
        maxDeductions % (integer): Maximum number of deductions to store on the stack
        storeDeductions % (logical): Whether to store deductions
        deductions % (integer(2,\*)): Deduction stack containing pairs (alpha, x)
    end

    methods (Static)

        function r = cyclicallyReduce(r)
        end

        function c = cyclicConjugates(relator)
            r = replab.fp.CosetTable.cyclicallyReduce
        end

        function ct = fromRelatorsAndSubgroupGenerators(nGenerators, relators, subgroupGenerators)
            ct = replab.fp.CosetTable(nGenerators, relators, subgroupGenerators, 2^50);
            ct.cosetEnumerationR;
        end

        function ct = fromGroupAndSubgroup(group, subgroup)
            nG = group.nGenerators;
            ct = replab.fp.CosetTable(nG);
            lc = group / subgroup;
            la = lc.leftAction;
            nC = double(lc.nElements);
            ct.C = zeros(nC, nG*2);
            for i = 1:nG
                gen = group.generator(i);
                invGen = group.generatorInverse(i);
                ct.C(:,i+nG) = la.imageElement(gen); % invert due to right action vs left action
                ct.C(:,i) = la.imageElement(invGen);
            end
            ct.n = nC;
            ct.p = 1:nC;
            ct.standardize;
        end

        function relators = presentation(group, subgroup, s, subgroupRelators)
        % Creates a presentation for a group based on the presentation for a subgroup
        %
        % We require ``group.generators(1:s) == subgroup.generators``.
        %
        % Args:
        %   group (`+replab.FiniteGroup`): Group
        %   subgroup (`+replab.FiniteGroup`): Subgroup
        %   s (integer): Number of generators of ``group`` that match the generators of ``subgroup``
        %   subgroupRelators (cell(1,\*) of integer(1,\*)): Relators defining a presentation of subgroup
            r = group.nGenerators;
            s = subgroup.nGenerators;
            relators = subgroupRelators;
            C = replab.fp.CosetTable.fromGroupAndSubgroup(group, subgroup);
            % line 1
            if r == s
                return
            end
            % line 2
            n = double(group.order/subgroup.order);
            Chat = replab.fp.CosetTable(r, 2*n);
            % line 3
            for i = 1:s
                Chat.C(1,i) = 1;
                Chat.C(1,i+r) = 1;
            end
            % line 4
            for i = 2:n
                [i1, x] = min(C.C(i,:));
                Chat.define(i1, Chat.generatorInverse(x));
            end
            genIndex = [1:r -(1:r)];
            % line 5
            m = subgroup.abstractGroupIsomorphism;
            while true
                Chat.n
                [x beta] = find((C.C ~= Chat.C(1:n, :))', 1);
                if isempty(beta)
                    break
                end
                betax = C.C(beta, x);
                tbeta = C.transversal(beta);
                tbetax = C.transversal(betax);
                uhat_inG = [tbeta genIndex(x) -fliplr(tbetax)];
                u = group.imageLetters(-fliplr(uhat_inG));
                assert(subgroup.contains(u));
                uhat_inH = -fliplr(m.target.toLetters(m.imageElement(u)));
                newRelator = [tbeta genIndex(x) -fliplr(tbetax) -fliplr(uhat_inH)];
                assert(group.isIdentity(group.imageLetters(newRelator)));
                relators{1,end+1} = newRelator;
                Chat.scanAndFill(1, newRelator);
                alpha = 1; % since p(1) is always 1
                while alpha <= Chat.n % Line 4
                    if Chat.p(alpha) ~= alpha % only consider live cosets
                        alpha = alpha + 1;
                        continue
                    end
                    for i = 1:length(relators) % Line 5
                        w = relators{i};
                        % Lines 6-7
                        try
                            Chat.scanAndFill(alpha, w);
                        catch
                        end
                        if Chat.p(alpha) < alpha
                            break
                        end
                    end
                    if Chat.p(alpha) < alpha % Line 8
                        alpha = alpha + 1;
                        continue
                    end
                    alpha = alpha + 1;
                end
            end
        end

        function ct = cosetEnumerationR(nGenerators, relators, subgroupGenerators)
        % Enumerates cosets (relator based method)
        %
        % Taken from COSETENUMERATIONR, p. 164 of Holt
            ct = replab.fp.CosetTable(nGenerators);
            for i = 1:length(subgroupGenerators) % Line 3
                w = subgroupGenerators{i};
                % Enter the generators of the subgroup as relations (the generators are part of the identity coset!)
                ct.scanAndFill(1, w);
            end
            alpha = 1; % since p(1) is always 1
            while alpha <= ct.n % Line 4
                if ct.p(alpha) ~= alpha % only consider live cosets
                    alpha = alpha + 1;
                    continue
                end
                for i = 1:length(relators) % Line 5
                    w = relators{i};
                    % Lines 6-7
                    ct.scanAndFill(alpha, w);
                    if ct.p(alpha) < alpha
                        break
                    end
                end
                if ct.p(alpha) < alpha % Line 8
                    alpha = alpha + 1;
                    continue
                end
                for x = 1:nGenerators * 2 % Lines 9-10
                    if ct.C(alpha, x) < 1
                        ct.define(alpha, x);
                    end
                end
                alpha = alpha + 1;
            end
            % Compress and standardize coset table
            ct.compress;
            ct.standardize;
        end

        function ct = cosetEnumerationC(nGenerators, relators, subgroupGenerators)
            ct = replab.fp.CosetTable(nGenerators);
            ct.storeDeductions = true;
            relators = cellfun(@(r) replab.fp.Letters.cyclicallyReduce(r), relators, 'uniform', 0); % line 1
            relators = relators(~cellfun(@isempty, relators));
            relators = replab.fp.Letters.unique(relators);
            relatorsC = cell(1, 0); % line 2
            for i = 1:length(relators)
                relatorsC = horzcat(relatorsC, replab.fp.Letters.cyclicConjugates(relators{i}));
                relatorsC = horzcat(relatorsC, replab.fp.Letters.cyclicConjugates(-fliplr(relators{i})));
            end
            relatorsC1 = replab.fp.Letters.unique(relatorsC);
            groupedRelators = cell(1, nGenerators*2);
            for i = 1:nGenerators*2
                groupedRelators{i} = cell(1, 0);
            end
            for i = 1:length(relatorsC1)
                r = relatorsC1{i};
                f = ct.lettersToColumnIndices(r(1));
                groupedRelators{f} = horzcat(groupedRelators{f}, r);
            end
            for i = 1:length(subgroupGenerators)
                w = subgroupGenerators{i};
                ct.scanAndFill(1, w);
            end
            ct.processDeductions(relators, groupedRelators);
            alpha = 1;
            while alpha <= ct.n
                if ct.p(alpha) ~= alpha % only consider live cosets
                    alpha = alpha + 1;
                    continue
                end
                for x = 1:nGenerators * 2
                    if ct.p(alpha) ~= alpha
                        break
                    end
                    if ct.C(alpha, x) == 0
                        ct.define(alpha, x);
                        ct.processDeductions(relators, groupedRelators);
                    end
                end
                alpha = alpha + 1;
            end
            ct.compress;
            ct.standardize;
        end

    end


    methods

        function self = CosetTable(nGenerators)
        % Constructs an empty coset table
        %
        % Args:
        %   nGenerators (integer): Number of generators
        %   M (integer, optional): Upper bound on the number of cosets during enumeration (default value `+replab.+globals.maxCosets`)
            self.nGenerators = nGenerators;
            self.C = zeros(1, nGenerators*2);
            self.n = 1;
            self.p = [1];
            self.M = replab.globals.maxCosets;
            self.maxDeductions = replab.globals.maxDeductions;
            self.storeDeductions = false;
            self.deductions = zeros(2, 0);
        end

        function w = transversal(self, beta)
        % Returns a transversal element
        %
        % Requires that the coset table has been standardized.
        %
        % Args:
        %   beta (integer): Coset number
        %
        % Returns:
        %   integer(1,\*): Word such that ``1^word = beta``
            w = zeros(1, 0);
            i = beta;
            while i ~= 1
                [i1, x] = min(self.C(i,:));
                w = [self.generatorInverse(x) w];
                i = i1;
            end
            nG = self.nGenerators;
            w(w > nG) = -(w(w > nG) - nG);
        end


        function invind = generatorInverse(self, ind)
        % Returns the index of the inverse of the given generator index
        %
        % Args:
        %   ind (integer): Index of a generator or its inverse, ``1 <= ind <= 2*nGenerators``
        %
        % Returns:
        %   integer: The index of the inverse of that generator
            nG = self.nGenerators;
            if ind <= nG
                invind = ind + nG;
            else
                invind = ind - nG;
            end
        end

        function compress(self)
        % Removes repeated cosets from table
        %
        % COMPRESS, Holt p. 167
            g = 0;
            for alpha = 1:self.n
                if self.p(alpha) == alpha
                    g = g + 1;
                    if g ~= alpha
                        for x = 1:self.nGenerators*2
                            beta = self.C(alpha, x);
                            if beta == alpha
                                beta = g;
                            end
                            self.C(g, x) = beta;
                            self.C(beta, self.generatorInverse(x)) = g;
                        end
                    end
                end
            end
            self.n = g;
            self.p = 1:self.n;
            self.C = self.C(1:self.n, :);
        end

        function switchElmts(self, beta, g)
        % Switch the elements beta and g
        %
        % SWITCH, Holt p. 167
            for x = 1:self.nGenerators*2
                z = self.C(g, x);
                self.C(g, x) = self.C(beta, x);
                self.C(beta, x) = z;
                for alpha = 1:self.n
                    if self.C(alpha, x) == beta
                        self.C(alpha, x) = g;
                    elseif self.C(alpha, x) == g
                        self.C(alpha, x) = beta;
                    end
                end
            end
        end

        function standardize(self)
        % Standardizes coset table such that elements appear in order they
        % would with no coincidences
        %
        % STANDARDIZE, Holt, p. 168
            g = 2;
            for alpha = 1:self.n
                for x = 1:self.nGenerators*2
                    beta = self.C(alpha, x);
                    if beta >= g
                        if beta > g
                            self.switchElmts(g, beta)
                        end
                        g = g + 1;
                        if g == self.n
                            return
                        end
                    end
                end
            end
        end

        function define(self, alpha, x)
        % Defines a new coset
        %
        % DEFINE, Holt p. 153
        %
        % Args:
        %   alpha (integer): next unused coset number
        %   x (integer): generator or inverse generator
            if self.n == self.M
                error('Error: exceeded maximum allowed number of cosets')
            end
            nG = self.nGenerators;
            self.n = self.n + 1;
            beta = self.n;
            self.p(beta) = beta;
            self.C(alpha, x) = beta;
            xinv = self.generatorInverse(x);
            self.C(beta, xinv) = alpha;
            if self.storeDeductions
                self.deductions(:,end+1) = [alpha; x];
            end
        end

        function lambda = rep(self, k)
        % Replaces a coset in p
        %
        % REP, Holt p. 157
        %
        % Args:
        %   k (integer): coset number to insert
        %
        % Returns:
        %   integer: coset number that was replaced
            lambda = k;
            rho = self.p(lambda);
            while rho ~= lambda
                lambda = rho;
                rho = self.p(lambda);
            end
            mu = k;
            rho = self.p(mu);
            while rho ~= lambda
                self.p(mu) = lambda;
                mu = rho;
                rho = self.p(mu);
            end
        end

        function l = lettersToColumnIndices(self, l)
            m = l < 0;
            l(m) = self.nGenerators + (-l(m));
        end

        function l = columnIndicesToLetters(self, l)
            m = l > self.nGenerators;
            l(m) = -(l(m) - self.nGenerators);
        end

        function scan(self, alpha, w)
        % Scans through a relator
        %
        % SCAN, Holt p. 155
        %
        % Args:
        %   alpha (integer): Initial coset number
        %   w (integer(1,\*)): Word letters
            w = self.lettersToColumnIndices(w); % line 1
            f = alpha; % line 2
            r = length(w);
            i = 1;
            while i <= r && self.C(f, w(i)) > 0 % line 3
                f = self.C(f, w(i)); % line 4
                i = i + 1;
            end
            if i > r % line 5
                if f ~= alpha % line 6
                    self.coincidence(f, alpha);
                end
                return % line 7
            end
            % Forward scan incomplete, scan backwards
            b = alpha; % line 8
            j = r;
            while j >= i && self.C(b, self.generatorInverse(w(j))) > 0 % line 9
                b = self.C(b, self.generatorInverse(w(j))); % line 10
                j = j - 1;
            end
            if j < i % line 11
                self.coincidence(f, b); % line 12
            elseif j == i % line 13
                self.C(f, w(i)) = b; % line 14
                self.C(b, self.generatorInverse(w(i))) = f;
                if self.storeDeductions
                    self.deductions(:,end+1) = [f; w(i)];
                end
            end
            % otherwise j > i, scan is incomplete and yields no information
        end

        function lookahead(self, relators)
        % Scan all cosets under all relators, without making new definitions
        %
        % LOOKAHEAD, Holt, p. 164
        %
        % Args:
        %   relators (cell(1,\*) of integer(1,\*)): Relators given as letters
            n = self.n;
            for beta = 1:n
                if self.p(beta) == beta % in Omega, live
                    for i = 1:length(relators)
                        w = relators{i};
                        self.scan(beta, w);
                        if self.p(beta) < beta
                            break
                        end
                    end
                end
            end
        end

        function processDeductions(self, relators, groupedRelators)
        % Process the deduction table
        %
        % PROCESSDEDUCTIONS, Holt, p. 165
        %
        % Args:
        %   relators (cell(1,\*) of integer(1,\*)): Relators given as letters
        %   groupedRelators (cell(1,\*) of cell(1,\*) of integer(1,\*)): Relators and their conjugates, grouped by their first letter (in column convention)
            while size(self.deductions, 2) > 0
                if size(self.deductions, 2) >= self.maxDeductions
                    self.lookahead(relators);
                    self.deductions = zeros(2, 0);
                else
                    alpha = self.deductions(1, end);
                    x = self.deductions(2, end);
                    self.deductions = self.deductions(:, 1:end-1);
                    if self.p(alpha) == alpha % is active?
                        rels = groupedRelators{x};
                        for i = 1:length(rels)
                            w = rels{i};
                            self.scan(alpha, w);
                            if self.p(alpha) < alpha
                                break
                            end
                        end
                    end
                    beta = self.C(alpha, x);
                    if beta ~= 0 && self.p(beta) == beta % is active
                        xinv = self.generatorInverse(x);
                        rels = groupedRelators{xinv};
                        for i = 1:length(rels)
                            w = rels{i};
                            self.scan(beta, w);
                            if self.p(beta) < beta
                                break
                            end
                        end
                    end
                end
            end
        end

        function scanAndFill(self, alpha, w)
        % Scans through a relator and fills in table
        %
        % SCANANDFILL, Holt p. 163
        %
        % Args:
        %   alpha (integer): initial coset number
        %   w (integer(1,\*)): word letters
            w = self.lettersToColumnIndices(w);
            r = length(w);
            f = alpha;
            i = 1;
            b = alpha;
            j = r;
            while true
                while i <= r && self.C(f, w(i)) > 0
                    f = self.C(f, w(i));
                    i = i + 1;
                end
                if i > r
                    if f ~= alpha
                        self.coincidence(f, alpha);
                    end
                    return
                end
                while j >= i && self.C(b, self.generatorInverse(w(j))) > 0
                    b = self.C(b, self.generatorInverse(w(j)));
                    j = j - 1;
                end
                if j < i
                    self.coincidence(f, b);
                    f = alpha;
                    i = 1;
                    b = alpha;
                    j = r;
                elseif j == i
                    self.C(f, w(i)) = b;
                    self.C(b, self.generatorInverse(w(j))) = f;
                    if self.storeDeductions
                        self.deductions(:,end+1) = [f; w(i)];
                    end
                    return
                else
                    self.define(f, w(i));
                end
            end
        end

        function q = merge(self, k, lambda, q)
        % Removes the smaller of k and lambda from active cosets
        %
        % MERGE, p. 157
        %
        % Args:
        %   k (integer): coset coincident with lambda
        %   lambda (integer): coset coincident with k
        %   q (integer(1,\*)): coset numbers that are replaced
            lam1 = self.rep(k);
            lam2 = self.rep(lambda);
            if lam1 ~= lam2
                mu = min([lam1, lam2]);
                v = max([lam1, lam2]);
                self.p(v) = mu;
                q(1,end+1) = v;
            end
        end

        function coincidence(self, alpha, beta)
        % Deals with the same coset being represented with different values
        %
        % COINCIDENCE p. 158
        %
        % Args:
        %   alpha (integer): coincident coset (arrive to this coset from one direction in the word)
        %   beta (integer): coset coicident with alpha (arrive to this coset from the other direction)
            q = [];
            q = self.merge(alpha, beta, q);
            while ~isempty(q)
                gamma = q(1);
                q = q(2:end);
                for x = 1:self.nGenerators*2
                    delta = self.C(gamma, x);
                    if delta > 0
                        xinv = self.generatorInverse(x);
                        self.C(delta, xinv) = 0;
                        self.deductions(:,end+1) = [delta; xinv];
                        mu = self.rep(gamma);
                        nu = self.rep(delta);
                        if self.C(mu, x) > 0
                            q = self.merge(nu, self.C(mu, x), q);
                        elseif self.C(nu, xinv) > 0
                            q = self.merge(mu, self.C(nu, xinv), q);
                        else
                            self.C(mu, x) = nu;
                            self.C(nu, xinv) = mu;
                        end
                    end
                end
            end
        end

    end

end
