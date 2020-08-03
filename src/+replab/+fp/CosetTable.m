classdef CosetTable < replab.Str
% Generates a coset table from the generators and relators of a finitely presented group and subset
%
% Based on Todd-Coxeter algorithm for coset enumeration from:
% Holt, Derek. "Coset Enumeration." Handbook of Computational Group Theory, Chapman & Hall/CRC, 2004, pp. 149–198
%
% Note that the action on cosets is a right action, contrary to the conventions used everywhere else in RepLAB.
% We usually indicate when we convert between conventions.
%
% The coset table is a mutable object.
%
% This class is used for both relators-based coset enumeration, and coset-table based coset enumeration. In the
% latter case, we need to keep track of deductions -- this happens when the flag `.storeDeductions` is true,
% this flag can be changed manually after construction of the coset table.
%
% The bounds `.maxCosets` and `.maxDeductions` are set to the default values provided by the RepLAB globals. These values
% can be modified after construction of the coset table.
%
% Words (for relators, subgroup generators, ...) are always written as integer row vectors, where ``1 ... nGenerators``
% represent the group generators. The generator inverses can either be written using ``-1 ... -nGenerators``, which is
% an external convention. However, in the coset table, the columns are numbered with generators and then their inverses.
% We append a ``L`` or ``C`` to distinguish between the former and latter conventions.

    properties
        nGenerators % (integer): Number of generators
        C % (integer(\*,\*)): Coset table, rows are cosets and columns are generators with their inverses
        n % (integer): Largest coset
        p % (integer(1,n)): Map of coincidences, ``p(i) <= i`` gives the canonical coset in case of coincidences
        maxCosets % (integer): Upper bound on the number on the cosets
        columnInverse % (integer(1,\*)): Inverse for each column
        maxDeductions % (integer): Maximum number of deductions to store on the stack
        storeDeductions % (logical): Whether to store deductions
        deductions % (integer(2,\*)): Deduction stack containing pairs (alpha, x)
    end

    properties (Access = protected)
        toMerge % (integer(1,\*)): Array of elements to merge, corresponds to the ``q`` variable in Holt
        toMergeFirst % (integer): Index of first element in queue
        toMergeLast % (integer): Index of last element in queue
    end

    methods (Static)

        function ct = fromGroupAndSubgroup(group, subgroup)
        % Constructs a coset table from a given group and a given subgroup
        %
        % Args:
        %   group (`+replab.FiniteGroup`): Group
        %   subgroup (`+replab.FiniteGroup`): Subgroup of ``group``
        %
        % Returns:
        %   `.CosetTable`: A standardized coset table corresponding to the right cosets of ``subgroup`` in ``group``
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

        function relators = presentation(group, subgroup, subgroupRelators)
        % Creates a presentation for a group based on the presentation for a subgroup
        %
        % We require ``group.generators(1:s) == subgroup.generators``, where ``s = subgroup.nGenerators``.
        %
        % PRESENTATION from Holt, p. 202
        %
        % Args:
        %   group (`+replab.FiniteGroup`): Group
        %   subgroup (`+replab.FiniteGroup`): Subgroup
        %   subgroupRelators (cell(1,\*) of integer(1,\*)): Relators defining a presentation of subgroup
        %
        % Returns:
        %   cell(1,\*) of integer(1,\*): Relators for ``group``
            r = group.nGenerators;
            s = subgroup.nGenerators;
            relators = subgroupRelators;
            rc = replab.fp.RelatorConjugates.fromRelators(r, relators);
            C = replab.fp.CosetTable.fromGroupAndSubgroup(group, subgroup);
            % line 1
            if r == s
                return
            end
            % line 2
            n = double(group.order/subgroup.order);
            Chat = replab.fp.CosetTable(r);
            Chat.storeDeductions = true;
            Chat.maxCosets = n;
            % line 3
            for i = 1:s
                Chat.C(1,i) = 1;
                Chat.C(1,i+r) = 1;
            end
            % line 4
            for i = 2:n
                [i1, x] = min(C.C(i,:));
                Chat.define(i1, Chat.columnInverse(x));
            end
            genIndex = [1:r -(1:r)];
            % line 5
            m = subgroup.abstractGroupIsomorphism;
            while true
                % transpose so that it is the smallest beta
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
                uhat_inH = -fliplr(m.target.toLetters(m.imageElement(u))); % flip due to inverse convention
                newRelator = [tbeta genIndex(x) -fliplr(tbetax) -fliplr(uhat_inH)];
                assert(group.isIdentity(group.imageLetters(newRelator)));
                relators{1,end+1} = newRelator;
                rc = rc.updated({newRelator});
                newRelatorC = Chat.lettersToColumnIndices(newRelator);
                Chat.scan(1, newRelatorC);
                % run coset-table based enumeration
                relatorsC = cellfun(@(w) Chat.lettersToColumnIndices(w), relators, 'uniform', 0);
                groupedRelatorsC = cellfun(@(list) cellfun(@(w) Chat.lettersToColumnIndices(w), list, 'uniform', 0), rc.groupedRelators, 'uniform', 0);
                Chat.processDeductions(relatorsC, groupedRelatorsC);
            end
        end

        function ct = cosetEnumerationR(nGenerators, relators, subgroupGenerators)
        % Enumerates cosets (relator based method)
        %
        % Taken from COSETENUMERATIONR, p. 164 of Holt
        %
        % Args:
        %   nGenerators (integer): Number of generators
        %   relators (cell(1,\*) of integer(1,\*)): Relators given as letter arrays
        %   subgroupGenerators (cell(1,\*) of integer(1,\*)): Subgroup generators given as letter arrays
        %
        % Returns:
        %   `.CosetTable`: The enumerated, standardized coset table
            ct = replab.fp.CosetTable(nGenerators);
            for i = 1:length(subgroupGenerators) % Line 3
                wC = ct.lettersToColumnIndices(subgroupGenerators{i});
                % Enter the generators of the subgroup as relations (the generators are part of the identity coset!)
                ct.scanAndFill(1, wC);
            end
            relatorsC = cellfun(@(w) ct.lettersToColumnIndices(w), relators, 'uniform', 0);
            alpha = 1; % since p(1) is always 1
            while alpha <= ct.n % Line 4
                if ct.p(alpha) ~= alpha % only consider live cosets
                    alpha = alpha + 1;
                    continue
                end
                for i = 1:length(relators) % Line 5
                    wC = relatorsC{i};
                    % Lines 6-7
                    ct.scanAndFill(alpha, wC);
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
        % Enumerates cosets (coset-table based method)
        %
        % Taken from COSETENUMERATIONC, p. 166 of Holt
        %
        % Args:
        %   nGenerators (integer): Number of generators
        %   relators (cell(1,\*) of integer(1,\*)): Relators given as letter arrays
        %   subgroupGenerators (cell(1,\*) of integer(1,\*)): Subgroup generators given as letter arrays
        %
        % Returns:
        %   `.CosetTable`: The enumerated, standardized coset table
            ct = replab.fp.CosetTable(nGenerators);
            ct.storeDeductions = true;
            rc = replab.fp.RelatorConjugates.fromRelators(nGenerators, relators);
            relatorsC = cellfun(@(w) ct.lettersToColumnIndices(w), relators, 'uniform', 0);
            groupedRelatorsC = cellfun(@(list) cellfun(@(w) ct.lettersToColumnIndices(w), list, 'uniform', 0), rc.groupedRelators, 'uniform', 0);
            for i = 1:length(subgroupGenerators)
                wC = ct.lettersToColumnIndices(subgroupGenerators{i});
                ct.scanAndFill(1, wC);
            end
            ct.processDeductions(relatorsC, groupedRelatorsC);
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
                        ct.processDeductions(relatorsC, groupedRelatorsC);
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
            self.maxCosets = replab.globals.maxCosets;
            self.maxDeductions = replab.globals.maxDeductions;
            self.storeDeductions = false;
            self.deductions = zeros(2, 0);
            self.columnInverse = [nGenerators+(1:nGenerators) 1:nGenerators];
            % initialize the toMerge queue
            m = 16;
            self.toMerge = zeros(1, m);
            self.toMergeFirst = 0;
            self.toMergeLast = 0;
        end

    end

    methods % Conversion between letters and columns
        function l = lettersToColumnIndices(self, l)
        % Converts a word given in letters to the convention used to index the columns of the coset table
        %
        % Args:
        %   l (integer(1,\*)): Letters in ``-n,...,-1`` and ``1,...,n`` with ``n`` the number of generators
        %
        % Returns:
        %   integer(1,\*): Expression referencing the columns of the coset table
            m = l < 0;
            l(m) = self.columnInverse(-l(m));
        end

        function l = columnIndicesToLetters(self, l)
        % Converts a word given using column indices to standard letters
        %
        % Args:
        %   l (integer(1,\*)): Expression referencing the columns of the coset table
        %
        % Returns:
        %   integer(1,\*): Letters in ``-n,...,-1`` and ``1,...,n`` with ``n`` the number of generators
            m = l > self.nGenerators;
            l(m) = -self.columnInverse(l(m));
        end

    end

    methods % Operations that do not modify the table

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
                w = [self.columnInverse(x) w];
                i = i1;
            end
            nG = self.nGenerators;
            w(w > nG) = -(w(w > nG) - nG);
        end

    end

    methods

        function newIndex = compress(self)
        % Removes repeated cosets from table
        %
        % COMPRESS, Holt p. 167
            g = 0;
            newIndex = zeros(1, self.n);
            for alpha = 1:self.n
                if self.p(alpha) == alpha
                    g = g + 1;
                    newIndex(alpha) = g;
                    if g ~= alpha
                        for x = 1:self.nGenerators*2
                            beta = self.C(alpha, x);
                            if beta == alpha
                                beta = g;
                            end
                            self.C(g, x) = beta;
                            if beta ~= 0
                                self.C(beta, self.columnInverse(x)) = g;
                            end
                        end
                    end
                end
            end
            self.n = g;
            self.p = 1:self.n;
            self.C = self.C(1:self.n, :);
        end

        function switchElmtsSlow(self, beta, g)
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

        function switchElmts(self, beta, gamma)
        % Switch the given coset numbers
        %
        % Faster implementation of `.switchElmtsSlow`
        %
        % Args:
        %   beta (integer): First coset number to switch
        %   gamma (integer): Second coset number to switch
            beta1 = self.C(gamma, :);
            gamma1 = self.C(beta, :);
            beta1(self.C(gamma, :) == beta) = gamma;
            beta1(self.C(gamma, :) == gamma) = beta;
            gamma1(self.C(beta, :) == beta) = gamma;
            gamma1(self.C(beta, :) == gamma) = beta;
            self.C(gamma, :) = gamma1;
            self.C(beta, :) = beta1;
            for x = 1:self.nGenerators*2
                betax = beta1(x);
                gammax = gamma1(x);
                self.C(betax, self.columnInverse(x)) = beta;
                self.C(gammax, self.columnInverse(x)) = gamma;
            end
        end

        function standardize(self)
        % Standardizes coset table such that elements appear in order they would with no coincidences
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
            if self.n == self.maxCosets
                error('Error: exceeded maximum allowed number of cosets')
            end
            nG = self.nGenerators;
            self.n = self.n + 1;
            if self.n > size(self.C, 1) % grow the table if needed
                self.C = [self.C; zeros(size(self.C, 1), self.nGenerators * 2)];
                self.p = [self.p zeros(1, length(self.p))];
            end
            beta = self.n;
            self.p(beta) = beta;
            self.C(alpha, x) = beta;
            xinv = self.columnInverse(x);
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

        function merge(self, k, lambda)
        % Removes the smaller of k and lambda from active cosets
        %
        % MERGE, p. 157
        %
        % Args:
        %   k (integer): coset coincident with lambda
        %   lambda (integer): coset coincident with k
            lam1 = self.rep(k);
            lam2 = self.rep(lambda);
            if lam1 ~= lam2
                mu = min([lam1, lam2]);
                v = max([lam1, lam2]);
                self.p(v) = mu;
                m = length(self.toMerge);
                if self.toMergeLast == m
                    if self.toMergeFirst > m/2
                        self.toMerge(1:m/2) = self.toMerge(m/2+1:m);
                        self.toMergeFirst = self.toMergeFirst - m/2;
                        self.toMergeLast = self.toMergeLast - m/2;
                    else
                        self.toMerge = [self.toMerge zeros(1, m)];
                    end
                end
                ind = self.toMergeLast + 1;
                self.toMergeLast = ind;
                self.toMerge(ind) = v;
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
            self.toMergeFirst = 1;
            self.toMergeLast = 0;
            self.merge(alpha, beta);
            while self.toMergeFirst <= self.toMergeLast
                gamma = self.toMerge(self.toMergeFirst);
                self.toMergeFirst = self.toMergeFirst + 1;
                for x = 1:self.nGenerators*2
                    delta = self.C(gamma, x);
                    if delta > 0
                        xinv = self.columnInverse(x);
                        self.C(delta, xinv) = 0;
                        if self.storeDeductions
                            self.deductions(:,end+1) = [delta; xinv];
                        end
                        mu = self.rep(gamma);
                        nu = self.rep(delta);
                        Cmux = self.C(mu, x);
                        if Cmux > 0
                            self.merge(nu, self.C(mu, x));
                        else
                            Cnuxinv = self.C(nu, xinv);
                            if Cnuxinv > 0
                                self.merge(mu, Cnuxinv);
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

    methods % Methods that work with relators

        function scan(self, alpha, wC)
        % Scans through a relator
        %
        % SCAN, Holt p. 155
        %
        % Args:
        %   alpha (integer): Initial coset number
        %   wC (integer(1,\*)): Word to scan, whose letters are column indices
            f = alpha; % line 2
            r = length(wC);
            colinv = self.columnInverse;
            i = 1;
            while i <= r
                CfwC = self.C(f, wC(i)); % line 3
                if CfwC == 0
                    break
                end
                f = CfwC; % line 4
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
            while j >= i
                CbwinvC = self.C(b, colinv(wC(j)));
                if CbwinvC == 0 % line 9
                    break
                end
                b = CbwinvC; % line 10
                j = j - 1;
            end
            if j < i % line 11
                self.coincidence(f, b); % line 12
            elseif j == i % line 13
                self.C(f, wC(i)) = b; % line 14
                self.C(b, colinv(wC(i))) = f;
                if self.storeDeductions
                    self.deductions(:,end+1) = [f; wC(i)];
                end
            end
            % otherwise j > i, scan is incomplete and yields no information
        end

        function lookahead(self, relatorsC)
        % Scan all cosets under all relators, without making new definitions
        %
        % LOOKAHEAD, Holt, p. 164
        %
        % Args:
        %   relatorsC (cell(1,\*) of integer(1,\*)): Relators given as letters in the column indices
            n = self.n;
            for beta = 1:n
                if self.p(beta) == beta % in Omega, live
                    for i = 1:length(relatorsC)
                        wC = relatorsC{i};
                        self.scan(beta, wC);
                        if self.p(beta) < beta
                            break
                        end
                    end
                end
            end
        end

        function processDeductions(self, relatorsC, groupedRelatorsC)
        % Process the deduction table
        %
        % PROCESSDEDUCTIONS, Holt, p. 165
        %
        % Arguments are given using columns indices as letters.
        %
        % Args:
        %   relatorsC (cell(1,\*) of integer(1,\*)): Relators given as letters in the column indices
        %   groupedRelatorsC (cell(1,\*) of cell(1,\*) of integer(1,\*)): Relators and their conjugates, grouped by their first letter
            while size(self.deductions, 2) > 0
                % If the deduction table is full, use lookahead and clear the table
                if size(self.deductions, 2) >= self.maxDeductions
                    self.lookahead(relatorsC);
                    self.deductions = zeros(2, 0);
                else
                    alpha = self.deductions(1, end);
                    x = self.deductions(2, end);
                    self.deductions = self.deductions(:, 1:end-1);
                    if self.p(alpha) == alpha % is active?
                        rels = groupedRelatorsC{x};
                        for i = 1:length(rels)
                            wC = rels{i};
                            self.scan(alpha, wC);
                            if self.p(alpha) < alpha
                                break
                            end
                        end
                    end
                    beta = self.C(alpha, x);
                    if beta ~= 0 && self.p(beta) == beta % is active
                        xinv = self.columnInverse(x);
                        rels = groupedRelatorsC{xinv};
                        for i = 1:length(rels)
                            wC = rels{i};
                            self.scan(beta, wC);
                            if self.p(beta) < beta
                                break
                            end
                        end
                    end
                end
            end
        end

        function scanAndFill(self, alpha, wC)
        % Scans through a relator and fills in table
        %
        % SCANANDFILL, Holt p. 163
        %
        % Args:
        %   alpha (integer): Initial coset number
        %   wC (integer(1,\*)): Word whose letters are column indices
            r = length(wC);
            colinv = self.columnInverse;
            f = alpha;
            i = 1;
            b = alpha;
            j = r;
            while true
                while i <= r
                    CfwC = self.C(f, wC(i)); % line 3
                    if CfwC == 0
                        break
                    end
                    f = CfwC; % line 4
                    i = i + 1;
                end
                if i > r
                    if f ~= alpha
                        self.coincidence(f, alpha);
                    end
                    return
                end
                while j >= i
                    CbwinvC = self.C(b, colinv(wC(j)));
                    if CbwinvC == 0 % line 9
                        break
                    end
                    b = CbwinvC; % line 10
                    j = j - 1;
                end
                if j < i
                    self.coincidence(f, b);
                    return % not in Holt, but is needed
                elseif j == i
                    self.C(f, wC(i)) = b;
                    self.C(b, colinv(wC(j))) = f;
                    if self.storeDeductions
                        self.deductions(:,end+1) = [f; wC(i)];
                    end
                    return
                else
                    self.define(f, wC(i));
                end
            end
        end

    end

end
