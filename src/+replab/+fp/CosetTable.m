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
        relators % (cell(1,\*) of integer(1,\*)): Relators given as reduced words
        subgroupGenerators % (cell(1,\*) of integer(1,\*)): Words representing the generators of the subgroup
        C % (integer(\*,\*)): Coset table, rows are cosets and columns are generators with their inverses
        n % (integer): Largest coset
        p % (integer(1,n)): Map of coincidences, ``p(i) <= i`` gives the canonical coset in case of coincidences
        M % (integer): Upper bound on the number on the cosets
    end


    methods

        function self = CosetTable(nGenerators, relators, subgroupGenerators, M)
        % Constructs an empty coset table
        %
        % Args:
        %   nGenerators (integer): Number of generators
        %   relators (cell(1,\*) of integer(1,\*)): Relators written as letter arrays
        %   subgroupGenerators (cell(1,\*) of integer(1,\*)): Subgroup generators written as letter arrays
        %   M (integer): Upper bound on the number of cosets during enumeration
            self.nGenerators = nGenerators;
            self.relators = relators;
            self.subgroupGenerators = subgroupGenerators;
            self.C = zeros(1, nGenerators*2);
            self.n = 1;
            self.p = [1];
            self.M = M;
        end

        function cosetEnumerationR(self)
        % Enumerates cosets (relator based method)
        %
        % Taken from COSETENUMERATIONR, p. 164 of Holt
            for i = 1:length(self.subgroupGenerators) % Line 3
                w = self.subgroupGenerators{i};
                % Enter the generators of the subgroup as relations (the generators are part of the identity coset!)
                self.scanAndFill(1, w);
            end
            alpha = 1; % since p(1) is always 1
            while alpha <= self.n % Line 4
                if self.p(alpha) ~= alpha % only consider live cosets
                    alpha = alpha + 1;
                    continue
                end
                for i = 1:length(self.relators) % Line 5
                    w = self.relators{i};
                    % Lines 6-7
                    self.scanAndFill(alpha, w);
                    if self.p(alpha) < alpha
                        break
                    end
                end
                if self.p(alpha) < alpha % Line 8
                    alpha = alpha + 1;
                    continue
                end
                for x = 1:self.nGenerators * 2 % Lines 9-10
                    if self.C(alpha, x) < 1
                        self.define(alpha, x);
                    end
                end
                alpha = alpha + 1;
            end
            % Clean up the coset table
            f = find(self.p - (1:self.n) == 0); % @jvdv37: what's the rationale here?
            self.C = self.C(1:length(f), :);
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

        function define(self, alpha, x)
        % Defines a new coset
        %
        % DEFINE, Holt p. 153
        %
        % Args:
        %   alpha (integer): TODO
        %   x (integer): TODO
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
        end

        function lambda = rep(self, k)
        % Replaces a coset in p
        %
        % REP, Holt p. 157
        %
        % Args:
        %   k (integer): TODO
        %
        % Returns:
        %   integer: TODO
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

        function scanAndFill(self, alpha, w)
        % Scans through a relator and fills in table
        %
        % Args:
        %   alpha (integer): TODO
        %   w (integer(1,\*)): Word letters
            s = sign(w);
            ngens = self.nGenerators;
            f = find(w < 0);
            w(f) = abs(w(f)) + ngens; % @jvdv37 optimized here!
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
                while j >= i && self.C(b, w(j) + s(j)*ngens) > 0
                    b = self.C(b, w(j) + s(j)*ngens);
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
                    self.C(b, w(j) + s(j)*ngens) = f;
                    return
                else
                    self.define(f, w(i));
                end
            end
        end

        function [q, l] = merge(self, k, lambda, q, l)
        % Removes the smaller of k and lambda from active cosets
        %
        % MERGE, p. TODO
        %
        % Args:
        %   k (integer): TODO
        %   lambda (integer): TODO
        %   q (integer(1,\*)): TODO
        %   l (integer): TODO
            lam1 = self.rep(k);
            lam2 = self.rep(lambda);
            if lam1 ~= lam2
                mu = min([lam1, lam2]);
                v = max([lam1, lam2]);
                self.p(v) = mu;
                l = l + 1;
                q(l) = v;
            end
        end

        function coincidence(self, alpha, beta)
        % Deals with the same coset being represented with different values
        %
        % COINCIDENCE p. TODO
        %
        % Args:
        %   alpha (integer): TODO
        %   beta (integer): TODO
            ngens = self.nGenerators;
            l = 0;
            q = [];
            [q, l] = self.merge(alpha, beta, q, l);
            i = 1;
            while i <= l
                g = q(i);
                i = i + 1;
                gx = self.C(g, :);
                for j = 1:self.nGenerators*2
                    if gx(j) > 0
                        delta = gx(j);
                        if j <= ngens
                            jinv = j + ngens;
                        else
                            jinv = j - ngens;
                        end
                        self.C(delta, jinv) = 0;
                        mu = self.rep(g);
                        nu = self.rep(delta);
                        if self.C(mu, j) > 0
                            [q, l] = self.merge(nu, self.C(mu, j), q, l);
                        elseif self.C(nu, jinv) > 0
                            [q, l] = self.merge(mu, self.C(nu, jinv), q, l);
                        else
                            self.C(mu, j) = nu;
                            self.C(nu, jinv) = mu;
                        end
                    end
                end
            end
        end

    end

end
