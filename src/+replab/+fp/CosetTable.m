classdef CosetTable < replab.Str
% Generates a coset table from the generators and relators of a finitely presented group and subset
%
% Based on Todd-Coxeter algorithm for coset enumeration from:
% Holt, Derek. “Coset Enumeration.” Handbook of Computational Group Theory,
% Chapman & Hall/CRC, 2004, pp. 149–198
%
% Example:
%   >>> ct = replab.fp.CosetTable.original({'x' 'y'}, {'x^2','y^3','(x*y)^3'}, {'x*y'});
%   >>> ct.table
%         | x  y  inv(x)  inv(y)
%       --------------------------
%       1 | 2  3     2       2
%       2 | 1  1     1       3
%       3 | 4  2     4       1
%       4 | 3  4     3       4

    properties
        generatorNames % (cell(1,\*) of charstring): Names of the generators
        nGenerators % (integer): Number of generators
        relators % (cell(1,\*) of integer(1,\*)): Relators given as reduced words
        subgroupGenerators % (cell(1,\*) of integer(1,\*)): Words representing the generators of the subgroup
        C % (integer(\*,\*)): Coset table, rows are cosets and columns are generators with their inverses
        n % (integer): Largest coset
        p % (integer(1,n)): Map of coincidences, ``p(i) <= i`` gives the canonical coset in case of coincidences
        M % (integer): Upper bound on the number on the cosets

    end

    methods (Static)

        function ct = original(generators, relators, y)
        % Enumerate cosets
        %
        % Args:
        %   generators (cell(1,\*) of charstring): group generators
        %   relators (cell(1,\*) of charstring): group relators to parse
        %   y (cell(1, \*) of charstring): relators for a finite subset of the group
            ct = replab.fp.CosetTable(generators, 2^51 - 1); % TODO: estimate memory use
            ngens = length(generators);
            for i = 1:length(y)
                w = y{i};
                ct.scanAndFill(1, w);
            end
            alpha = 1; % since p(1) is always 1
            while alpha <= ct.n
                if ct.p(alpha) ~= alpha
                    alpha = alpha + 1;
                    continue
                end
                for i = 1:length(relators)
                    w = relators{i};
                    ct.scanAndFill(alpha, w);
                    if ct.p(alpha) < alpha
                        break
                    end
                end
                if ct.p(alpha) < alpha
                    alpha = alpha + 1;
                    continue
                end
                for x = 1:ngens * 2
                    if ct.C(alpha, x) < 1
                        ct.define(alpha, x);
                    end
                end
                alpha = alpha + 1;
            end
            ct.C
            f = find(ct.p - (1:ct.n) == 0);
            ct.C = ct.C(1:length(f), :);
        end

    end

    methods

        function self = CosetTable(generatorNames, M)
            nG = length(generatorNames);
            self.generatorNames = generatorNames;
            self.nGenerators = nG;
            self.C = zeros(1, nG*2);
            self.n = 1;
            self.p = [1];
            self.M = M;
        end

        function t = table(self)
            t = replab.str.Table(self.C);
            invNames = cellfun(@(x) ['inv(', x, ')'], self.generatorNames, 'UniformOutput', false);
            t.addColumnNames([self.generatorNames, invNames]);
            t.addRowNames(num2cell(1:size(self.C, 1)));
            t.setRowSep(1, '-');
            t.setColSep(1, ' | ');
        end

        function invind = generatorInverse(self, ind)
        % Returns the index of the inverse of the given generator index
            nG = self.nGenerators;
            if i <= nG
                invind = ind + nG;
            else
                invind = ind - nG;
            end
        end

        function define(self, alpha, x)
        % Defines a new coset
        %
        % DEFINE, Holt p. 153
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
            letters = replab.fp.parseLetters(w, self.generatorNames);
            s = sign(letters);
            ngens = self.nGenerators;
            f = find(letters < 0);
            letters(f) = arrayfun(@(x) abs(x) + ngens, letters(f));
            r = length(letters);
            f = alpha;
            i = 1;
            b = alpha;
            j = r;
            while true
                while i <= r && self.C(f, letters(i)) > 0
                    f = self.C(f, letters(i));
                    i = i + 1;
                end
                if i > r
                    if f ~= alpha
                        self.coincidence(f, alpha);
                    end
                    return
                end
                while j >= i && self.C(b, letters(j) + s(j)*ngens) > 0
                    b = self.C(b, letters(j) + s(j)*ngens);
                    j = j - 1;
                end
                if j < i
                    self.coincidence(f, b);
                    f = alpha;
                    i = 1;
                    b = alpha;
                    j = r;
                elseif j == i
                    self.C(f, letters(i)) = b;
                    self.C(b, letters(j) + s(j)*ngens) = f;
                    return
                else
                    self.define(f, letters(i));
                end
            end
        end

        function [q, l] = merge(self, k, lambda, q, l)
        % Removes the smaller of k and lambda from active cosets
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

        function  coincidence(self, alpha, beta)
        % Deals with the same coset being represented with different values
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
