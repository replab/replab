classdef CosetTable
% Generates a coset table from the generators and relators of a finitely presented group and subset
%
% Based on Todd-Coxeter algorithm for coset enumeration from:
% Holt, Derek. “Coset Enumeration.” Handbook of Computational Group Theory,
% Chapman & Hall/CRC, 2004, pp. 149–198
%
%
% Returns:
%   table (`replab.str.Table`): coset table
%
% Example:
%   >>> replab.fp.cosetEnumeration({'x','y'}, {'x^2','y^3','(x*y)^3'}, {'x*y'})
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
            for i = 1:length(y)
                w = y{i};
                [C, n, p] = scanAndFill(C, n, p, M, 1, w, generators);
            end
        end

    end

    methods

        function self = CosetTable(generatorNames, M)
            nG = length(generators);
            self.generatorNames = generatorNames;
            self.nGenerators = nG;
            self.C = zeros(1, nG*2);
            self.n = 1;
            self.p = [1];
            self.M = M;
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
            if n == M
                error('Error: exceeded maximum allowed number of cosets')
            end
            nG = self.nGenerators;
            n = n + 1;
            beta = n;
            p(beta) = beta;
            C(alpha, x) = beta;
            xinv = self.generatorInverse(x);
            C(beta, xinv) = alpha;
        end



    end

end

function table = cosetEnumeration(generators, relators, y)

    alpha = 1; % since p(1) is always 1
    while alpha <= n
        if p(alpha) ~= alpha
            alpha = alpha + 1;
            continue
        end
        for i = 1:length(relators)
            w = relators{i};
            [C, n, p] = scanAndFill(C, n, p, M, alpha, w, generators);
            if p(alpha) < alpha
                break
            end
        end
        if p(alpha) < alpha
            alpha = alpha + 1;
            continue
        end
        for x = 1:ngens * 2
            if C(alpha, x) < 1
                [C, n, p] = define(C, M, n, p, alpha, x);
            end
        end
        alpha = alpha + 1;
    end
    f = find(p - (1:n) == 0);
    C = C(1:length(f), :);
    table = replab.str.Table(C);
    invNames = cellfun(@(x) ['inv(', x, ')'], generators, 'UniformOutput', false);
    table.addColumnNames([generators, invNames])
    table.addRowNames(num2cell(1:length(f)))
    table.setRowSep(1, '-')
    table.setColSep(1, ' | ')
end

function [C, n, p] = scanAndFill(C, n, p, M, alpha, w, gens)
% Scans through a relator and fills in table
    letters = replab.fp.parseLetters(w, gens);
    s = sign(letters);
    ngens = length(gens);
    f = find(letters < 0);
    letters(f) = arrayfun(@(x) abs(x) + ngens, letters(f));
    r = length(letters);
    f = alpha;
    i = 1;
    b = alpha;
    j = r;
    while true
        while i <= r && C(f, letters(i)) > 0
            f = C(f, letters(i));
            i = i + 1;
        end
        if i > r
            if f ~= alpha
                [C, p] = coincidence(C, p, f, alpha);
            end
            return
        end
        while j >= i && C(b, letters(j) + s(j)*ngens) > 0
            b = C(b, letters(j) + s(j)*ngens);
            j = j - 1;
        end
        if j < i
            [C, p] = coincidence(C, p, f, b);
            f = alpha;
            i = 1;
            b = alpha;
            j = r;
        elseif j == i
            C(f, letters(i)) = b;
            C(b, letters(j) + s(j)*ngens) = f;
            return
        else
            [C, n, p] = define(C, M, n, p, f, letters(i));
        end
    end
end


function [lambda, p] = rep(k, p)
% replaces a coset in p
    lambda = k;
    rho = p(lambda);
    while rho ~= lambda
        lambda = rho;
        rho = p(lambda);
    end
    mu = k;
    rho = p(mu);
    while rho ~= lambda
        p(mu) = lambda;
        mu = rho;
        rho = p(mu);
    end
end

function [p, q, l] = merge(k, lambda, p, q, l)
% Removes the smaller of k and lambda from active cosets
    lam1 = rep(k, p);
    lam2 = rep(lambda, p);
    if lam1 ~= lam2
        mu = min([lam1, lam2]);
        v = max([lam1, lam2]);
        p(v) = mu;
        l = l + 1;
        q(l) = v;
    end
end

function [C, p] = coincidence(C, p, alpha, beta)
% Deals with the same coset being represented with different values
    dim = size(C);
    ngens = dim(2) / 2; % unless we start removing self-inverses
    l = 0;
    q = [];
    [p, q, l] = merge(alpha, beta, p, q, l);
    i = 1;
    while i <= l
        g = q(i);
        i = i + 1;
        gx = C(g, :);
        for j = 1:dim(2)
            if gx(j) > 0
                delta = gx(j);
                if j <= ngens
                    jinv = j + ngens;
                else
                    jinv = j - ngens;
                end
                C(delta, jinv) = 0;
                [mu, p] = rep(g, p);
                [nu, p] = rep(delta, p);
                if C(mu, j) > 0
                    [p, q, l] = merge(nu, C(mu, j), p, q, l);
                elseif C(nu, jinv) > 0
                    [p, q, l] = merge(mu, C(nu, jinv), p, q, l);
                else
                    C(mu, j) = nu;
                    C(nu, jinv) = mu;
                end
            end
        end
    end
end
