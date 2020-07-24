classdef ConjugacyClasses
% Computations involving conjugacy classes

    methods (Static)

        function [h1 g] = representative(group, h)
        % Returns the lexicographically minimal representative of a conjugacy class
        %
        % We find the permutation ``g`` such that ``h1 = group.leftConjugate(g, h)`` is lex-minimal.
        %
        % For that, we first look for ``g = u1 g1`` with ``u1`` a transversal of the cosets
        % for the first base point ``beta = 1`` of the stabilizer chain of ``group``.
        %
        % Now, we look at ``g`` such that ``gInv(h(g(beta)))`` is minimal. We minimize
        % ``gInv1(uInv1(h(u1(g1(beta))))) = gInv1(uInv1(h(b)))`` by examining all
        % possibilities for ``b``. We write ``a = uInv1(h(b))`` and examine the
        % possibilities for ``gInv1(a)`` by considering the orbit of ``a`` in
        % the stabilizer of ``group`` by ``beta``. We write ``gInv1 = ua g1a``
        % where ``ua`` is a transversal of the cosets of this second stabilizer, and ``g1a``
        % is an element of ``group`` that stabilizes both ``1`` and ``a``. This provides
        % a candidate for ``g`` at the first level.
        %
        % At the first level, we collect all candidates ``h1 = gInv h g``, and work on the
        % second base point ``beta = 2``, now with ``group`` stabilizing both ``beta_1 = 1`` and ``gInv(h(g(1)))``.
        %
        % Args:
        %   group (`+replab.+bsgs.Permutatoin`): Group to explore
        %   h (permutation): Group element
        %
        % Returns
        % -------
        % h1:
        %   permutation: Lexicographically minimal element in the conjugacy class of ``h``
        % g:
        %   permutation: Permutation that left conjugates ``h`` to ``h1``
            n = group.domainSize;
            candidates = h; % candidates as rows
            conjugates = 1:n; % conjugating elements as rows
            sub = group.lexChain;
            beta = 1;
            while beta <= n
                candidates1 = zeros(0, n);
                conjugates1 = zeros(0, n);
                minImg = n + 1;
                stab2 = []; % corresponds to the result of sub.stabilizer(beta).stabilizer(minImg)
                orbit2 = [];
                iOrbit2 = [];
                U2 = [];
                Uinv2 = [];
                [stab1 orbit1 iOrbit1 U1 Uinv1] = sub.stabilizer(beta, true);
                for j = 1:length(orbit1)
                    b = orbit1(j);
                    u1 = U1(:,j)';
                    uinv1 = Uinv1(:,j)';
                    for i = 1:size(candidates, 1)
                        c = candidates(i,:);
                        ci = conjugates(i,:);
                        cc = group.leftConjugate(ci, h);
                        assert(all(c == cc));
                        a = uinv1(c(b));
                        if i == 1 || candidates(i,b) ~= candidates(i-1,b)
                            orbit2 = stab1.orbitUnderG(1, a);
                        end
                        m = min(orbit2);
                        if m <= minImg
                            if m < minImg
                                minImg = m;
                                candidates1 = zeros(0, n);
                                conjugates1 = zeros(0, n);
                                [stab2 orbit2 iOrbit2 U2 Uinv2] = stab1.stabilizer(minImg, true);
                            end
                            ind = iOrbit2(a);
                            u2 = U2(:,ind)';
                            uinv2 = Uinv2(:,ind)';
                            c1 = uinv2(uinv1(c(u1(u2))));
                            candidates1(end+1,:) = c1;
                            conjugates1(end+1,:) = uinv2(uinv1(ci));
                        end
                    end
                end
                sub = stab2;
                [candidates, I] = unique(candidates1, 'rows');
                conjugates = conjugates1(I,:);
                beta = beta + 1;
            end
            assert(size(candidates, 1) == 1);
            h1 = candidates(1,:);
            g = conjugates1(1,:);
            assert(all(group.leftConjugate(g, h) == h1));
        end

    end

end
