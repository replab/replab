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
            sub = group.lexChain;
            beta = 1;
            while beta <= n
                candidates1 = zeros(0, n);
                minImg = n + 1;
                [stab1 orbit1 iOrbit1 U1 Uinv1] = sub.stabilizer(beta, true);
                for j = 1:length(orbit1)
                    b = orbit1(j);
                    u1 = U1(:,j)';
                    uinv1 = Uinv1(:,j)';
                    for i = 1:size(candidates, 1)
                        h = candidates(i,:);
                        a = uinv1(h(b));
                        if i == 1 || candidates(i,b) ~= candidates(i-1,b)
                            [stab1a orbit1a iOrbit1a U1a Uinv1a] = stab1.stabilizer(a, true);
                        end
                        [m, ind] = min(orbit1a);
                        if m < minImg
                            candidates1 = zeros(0, n);
                            minImg = m;
                        end
                        if m <= minImg
                            u1a = U1a(:,ind)';
                            uinv1a = Uinv1a(:,ind)';
                            % ``rest u1a uinv1 h u1 uinv1a rest``
                            h1 = u1a(uinv1(h(u1(uinv1a))));
                            candidates1(end+1,:) = h1;
                        end
                    end
                end
                sub = stab1.stabilizer(minImg);
                candidates = unique(candidates1, 'rows');
                beta = beta + 1;
            end
            assert(size(candidates, 1) == 1);
            h1 = candidates(1,:);
            g = replab.bsgs.LeftConjugation(group, h, h1).find;
            assert(all(group.leftConjugate(g, h) == h1));
        end

    end

end
