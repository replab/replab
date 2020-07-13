classdef Cosets

    methods (Static)

        function T = leftTransversalAsMatrix(group, subgroup)
        % Returns, as a matrix, the transversal of the left cosets of a subgroup in a group
        %
        % Args:
        %   group (`+replab.+bsgs.Chain`): Chain with monotonic base
        %   subgroup (`+replab.+bsgs.Chain`): Chain with monotonic base
        %
        % Returns:
        %   integer(\*,\*): Transversal matrix with permutations as columns
            Tinv = replab.bsgs.rightTransversalMatrix(group, subgroup);
            T = zeros(size(Tinv));
            n = group.n;
            g = zeros(1, n);
            for i = 1:size(Tinv, 1)
                g(T(:,i)) = 1:n; % invert
                T(:,i) = replab.bsgs.leftRepresentative(subgroup, g);
            end
        end

        function T = rightTransversalMatrix(group, subgroup)
        % Returns, as a matrix, the transversal of the right cosets of a subgroup in a group
        %
        % Args:
        %   group (`+replab.+bsgs.Chain`): Chain with monotonic base
        %   subgroup (`+replab.+bsgs.Chain`): Chain with monotonic base
        %
        % Returns:
        %   integer(\*,\*): Transversal matrix with permutations as columns
            group = self.groupChain;
            subgroup = self.subgroupChain;
            n = group.n;
            L = group.length;
            stackG = zeros(n, L+1); % current group element
            stackG(:,1) = 1:n;
            stackS = cell(1, L+1);
            stackS{1} = subgroup; % subgroup chain
            stackB = cell(1, L+1);
            stackB{1} = sort(group.Delta{1}); % remaining points in orbit to check
            stackM = cell(1, L+1);
            stackM{1} = replab.bsgs.minimalMaskInOrbit(n, stackS{1}.strongGeneratorsForLevel(1));
            l = 1;
            T = [];
            while l >= 1
                g = stackG(:,l)';
                while 1
                    B = stackB{l};
                    if isempty(B)
                        l = l - 1;
                        break
                    end
                    b = B(1);
                    stackB{l} = B(2:end);
                    bg = g(b);
                    if stackM{l}(bg)
                        if l == L
                            T(:,end+1) = g(group.u(l, b));
                        else
                            stackS{l+1} = stackS{l}.stabilizer(bg);
                            gnext = g(group.u(l, b));
                            stackG(:,l+1) = gnext; % compose transversal element
                            orbit = sort(group.Delta{l+1}); % TODO: unnecessary sort
                            [~, I] = sort(gnext(orbit));
                            stackB{l+1} = orbit(I);
                            stackM{l+1} = replab.bsgs.minimalMaskInOrbit(n, stackS{l+1}.strongGeneratorsForLevel(1));
                            l = l + 1;
                            break
                        end
                    end
                end
            end
        end

        function t = leftRepresentative(subgroup, g)
        % Computes the representative of the left coset in which ``g`` belongs
        %
        % Args:
        %   subgroup (`+replab.+bsgs.Chain`): Chain with monotonic base
        %   g (permutation): Permutation to find a representative of
        %
        % Returns:
        %   permutation: The coset representative
            L = subgroup.length;
            n = subgroup.n;
            h = 1:n;
            % we are looking for an element of the form g h1 h2 ... hk
            % so finding the representative is easy, we just look for the
            % element hi so that ((g h1 h2 ... hi-1) hi)(i) is minimal
            %
            % at all times, h = h1 h2 ... h{i-1}
            for l = 1:L
                orbit = subgroup.Delta{l};
                [~, i] = min(g(h(orbit)));
                Ul = subgroup.U{l};
                h = h(Ul(:,i)');
            end
            t = g(h);
        end

        function t = rightRepresentative(subgroup, g)
        % Computes the representative of the right coset in which ``g`` belongs
        %
        % Args:
        %   subgroup (`+replab.+bsgs.Chain`): Chain with monotonic base
        %   g (permutation): Permutation to find a representative of
        %
        % Returns:
        %   permutation: The coset representative
            sub = subgroup.mutableCopy;
            n = subgroup.domainSize;
            % we build an element of the form h1 ... hn g
            % we iterate over the sequence
            % h1 g <- find h1
            % h1 h2 g <- find h2
            % (h1 h2) h3 g <- find h3
            % ...
            % At the i-th step, the subgroup chain fixes the points
            % g(1) ... g(i-1), because these images are minimal under
            % h1 ... h_{i-1}
            % Then we compute the stabilizer and orbits of the subgroup chain
            % with respect to g(i). The orbit gives the image of hi(g(i)) for
            % possible transversal elements hi. We then pick the one that provides
            % the minimal image ((h1 ... h_{i-1}) h1 g)(i)
            h = 1:n; % subgroup element, = h1 h2 ... h{i-1}
            for i = 1:n
                beta = g(i);
                % verifies if the stabilizer chain stabilizes beta already
                if all(sub.S(beta,:) == beta)
                    % do nothing
                else
                    [sub, orbit, ~, U, ~] = sub.stabilizer(beta);
                    [~, ind] = min(h(orbit));
                    h = h(U(:,ind)');
                end
            end
            t = h(g);
        end

    end

end
