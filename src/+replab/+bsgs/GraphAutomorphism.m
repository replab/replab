classdef GraphAutomorphism < replab.bsgs.Backtrack
% Computes the automorphism group of a graph

    properties
        graph
        candidate % permutation under consideration
        lastL     % last level considered
        matchesU  % potential matches on the left
        matchesV  % potential matches on the right
        Phi       % known permutation images
    end

    methods

        function self = GraphAutomorphism(graph, knownSubgroup, debug)
            if nargin < 3 || isempty(debug)
                debug = false;
            end
            if nargin < 2 || isempty(knownSubgroup)
                knownSubgroup = [];
            end
            group = replab.S(graph.nVertices);
            base = 1:graph.nVertices;
            self@replab.bsgs.Backtrack(group, base, knownSubgroup, knownSubgroup, debug);
            self.graph = graph;

            self.candidate = zeros();
            self.lastL = 0;
            self.matchesU = {};
            self.matchesV = {};
        end

        function ok = test(self, l, prev, ul)
            % TODO: optimize this condition

%            disp(num2str(l));
%            disp(num2str([prev; ul]));
            
            % Initialize variables
            ok = false;
            tolerance = 1e-10; % Since we allow false positives, it's ok if this number is significantly bigger than eps
            nVertices = self.graph.nVertices;
            Kt = self.graph.heatKernel;
            nbT = size(Kt,1); % number of considered timings
            p = 1:l;
            pTilde = candidate(1:l);
            
            % We only update the last possibly modified item
            self.candidate(l) = prev(ul(l));

            if l == 1
                % Initialize the potential matches
                matchesU = cell(1, nVertices);
                matchesV = cell(1, nVertices);
                for i = 1:nVertices
                    matchesU{i} = 1:nVertices;
                    matchesV{i} = 1:nVertices;
                end
            end
            if self.lastL > l
                % We erase the
                matchesU(matchesU > l) = l;
                matchesV(matchesV > l) = l;
            end
            
            
            self.lastL = l;

            return;
            
            

            % Here is the candidate permutation
            candidate = prev(ul);

            % Initialize variables
            tolerance = 1e-10;
            nVertices = self.graph.nVertices;
            nbT = size(self.graph.heatKernel,1); % number of considered timings
            ok = false;
            p = 1:l;
            pTilde = candidate(1:l);
            Kt = self.graph.heatKernel;

            % The partial
            Phi = zeros(1, nVertices);
            Phi(1:l) = candidate(1:l);
            
            % We initialize the possible matches of all vertices
            % When given just the Kernel at different times, we initialize Phi, and
            % scan all possible p and pTilde
            matchesU = cell(1, nVertices);
            matchesV = cell(1, nVertices);
            for i = 1:nVertices
                matchesU{i} = 1:nVertices;
                matchesV{i} = 1:nVertices;
            end
            
            % It seems that we should still check all case, so:
            % We compute the vectors corresponding to p and pTilde for all vertices
            Sku = reshape(Kt(:,p(end),:), nbT, nVertices);
            Skv = reshape(Kt(:,pTilde(end),:), nbT, nVertices);

            % We check that previous matches still match...
            for u = 1:nVertices
                if (Phi(u) ~= 0) && (sum(abs(Sku(:,u) - Skv(:,Phi(u)))) >= tolerance)
                    % Some assignment done previously doesn't work
                    ok = false;
                    return;
                end
            end
            
            % We refine the matches
            %candidatesV = setdiff(1:length(Phi), Phi); % too slow...
            candidatesV = 1:nVertices;
            candidatesV(Phi(Phi~=0)) = 0;
            candidatesV = candidatesV(candidatesV~=0);
            candidatesU = find(Phi == 0);

            nMatchesU = zeros(1,nVertices);
            for u = candidatesU
                differenceU = repmat(Sku(:,u), [1, length(matchesU{u})]) - Skv(:,matchesU{u});
                matchU = sum(abs(differenceU)) < tolerance;
                matchesU{u} = matchesU{u}(find(matchU));
                nMatchesU(u) = sum(matchU);
                if nMatchesU(u) == 0
                    % No match possible
                    ok = false;
                    return;
                end
            end
            
            nMatchesV = zeros(1,nVertices);
            for v = candidatesV
                differenceV = Sku(:,matchesV{v}) - repmat(Skv(:,v), [1, length(matchesV{v})]);
                matchV = sum(abs(differenceV)) < 1e-10;
                matchesV{v} = matchesV{v}(find(matchV));
                nMatchesV(v) = sum(matchV);
                if nMatchesV(v) == 0
                    % No match possible
                    allPhi = [];
                    success = false;
                    return;
                end
            end

            % matches should be unique in both directions
            for u = candidatesU
                if (nMatchesU(u) == 1) && (nMatchesV(matchesU{u}(1)) == 1)
                    Phi(u) = matchesU{u};
                end
            end
            
            % If everything is assigned, we exit successfully
            if min(Phi) > 0
                % If we made an inconsistant deduction, we abort this branch and produce a
                % warning
                if ~isequal(sort(Phi),[1:nVertices])
                    warning(['Incorrect permutation: ', num2str(Phi), ' found with the pair [', num2str(p), '], [', num2str(pTilde), '].']);
                    ok = false;
                    return;
                end

            %    disp(['With the pair [', num2str(p), '], [', num2str(pTilde), '], found permutation: ', num2str(Phi)]);
                ok = true;
                return;
            end

            % If a match is missing, we abort this branch
            if min(nMatchesU(candidatesU)) == 0
                % Most likely, p cannot be exchanged with pTilde...
            %    disp(['The following pair is not possible: [', num2str(p), '], [', num2str(pTilde), '].']);
                ok = false;
                return;
            end
            
            % If we reached here, we cannot exclude this partial assignment
            ok = true;
            
%             % If some freedom is left, we explore the different possible branches.
%             % Note that we can assume that elements in p increase, so no need to try
%             % elements smaller than those already in p.
%             toTry = p(end)+find(nMatchesU(p(end)+1:end) >= 2);
% 
%             allPhi = {};
%             co = 0;
%             for newP = toTry
%                 for newPTilde = matchesU{newP}
%                     [newPhi, success] = graphAutomorphisms_1(Kt, [p newP], [pTilde newPTilde], Phi, matchesU, matchesV);
%                     if success
%                         co = co + 1;
%                         allPhi{co} = newPhi;
%                     end
%                 end
%             end
%             allPhi = cat(1, allPhi{:});
%             allPhi = unique(allPhi, 'rows');
        end

        function ok = prop(self, g)
            ok = false;
            if numel(self.graph.colors) > 1
                ok = isequal(self.graph.colors(g), self.graph.colors);
            end
            if ~ok
                adj = self.graph.adjacencyMatrix;
                ok = isequal(adj, adj(g, g));
            end
        end

    end

end
