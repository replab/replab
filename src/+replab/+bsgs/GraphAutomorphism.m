classdef GraphAutomorphism < replab.bsgs.Backtrack
% Computes the automorphism group of a graph

    properties
        graph
        matchesU  % potential matches on the left
        matchesV  % potential matches on the right
        Phi       % stores the images which are already inferred
    end

    methods

        function self = GraphAutomorphism(graph, knownSubgroup, debug)
            if nargin < 3 || isempty(debug)
                debug = false;
            end
            if nargin < 2 || isempty(knownSubgroup)
                knownSubgroup = [];
            end
            
            % We compute a group 
            group = replab.S(graph.nVertices);
            
            % Try to make the group smaller by taking into account some
            % vertices invariants
            if graph.nVertices > 10
                invariant1 = graph.colors;
                invariant2 = graph.degrees;
                invariant3 = graph.secondOrderDegrees;

                if numel(invariant1) > 1
                    group = group.vectorStabilizer(invariant1(:).');
                end
                if length(unique(invariant2)) > 1
                    group = group.vectorStabilizer(invariant2(:).');
                end
                if length(unique(invariant3)) > 1
                    group = group.vectorStabilizer(invariant3(:).');
                end
            end
            
            base = 1:graph.nVertices;
            self@replab.bsgs.Backtrack(group, base, knownSubgroup, knownSubgroup, debug);
            self.graph = graph;

            % Initialize the potential matches at level 1
            self.matchesU = cell(graph.nVertices,graph.nVertices);
            self.matchesV = cell(graph.nVertices,graph.nVertices);
            for i = 1:graph.nVertices
                self.matchesU{1,i} = 1:graph.nVertices;
                self.matchesV{1,i} = 1:graph.nVertices;
            end
            
            % Initialize the known permutation at level 1
            self.Phi = cell(1,graph.nVertices);
            self.Phi{1} = zeros(1,graph.nVertices);
        end

        function ok = test(self, l, prev, ul)
            % Here is the candidate permutation
            candidate = prev(ul);

            % We quickly check whether the we can answer directly
            Phi = self.Phi{l};
            
            % Quickly check if we already excluded this case
            if (Phi(l) ~= 0) && (Phi(l) ~= candidate(l))
                ok = false;
                return;
            end
            
            % Check if anything still needs to be computed
            if min(Phi) > 0
                self.Phi{l+1} = Phi;
                ok = true;
                return;
            end            

            
            %% Ok, we need to compute seomthing to answer the question
            
            % Initialize variables
            ok = false;
            tolerance = 1e-10; % Since we allow false positives, it's ok if this number is significantly bigger than eps
            nVertices = self.graph.nVertices;
            Kt = self.graph.heatKernel;
            nbT = size(Kt,1); % number of considered timings
            p = 1:l;
            pTilde = candidate(1:l);
            
            % Loads the previous coarse grained matches
            matchesU = self.matchesU(l,:);
            matchesV = self.matchesV(l,:);
            
            % We compute the vectors corresponding to p and pTilde for all vertices
            Sku = reshape(Kt(:,p(end),:), nbT, nVertices);
            Skv = reshape(Kt(:,pTilde(end),:), nbT, nVertices);
            
            % Now we refine the matches
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
                    ok = false;
                    return;
                end
            end
            
            % We make all possible inferences.
            % Matches should be unique in both directions
            for u = candidatesU
                if (nMatchesU(u) == 1) && (nMatchesV(matchesU{u}(1)) == 1)
                    Phi(u) = matchesU{u};
                end
            end
            
            % We save the result obtained for the next step
            self.matchesU(l+1,:) = matchesU;
            self.matchesV(l+1,:) = matchesV;
            self.Phi{l+1} = Phi;
            
            % If everything is assigned, we exit successfully
            if min(Phi) > 0
                ok = true;
                return;
            end

            % If a match is missing, we abort this branch
            if min(nMatchesU(candidatesU)) == 0
                ok = false;
                return;
            end
            
            % If we reached here, we cannot exclude this partial assignment
            ok = true;
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
