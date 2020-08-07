classdef AutomorphismSearch < replab.bsgs.Backtrack
% Computes the automorphism group of a graph

    properties
        graph
    end

    methods

        function self = AutomorphismSearch(graph, knownSubgroup, debug)
            if nargin < 3 || isempty(debug)
                debug = false;
            end
            if nargin < 2 || isempty(knownSubgroup)
                knownSubgroup = [];
            end
            group = replab.S(graph.nVertices);
            base = group.lexChain.base;
            self@replab.bsgs.Backtrack(group, base, knownSubgroup, knownSubgroup, debug);
            self.graph = graph;
        end

        function ok = test(self, l, prev, ul)
            % TODO: refine this condition
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
