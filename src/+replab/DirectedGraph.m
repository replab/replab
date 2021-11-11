classdef DirectedGraph < replab.graph.Graph
% Describes an immutable directed graphs

    methods (Access = protected)

        function self = DirectedGraph(nVertices, edges, colors, weights)
        % Construct a directed graph
        %
        % Do not use this function directly, rather use another
        % constructor such as `.fromAdjacencyMatrix` .
            self@replab.graph.Graph(nVertices, edges, colors, weights);
        end

    end

    methods (Static) % Constructors

        function self = fromAdjacencyMatrix(adj, colors)
        % Constructs a directed graph from an adjacency matrix
        %
        % The element (i,j) of an adjacency matrix is 1 iff vertex i is
        % linked to vertex j. Alternatively, the value of the matrix
        % element is the weight associated to this edge.
        %
        % Args:
        %   adj (double (\*,\*)): the adjacency matrix
        %   colors (double (\*,1), optional): coloring of the vertices
        %
        % Returns:
        %   graph (`.DirectedGraph`)
        %
        % Example:
        %   >>> replab.DirectedGraph.fromAdjacencyMatrix([0 1 1; 1 0 1; 1 1 0])
        %     Directed graph with 3 vertices and 6 edges
        %     edges: [2, 1; 3, 1; 1, 2; 3, 2; 1, 3; 2, 3]

            assert(size(adj,1) == size(adj,2), 'Adjacency matrix should be square');

            [edges, nVertices, weights] = replab.graph.adj2edge(adj);

            if nargin < 2
                colors = 0;
            end

            self = replab.DirectedGraph(nVertices, edges, colors, weights);
        end

        function self = fromEdges(edges, nVertices, colors, weights)
        % Constructs a directed graph from a liste of edges
        %
        % Duplicated edges are merged. Their weights are added up iff a
        % weight was defined individually for each edge.
        %
        % Args:
        %     edges (integer (\*,2)): array of vertices linked by an edge
        %     nVertices (integer, optional): number of vertices
        %     colors (double (\*,1), optional): coloring of the vertices
        %     weights (double (\*,1), optional): weight associated to each edge
        %
        % Returns:
        %   graph (`.DirectedGraph`)
        %
        % Example:
        %   >>> replab.DirectedGraph.fromEdges([1 2; 2 3; 3 1], 3)
        %     Directed graph with 3 vertices and 3 edges
        %     edges: [1, 2; 2, 3; 3, 1]

            if nargin < 2
                if isempty(edges)
                    nVertices = 0;
                else
                    nVertices = max(max(edges));
                end
            end

            if (nargin < 3)
                colors = [];
            end

            if (nargin < 4)
                weights = [];
            end

            % Remove duplicated edges
            if (numel(weights) > 1) && (numel(weights) == size(edges,1))
                % We add up the weights associated to identical edges
                [edges, IA, IC] = unique(edges, 'rows', 'first');
                weights = full(sparse(IC, 1, weights));
            else
                % The edge weight remains uniform
                edges = unique(edges, 'rows');
            end

            self = replab.DirectedGraph(nVertices, edges, colors, weights);
        end

        function self = fromBiadjacencyMatrix(biadj)
        % Constructs a directed graph from a biadjacency matrix
        %
        % A bipartite directed graph is one in which vertices can
        % be divided into two sets, such that edges only connect
        % vertices belonging to different sets. Moreover all edges are
        % directed from the first set to the second one.
        %
        % The element (i,j) of an biadjacency matrix is 1 for an
        % directed bipartite graph iff vertex i of the first set of
        % verices is linked to vertex j in the second set of vertices.
        % Alternatively, the value of the matrix element is the weight
        % associated to this edge.
        %
        % Args:
        %   biadj(double (\*,\*)): array of vertices linked by an edge
        %
        % Returns:
        %   graph (`.DirectedGraph`)
        %
        % Example:
        %   >>> replab.DirectedGraph.fromBiadjacencyMatrix([1 0 1; 1 1 0])
        %     Directed graph with 5 vertices and 4 edges
        %     edges: [1, 3; 2, 3; 2, 4; 1, 5]

            % Construct the associated full adjacency matrix
            adj = [zeros(size(biadj,1)*[1 1]), biadj;
                   zeros(size(biadj,2), sum(size(biadj)))];

            self = replab.DirectedGraph.fromAdjacencyMatrix(adj);
        end

        function self = fromUndirectedGraph(graph)
        % Define an directed graph from an undirected one
        %
        % This function returns an directed graph with the same
        % connectivity as the input undirected graph.
        %
        % Args:
        %   graph (`.UndirectedGraph`)
        %
        % Returns:
        %   graph (`.DirectedGraph`)

            assert(isa(graph, 'replab.UndirectedGraph'), 'Input is not an undirected graph');

            % We duplicate edges
            edges = graph.edges;
            sel = (diff(edges,1,2) ~= 0); % only duplicate links to distinct vertices
            edges = [edges; fliplr(edges(sel,:))];
            weights = graph.weights;
            if numel(weights) ~= 1
                weights = [weights; weights(sel)];
            end

            self = replab.DirectedGraph.fromEdges(edges, graph.nVertices, graph.colors, weights);
        end

        function self = fromFile(fileName)
        % Loads a graph from a Bliss file
        %
        % Following the specifications from
        % http://www.tcs.hut.fi/Software/bliss/fileformat.shtml
        %
        % Args:
        %   fileName (charstring) : Name of the file to be loaded
        %
        % Returns:
        %   graph (`.UndirectedGraph`)

            [nVertices, edges, colors] = replab. graph.parseBlissFile(fileName);

            self = replab.DirectedGraph.fromEdges(edges, nVertices, colors);
        end

    end

    methods

        function s = headerStr(self)
        % Header string representing the object
        %
        % Returns:
        %     charstring: description of the object

            s = sprintf('Directed graph with %d vertices and %d edges', self.nVertices, size(self.edges,1));
        end

    end

    methods % Methods

        function graph = toUndirectedGraph(self)
        % Removes the graph directionality
        %
        % This function returns an undirected graph with the same
        % connectivity as the current graph.
        %
        % Returns:
        %   graph (`.DirectedGraph`)

            graph = replab.UndirectedGraph.fromDirectedGraph(self);
        end

        function adj = computeAdjacencyMatrix(self)
        % Computes the adjacency matrix

            adj = replab.graph.edge2adj(self.edges, self.nVertices, self.weights);
            adj = full(adj);
        end

        function deg = degrees(self)
        % Returns the degrees of all vertices
        %
        % Returns:
        %   deg (double (1,\*)): list of degrees
        %
        % Example:
        %   >>> replab.DirectedGraph.fromEdges([1 3]).degrees
        %     1     0     1

            deg = sum(self.adjacencyMatrix()) + sum(self.adjacencyMatrix().');
        end

        function deg = degree(self, v)
        % Returns the degrees of vertex v
        %
        % Args:
        %   v (integer) : vertex number
        %
        % Returns:
        %   deg (double (1,\*)): list of degrees
        %
        % Example:
        %   >>> graph = replab.DirectedGraph.fromEdges([1 3; 1 4]);
        %   >>> graph.degree(1)
        %     2

            assert(all(v >= 0) && all(v <= self.nVertices) && isequal(v, round(v)), ...
                ['No vertex number ', num2str(v)]);

            adj = self.adjacencyMatrix();
            deg = sum(adj(v,:)) + sum(adj(:,v).');
        end

        function deg2 = secondOrderDegree(self, v)
        % Returns the number of vertices at a distance 2 of all vertices
        %
        % Args:
        %   v (integer) : vertex number
        %
        % Returns:
        %   deg (double (1,\*)): list of degrees
        %
        % Example:
        %   >>> replab.DirectedGraph.fromEdges([1 3]).secondOrderDegrees
        %     1     0     1

            assert((numel(v) == 1) && (v >= 0) && (v <= self.nVertices) && isequal(v, round(v)), ...
                ['No vertex number ', num2str(v)]);

            adj = self.adjacencyMatrix();
            sel = (adj(v,:) ~= 0) | (adj(:,v).' ~= 0);
            deg2 = sum(sum(adj(sel,:) + adj(:,sel).') ~=0);
        end

        function L = computeLaplacian(self)
        % Computes the graph Laplacian
        %
        % The laplacian should be positive semi-definite, so we define it
        % from the equivalent un-directed graph

%             % Here we compute a hermitian Laplacian
%             M = self.adjacencyMatrix;
%             selSym = (M == M.');
%             MSym = M;
%             MSym(~selSym) = 0;
%             MAsym = M;
%             MAsym(selSym) = 0;
%             MH = MSym + 1i*(MAsym - MAsym.');
%
%             D = diag(sum(M + M.'));
%
%             colors = abs(self.colors);
%             if numel(colors) == 1
%                 colors = colors*ones(1,self.nVertices);
%             end
%
%             L = D - MH + diag(colors);
%
%             return;

            M = self.toUndirectedGraph.adjacencyMatrix;
            D = diag(sum(M));

            colors = abs(self.colors);
            if numel(colors) == 1
                colors = colors*ones(1,self.nVertices);
            end

            L = D - M + diag(colors);
        end

    end

end
