classdef DirectedGraph < replab.graph.Graph
% Describes an immutable directed graphs

    methods (Access = protected)
        
        function self = DirectedGraph(nVertices, edges, colors, weights)
        % Construct a directed graph
        %
        % Do not use this function direclty, rather use another
        % constructor such as ``.fromBlocks``.
        %
        % See also:
        %   `.fromBlocks`
        %   `.check`

            self@replab.graph.Graph(nVertices, edges, colors, weights)
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
        %   >>> replab.DirectedGraph.fromAdjacencyMatrix([0 1 1; 1 0 1; 1 1 0]);
            
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
        %   >>> replab.DirectedGraph.fromEdges([1 2; 2 3; 3 1], 3);

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
                [edges, IA, IC] = unique(edges, 'rows');
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
        %   >>> replab.DirectedGraph.fromBiadjacencyMatrix([1 0 1; 1 1 0]);

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
            % Args:
            %   graph (`.DirectedGraph`)
            %
            % Returns:
            %   graph (`.UndirectedGraph`)
            
            graph = replab.UndirectedGraph.fromDirectedGraph(self);
        end
        
        function adj = computeAdjacencyMatrix(self)
        % Computes the adjacency matrix

            adj = replab.graph.edge2adj(self.edges, self.nVertices, self.weights);
            adj = full(adj);
        end
        
        function L = computeLaplacian(self)
        % Computes the graph Laplacian
        %
        % The laplacian should be positive semi-definite, so we define it
        % from the equivalent un-directed graph
        
            M = self.toUndirectedGraph.adjacencyMatrix();
            D = diag(sum(M));
            
            colors = abs(self.colors);
            if numel(colors) == 1
                colors = colors*ones(1,self.nVertices);
            end
            
            L = D - M + diag(colors);
        end
        
    end
    
end
