classdef UndirectedGraph < replab.DirectedGraph
% Describes an immutable undirected graph
%
% Vertices of the graph are numbered continuously from 1 to nVertices

    methods (Access = public)
        
        function self = UndirectedGraph(nVertices, edges, colors, weights)
        % Construct an undirected graph
        %
        % Do not use this function direclty, rather use another
        % constructor such as ``.fromBlocks``.
        %
        % See also:
        %   `.fromBlocks`
        %   `.check`

            self@replab.DirectedGraph(nVertices, edges, colors, weights);
        end

    end
    
    methods (Static) % Constructors

        function self = fromAdjacencyMatrix(adj, colors)
        % Constructs an undirected graph from an adjacency matrix
        %
        % The element (i,j) of an adjacency matrix is 1 only if vertex i is
        % connected to vertex j. Alternatively, the value of the matrix
        % element is the weight associated to this edge.
        %
        % Args:
        %   adj (double (\*,\*)): the adjacency matrix
        %   colors (double (\*,1), optional): coloring of the vertices
        %
        % Returns:
        %   graph (`.UndirectedGraph`)
        %
        % Example:
        %   >>> replab.UndirectedGraph.fromAdjacencyMatrix([0 1 1; 1 0 1; 1 1 0]);
            
            assert(size(adj,1) == size(adj,2), 'Adjacency matrix should be square');
            
            % We symmetrize the adjacency matrix
            adj = max(adj, adj.');
            
            % But keep only a short description
            adj = triu(adj);
            
            [edges, nVertices, weights] = replab.graph.adj2edge(adj);

            if nargin < 2
                colors = 0;
            end
            
            self = replab.UndirectedGraph(nVertices, edges, colors, weights);
        end
        
        function self = fromEdges(edges, nVertices, colors, weights)
        % Constructs a undirected graph from a liste of edges
        %
        % Duplicated edges are merged. Their weights are added up iff a
        % weights was defined individually for each edge.
        %
        % Args:
        %     edges (integer (\*,2)): array of vertices linked by an edge
        %     nVertices (integer, optional): number of vertices
        %     colors (double (\*,1), optional): coloring of the vertices
        %     weights (double (\*,1), optional): weight associated to each edge
        %
        % Returns:
        %   graph (`.UndirectedGraph`)
        %
        % Example:
        %   >>> replab.UndirectedGraph.fromEdges([1 2; 2 3; 3 1], 3);

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
            
            % Symmetrize the connections
            edges = sort(edges, 2);
            
            % Remove duplicated edges
            if (numel(weights) > 1) && (numel(weights) == size(edges,1))
                % We add up the weights associated to identical edges
                [edges, IA, IC] = unique(edges, 'rows');
                weights = full(sparse(IC, 1, weights));
            else
                % The edge weight remains uniform
                edges = unique(edges, 'rows');
            end
            
            self = replab.UndirectedGraph(nVertices, edges, colors, weights);
        end

        function self = fromBiadjacencyMatrix(biadj)
        % Constructs an undirected graph from a biadjacency matrix
        %
        % A bipartite undirected graph is one in which vertices can
        % be divided into two sets, such that edges only connect
        % vertices belonging to different sets.
        %
        % The element (i,j) of an biadjacency matrix is 1 for an
        % undirected bipartite graph iff vertex i of the first set of
        % vertices is linked with vertex j in the second set of vertices.
        % Alternatively, the value of the matrix element is the weight
        % associated to this edge. 
        %
        % Args:
        %   biadj(double (\*,\*)): array of vertices linked by an edge
        %
        % Returns:
        %   graph (`.UndirectedGraph`)
        %
        % Example:
        %   >>> replab.UndirectedGraph.fromBiadjacencyMatrix([1 0 1; 1 1 0]);

            % Construct the associated full adjacency matrix
            adj = [zeros(size(biadj,1)*[1 1]), biadj;
                   zeros(size(biadj,2), sum(size(biadj)))];
               
            self = replab.UndirectedGraph.fromAdjacencyMatrix(adj);
        end

        function self = fromDirectedGraph(graph)
            % Define an undirected graph from a directed one
            %
            % This function returns an undirected graph with the same
            % connectivity as the input directed graph.
            %
            % Args:
            %   graph (`.DirectedGraph`)
            %
            % Returns:
            %   graph (`.UndirectedGraph`)

            assert(isa(graph, 'replab.DirectedGraph'), 'Input is not a directed graph');

            self = replab.UndirectedGraph.fromEdges(graph.edges, graph.nVertices, graph.colors, graph.weights);
        end
        
    end

    methods
        
        function s = headerStr(self)
        % Header string representing the object
        %
        % Returns:
        %     charstring: description of the object

            s = sprintf('Undirected graph with %d vertices and %d edges', self.nVertices, size(self.edges,1));
        end
        
    end
    
    methods % Methods

        function graph = toDirectedGraph(self)
            % Adds directionality to the graph's edges
            %
            % This function returns an directed graph with the same
            % connectivity as the current graph.
            %
            % Args:
            %   graph (`.UndirectedGraph`)
            %
            % Returns:
            %   graph (`.DirectedGraph`)
            
            graph = replab.DirectedGraph.fromUndirectedGraph(self);
        end
        
        function adj = adjacencyMatrix(self)
        % Returns the adjacency matrix of a graph
        %
        % Args:
        %   graph (`.UndirectedGraph`)
        %
        % Returns:
        %   adj (double (\*,\*)): adjacency matrix

            adj = adjacencyMatrix@replab.DirectedGraph(self);
            adj = adj + adj.' - diag(diag(adj));
            adj = full(adj);
        end
        
        function L = laplacian(self)
        % Returns the graph Laplacian
        %
        % This is a generalized laplacian which also takes into account the
        % graph coloring. It reduces to the standard laplacian when
        % self.color = 0. A graph laplacian is semidefinite positive
        %
        % Args:
        %   graph (`.UndirectedGraph`)
        %
        % Returns:
        %   L (double (\*,\*)): Laplacian of the graph
            L = self.cached('laplacian', @() self.computeLaplacian);
        end
        
        function L = computeLaplacian(self)
        % Computes the graph Laplacian
        
            M = self.adjacencyMatrix();
            D = diag(sum(M));
            
            colors = abs(self.colors);
            if numel(colors) == 1
                colors = colors*ones(1,self.nVertices);
            end
            
            L = D - M + diag(colors);
        end
        
        function Kt = heatKernel(self, nbTimes)
        % The heat kernel
        %
        % Evaluates the heat kernel of the graph at several times.
        %
        % Args:
        %   graph (`.UndirectedGraph`)
        %   nbTimes (integer, optional): the number of distinct times at 
        %     which to evaluate the kernel (default is 10)
        %
        % Returns:
        %   Kt (double (\*,\*,\*)): The heat kernel
            
            defaultNbTimes = 3;
            
            if nargin < 2
                nbTimes = defaultNbTimes;
            else
                nbTimes = round(nbTimes);
                assert(nbTimes > 0, 'The number of time steps should be a positive integer');
            end
            
            if (nbTimes == defaultNbTimes)
                Kt = self.cached('heatKernel', @() self.computeHeatKernel(nbTimes));
            else
                Kt = self.computeHeatKernel(nbTimes);
            end
        end
        
        function Kt = computeHeatKernel(self, nbTimes)
        % Computes the heat kernel
        
            [phi, lambda] = eig(self.laplacian());
            lambda = diag(lambda);

            co = 0;
            Kt = zeros(nbTimes, self.nVertices, self.nVertices);
            for t = 1/nbTimes:1/nbTimes:1
                co = co + 1;
                Kt(co,:,:) = exp(-lambda(1)*t)*phi(:,1)*phi(:,1)';
                for i = 2:self.nVertices
                    Kt(co,:,:) = Kt(co,:,:) + permute(exp(-lambda(i)*t)*phi(:,i)*phi(:,i)', [3 1 2]);
                end
            end
        end
        
    end
    
end
