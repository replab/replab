classdef DirectedGraph < replab.Obj
% Describes an immutable directed graph
%
% Vertices of the graph are numbered continuously from 1 to nVertices

    properties (SetAccess = protected)
        nVertices % (integer): number of vertices in the graph
        edges % (integer (\*,2)): list of edges between vertices
        colors % (double (\*,1)): color of each all vertices; if a scalar, identical for all vertices
        weights % (double (\*,1)): weights associated to the edges; if a scalar, identical for all edges
    end

    methods (Access = protected)
        
        function self = DirectedGraph(nVertices, edges, colors, weights)
        % Construct a directed graph
        %
        % Do not use this function directly, rather use another
        % constructor such as ``.fromBlocks``.
        %
        % See also:
        %   `.fromBlocks`
        %   `.check`

            if isempty(edges)
                edges = zeros(0,2);
            end
            
            if isempty(colors)
                colors = 0;
            end
            
            if isempty(weights)
                weights = 1;
            end
            
            assert(numel(nVertices) == 1, 'Number of elements should be a scalar');
            assert(isempty(edges) || (nVertices >= max(max(edges))), 'Insufficient number of vertices');
            assert(isequal(size(edges,2), 2), 'Incorrect size for edge argument');
            assert(length(size(edges)) <= 2, 'Edge argument should be a 2-dimensional array');

            assert((numel(colors) == 1) || (numel(colors) == nVertices), 'Number of vertices and of colors incompatible');
            assert(isvector(colors), 'Colors argument should be a vector');
            assert((numel(weights) == 1) || (numel(weights) == size(edges,1)), 'Number of edges and of weights incompatible');
            assert(isvector(weights), 'Weights argument should be a vector');

            self.nVertices = nVertices;
            self.edges = edges;
            self.colors = colors(:);
            self.weights = weights(:);
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
        
    end

    methods
        
        function s = headerStr(self)
        % Header string representing the object
        %
        % Returns:
        %     charstring: description of the object

            s = sprintf('Directed graph with %d vertices and %d edges', self.nVertices, size(self.edges,1));
        end
        
        function names = hiddenFields(self)
        % Overload of `+replab.Str.hiddenFields`

            names = hiddenFields@replab.Str(self);
            names{1, end+1} = 'nVertices';
            names{1, end+1} = 'colors';
            names{1, end+1} = 'weights';
        end

        function [names, values] = additionalFields(self)
        % Overload of `+replab.Str.additionalFields`

            [names, values] = additionalFields@replab.Str(self);

            % We only display nontrivial colorings
            if ~isequal(self.colors, 0)
                names{1, end+1} = 'colors';
                values{1, end+1} = self.colors;
            end
            
            % We only display nontrivial weights
            if ~isequal(self.weights, 1)
                names{1, end+1} = 'weights';
                values{1, end+1} = self.weights;
            end
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
        
        function adj = adjacencyMatrix(self)
        % Returns the adjacency matrix of a graph
        %
        % Args:
        %   graph (`.DirectedGraph`)
        %
        % Returns:
        %   adj (double (\*,\*)): adjacency matrix
            
            adj = self.cached('adjacencyMatrix', @() self.computeAdjacencyMatrix);
        end
        
        function adj = computeAdjacencyMatrix(self)
        % Computes the adjacency matrix

            adj = replab.graph.edge2adj(self.edges, self.nVertices, self.weights);
            adj = full(adj);
        end
        
        function ok = isBipartite(self)
        % Tests if a graph is bipartite
        %
        % Args:
        %   graph (`.DirectedGraph`)
        %
        % Returns:
        %   ok (bool): true iff the graph is bipartite
            
            ok = replab.graph.graphIsBipartite(self.edges);
        end

        function P = connectedComponents(self)
        % Returns the sets of vertices corresponding to connected components
        %
        % For edges = [1 3] and n = 3, it returns the partition {[1 3] [2]}
        %
        % Note: Graph connectedness does not take into account the graph
        % directionality.
        %
        % Args:
        %   graph (`.DirectedGraph`)
        %
        % Returns:
        %   P (`.Partition`): Partitioning of the vertices

            if isempty(self.edges)
                % Trivial case
                blockIndex = 1:self.nVertices;
                blocks = num2cell(1:self.nVertices, 1);
            else
                [blocks, blockIndex] = replab.graph.connectedComponents(self.edges);

                % If some elements are isolated, we add them
                connectedVertices = [blocks{:}];
                isolatedVertices = setdiff(1:self.nVertices, connectedVertices);
                nbConnectedSets = length(blocks);

                if length(isolatedVertices) >= 1
                    % allocate memory
                    blocks{nbConnectedSets + length(isolatedVertices)} = [];
                    if length(blockIndex) < self.nVertices
                        blockIndex(self.nVertices) = 0;
                    end

                    % assign values
                    co = nbConnectedSets;
                    for i = 1:length(isolatedVertices)
                        co = co + 1;
                        blocks{co} = isolatedVertices(i);
                        blockIndex(isolatedVertices(i)) = co;
                    end
                end
            end

            % Construct the Partition object
            P = replab.Partition(blockIndex, blocks);
        end
        
        function autoG = automorphismGroup(self)
        % Returns the automorphism group of the graph
        %
        % Args:
        %   graph (`.DirectedGraph`)
        %
        % Returns:
        %   autoG (`.PemutationGroup`)
        
            autoG = self.cached('automorphismGroup', @() replab.graph.AutomorphismSearch(self).subgroup);
        end
        
    end
    
end
