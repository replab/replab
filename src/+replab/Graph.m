classdef Graph < replab.Obj
% Describes a graph
%
% Vertices of the graph are numbered continuously from 1 to nVertices

    properties (SetAccess = protected)
        nVertices % (integer): number of vertices in the graph
        edges % (integer (\*,2)): list of edges between vertices
        weights % (double (\*,1)): weights associated to the edges
    end

    methods (Access = protected)
        
        function self = Graph(nVertices, edges, weights)
        % Construct a Graph
        %
        % Do not use this function direclty, rather use another
        % constructor such as ``.fromBlocks``.
        %
        % See also:
        %   `.fromBlocks`
        %   `.check`

            self.nVertices = nVertices;
            self.edges = edges;
            self.weights = weights;
        end

    end
    
    methods (Static) % Constructors

        function self = fromAdjacencyMatrix(adj)
        % Constructs a Graph from an adjacency matrix
        %
        % The element (i,j) of an adjacency matrix is 1 iff vertex i is
        % linked to vertex j. Alternatively, the value of the matrix
        % element is the weight associated to this edge.
        %
        % Args:
        %   adj (double (\*,\*)):
        %
        % Returns:
        %   graph (`.Graph`)
        %
        % Example:
        %   >>> replab.Graph.fromAdjacencyMatrix([0 1 1; 1 0 1; 1 1 0]);
            
            [edges, nVertices, weights] = replab.graph.adj2edge(adj);
            self = replab.Graph(nVertices, edges, weights);
        end

        function self = fromEdges(edges, nVertices, weights)
        % Constructs a Graph from a liste of edges
        %
        % The element (i,j) of an adjacency matrix is 1 iff vertex i is
        % linked to vertex j. Alternatively, the value of the matrix
        % element is the weight associated to this vertex.
        %
        % Args:
        %     edges (integer (\*,2)): array of vertices linked by an edge
        %     nVertices (integer, optional): number of vertices
        %     weights (double (\*,1), optional): weight associated to each edge
        %
        % Returns:
        %   graph (`.Graph`)
        %
        % Example:
        %   >>> replab.Graph.fromEdges([1 2; 2 3; 3 1], 3);

            if nargin < 2
                if isempty(edges)
                    nVertices = 0;
                else
                    nVertices = max(max(edges));
                end
            end
            
            if nargin < 3
                if isempty(edges)
                    weights = [];
                else
                    weights = 1;
                end
            end
            
            if isempty(edges)
                edges = zeros(0,2);
            else
                assert(nVertices >= max(max(edges)), 'Number of vertices insufficient');
                assert(isequal(size(edges,2), 2), 'Incorrect size for edge argument');
                assert(length(size(edges)) <= 2, 'Edge argument should be a 2-dimensional array');
                assert(isvector(weights), 'Weight argument should be a vector');
                weights = weights(:); % Force verticality
                assert((size(weights,1) == size(edges,1)) || (size(weights,1) == 1), 'Incompatible number of weight and edges');
            end

            self = replab.Graph(nVertices, edges, weights);
        end

        function self = fromBiadjacencyMatrix(biadj)
        % Constructs a Graph from a biadjacency matrix
        %
        % A bipartite (undirected) graph is one in which vertices can
        % be divided into two sets, such that edges only connect
        % vertices belonging to different sets.
        %
        % The element (i,j) of an biadjacency matrix is 1 for an
        % undirected bipartite graph iff vertex i of the first set of
        % verices is linked to vertex j in the second set of vertices.
        % Alternatively, the value of the matrix element is the weight
        % associated to this edge. 
        %
        % Args:
        %   biadj(double (\*,\*)): array of vertices linked by an edge
        %
        % Returns:
        %   graph (`.Graph`)
        %
        % Example:
        %   >>> replab.Graph.fromBiadjacencyMatrix([1 0 1; 1 1 0]);

            % Construct the associated full adjacency matrix
            adj = [zeros(size(biadj,1)*[1 1]), biadj;
                   biadj.' zeros(size(biadj,2)*[1 1])];
            [edges, nVertices, weights] = replab.graph.adj2edge(adj);
            self = replab.Graph(nVertices, edges, weights);
        end

    end

    methods % Methods

        function adj = adjacencyMatrix(self)
        % Returns the adjacency matrix of a graph
        %
        % Args:
        %   graph (`.Graph`)
        %
        % Returns:
        %   adj (double (\*,\*)): adjacency matrix

            adj = replab.graph.edge2adj(self.edges, self.nVertices, self.weights);
            adj = full(adj);
        end
        
        function ok = isBipartite(self)
        % Tests if a graph is bipartite
        %
        % Args:
        %   graph (`.Graph`)
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
        % Args:
        %   graph (`.Graph`)
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
        
    end
    
end
