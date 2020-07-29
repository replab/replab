classdef Graph < replab.Obj
% Describes a graph

    properties (SetAccess = protected)
        nVertices % (integer): number of vertices in the graph
        edges % (integer (\*,2)): list of edges between vertices
        weights % (double (\*,1)): weights associated to the edges
    end

    methods (Static) % Constructors

        function self = fromCoincidenceMatrix(adj)
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
            %   graph (``replab.Graph``)
            %
            % Example:
            %   >>> replab.Graph.fromCoincidenceMatrix([0 1 1; 1 0 1; 1 1 0]);
            
            self = replab.Graph;
            [self.edges, self.nVertices, self.weights] = replab.graph.adj2edge(adj);
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
            %     nVertices (integer): number of vertices
            %     weights (double (\*,1), optional): weight associated to each edge
            %
            % Returns:
            %   graph (``replab.Graph``)
            %
            % Example:
            %   >>> replab.Graph.fromEdges([1 2; 2 3; 3 1], 3);

            if nargin < 3
                weights = 1;
            end
            
            if isempty(edges)
                edges = zeros(0,2);
            end

            assert(nVertices >= max(max(edges)));
            assert(isequal(size(edges,2), 2));
            assert(length(size(edges)) <= 2);
            assert(isequal(size(weights,2), 1));
            assert((size(weights,1) == size(edges,1)) || (size(weights,1) == 1));

            self = replab.Graph;
            self.nVertices = nVertices;
            self.edges = edges;
            self.weights = weights;
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
            %   graph (``replab.Graph``)
            %
            % Example:
            %   >>> replab.Graph.fromBiadjacencyMatrix([1 0 1; 1 1 0]);

            % Construct the associated full adjacency matrix
            adj = [zeros(size(biadj,1)*[1 1]), biadj;
                   biadj.' zeros(size(biadj,2)*[1 1])];
            self = replab.Graph;
            [self.edges, self.nVertices, self.weights] = replab.graph.adj2edge(adj);
        end

    end

    methods % Methods

        function adj = adjacencyMatrix(self)
            % adjacency matrix of a graph
            %
            % Args:
            %   graph (``replab.Graph``)
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
            %   graph (``replab.Graph``)
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
        %   graph (``replab.Graph``)
        %
        % Returns:
        %   P (``replab.Partition``): Partitioning of the vertices

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
                    blockIndex(self.nVertices) = 0;

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
