classdef Graph < replab.Obj
% Abstract class for immutable graphs
%
% Vertices of the graph are numbered continuously from 1 to nVertices

    properties (SetAccess = protected)
        nVertices % (integer): number of vertices in the graph
        edges % (integer (\*,2)): list of edges between vertices
        colors % (double (\*,1)): color of each all vertices; if a scalar, identical for all vertices
        weights % (double (\*,1)): weights associated to the edges; if a scalar, identical for all edges
    end

    methods (Access = protected)
        
        function self = Graph(nVertices, edges, colors, weights)
        % Construct a graph
        %
        % See also:
        %   `+replab.DirectedGraph`
        %   `+replab.UndirectedGraph`

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
            error('abstract');
        end

        function self = fromEdges(edges, nVertices, colors, weights)
        % Constructs a directed graph from a liste of edges
            error('abstract');
        end

        function self = fromBiadjacencyMatrix(biadj)
        % Constructs a directed graph from a biadjacency matrix
            error('abstract');
        end
        
        function self = fromFile(fileName)
        % Loads a graph from a Bliss file
            error('abstract');
        end
   end

    methods
        
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
            error('abstract');
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
            error('abstract');
        end
        
        function Kt = heatKernel(self, nbTimes)
        % The heat kernel
        %
        % Evaluates the heat kernel of the graph at several times.
        % For graphs with more than 100 vertices, an approximate heat
        % kernel is computed.
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
        
            if self.nVertices > 100
                % We compute an approximate eigendecomposition
                [phi, lambda] = eigs(self.laplacian());
            else
                [phi, lambda] = eig(self.laplacian());
            end
            lambda = diag(lambda);

            co = 0;
            Kt = zeros(nbTimes, self.nVertices, self.nVertices);
            for t = 1/nbTimes:1/nbTimes:1
                co = co + 1;
                Kt(co,:,:) = exp(-lambda(1)*t)*phi(:,1)*phi(:,1)';
                for i = 2:length(lambda)
                    Kt(co,:,:) = Kt(co,:,:) + permute(exp(-lambda(i)*t)*phi(:,i)*phi(:,i)', [3 1 2]);
                end
            end
        end
        
        function autoG = automorphismGroup(self)
        % Returns the automorphism group of the graph
        %
        % Args:
        %   graph (`.DirectedGraph`)
        %
        % Returns:
        %   autoG (`.PemutationGroup`)
        
            autoG = self.cached('automorphismGroup', @() replab.bsgs.GraphAutomorphism(self).subgroup);
        end
        
    end
    
end
