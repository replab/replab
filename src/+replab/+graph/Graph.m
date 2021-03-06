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
        % Returns:
        %   adj (double (\*,\*)): adjacency matrix
        %
        % Example:
        %   >>> replab.UndirectedGraph.fromEdges([1 3]).adjacencyMatrix
        %     0     0     1
        %     0     0     0
        %     1     0     0
     
            adj = self.cached('adjacencyMatrix', @() self.computeAdjacencyMatrix);
        end
        
        function adj = computeAdjacencyMatrix(self)
        % Computes the adjacency matrix
            error('abstract');
        end
        
        function deg = degrees(self)
        % Returns the degrees of all vertices
            error('abstract');
        end
        
        function deg = degree(self, v)
        % Returns the degrees of vertex v
            error('abstract');
        end
        
        function deg2 = secondOrderDegrees(self)
        % Returns the number of vertices at a distance 2 of all vertices
        %
        % Returns:
        %   deg (double (1,\*)): list of degrees
        %
        % Example:
        %   >>> replab.UndirectedGraph.fromEdges([1 3]).secondOrderDegrees
        %     1     0     1

            deg2 = arrayfun(@(x) self.secondOrderDegree(x), [1:self.nVertices]);
        end
        
        function deg2 = secondOrderDegree(self, v)
        % Returns the number of vertices at a distance 2 of vertex v
            error('abstract');
        end
        
        function degN = degreesSequences(self)
        % Returns the degrees sequence for all vertices
        %
        % In each sequence, a vertex can only be counted once.
        %
        % Returns:
        %   degN (integer (\*,\*)): list of degrees sequences
        %
        % Example:
        %   >>> replab.UndirectedGraph.fromEdges([1 3]).degreesSequences
        %     1
        %     0
        %     1

            degNs = arrayfun(@(x) self.degreesSequence(x), [1:self.nVertices], 'UniformOutput', false);
            
            maxLength = max(cellfun(@(x) length(x), degNs));
            
            degN = zeros(self.nVertices, maxLength);
            for i = 1:self.nVertices
                degN(i,1:length(degNs{i})) = degNs{i};
            end
        end
        
        function degN = degreesSequence(self, v)
        % Returns the degrees sequence for all vertices
        %
        % Each vertex can only be counted once. For directed graphs, this
        % corresponds to the outgoing degree sequence.
        %
        % Returns:
        %   degN (integer (1,\*)): sequence of degrees
        %
        % Example:
        %   >>> graph = replab.DirectedGraph.fromEdges([1 3]);
        %   >>> graph.degreesSequence(1)
        %     1
        
            adj = self.adjacencyMatrix;

            co = 0;
            selBefore = false(1,self.nVertices);
            sel = ([1:self.nVertices] == v);
            degN = [];
            while sum(sel) > 0
                co = co + 1;
                selBefore = selBefore | sel;
                sel = (sum(adj(sel,:) ~= 0,1) ~= 0);
                sel = sel & (~selBefore);
                if any(sel)
                    degN(co) = sum(sel);
                end
            end
        end
        
        function ok = isBipartite(self)
        % Tests if a graph is bipartite
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
        % Returns:
        %   P (`+replab.Partition`): Partitioning of the vertices
        %
        % Example:
        %   >>> replab.UndirectedGraph.fromEdges([1 3]).connectedComponents
        %     Partition '13|2'
        %     blockIndex: [1, 2, 1]
        %         blocks: {[1, 3], 2}
        %              n: 3

            if isempty(self.edges)
                % Trivial case
                blocks = num2cell(1:self.nVertices);
            else
                [~, blocks] = replab.graph.connectedComponents(self.nVertices, self.edges);

                % If some elements are isolated, we add them
                connectedVertices = [blocks{:}];
                isolatedVertices = setdiff(1:self.nVertices, connectedVertices);
                if length(isolatedVertices) >= 1
                    blocks = [blocks, num2cell(isolatedVertices)];
                end
            end

            % Construct the Partition object
            P = replab.Partition.fromBlocks(blocks);
        end
        
        function L = laplacian(self)
        % Returns the graph Laplacian
        %
        % This is a generalized laplacian which also takes into account the
        % graph coloring. It reduces to the standard laplacian when
        % self.color = 0. A graph laplacian is semidefinite positive
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
                for i = 2:length(lambda)
                    Kt(co,:,:) = Kt(co,:,:) + permute(exp(-lambda(i)*t)*phi(:,i)*phi(:,i)', [3 1 2]);
                end
            end
            %Kt = abs(Kt);
        end
        
        function autoG = automorphismGroup(self)
        % Returns the automorphism group of the graph
        %
        % Returns:
        %   autoG (`+replab.PermutationGroup`)
        
            autoG = self.cached('automorphismGroup', @() replab.bsgs.GraphAutomorphism(self).subgroup);
        end
        
    end
    
end
