function [ok, coloring] = graphIsBipartite(pairs)
% Tests if a graph is bipartite
%
% Checks whether the graph described by a list of pairs of
% connected vertices is bipartite, i.e. whether the graph can be
% colored with 2 colors. Returns the coloring in case of success.
%
% Args:
%     pairs: nx2 array listing the (undirected) edges connecting
%         pairs of vertices. Should only contain strictly positive
%         integers.
%
% Returns:
%     ok: 1 if the graph is bipartite
%         0 if the graph requires more than 2 colors to be
%           colored
%     coloring: coloring of the vertices in case of success

    % We list all vertices
    vertices = unique(pairs(:))';

    % We assign a neutral color to every vertex
    colors = zeros(1,max(vertices));

    coloredVertices = zeros(size(vertices));
    coColoredVertices = 0;
    coV = 0;
    ok = 0;
    while coColoredVertices < length(vertices)
        % We assign color 1 to one new un-reached vertex
        newVertex = vertices(find(setdiff(vertices, coloredVertices), 1, 'first'));
        colors(newVertex) = 1;
        coColoredVertices = coColoredVertices + 1;
        coloredVertices(coColoredVertices) = newVertex;

        % We assign alternating colors to neighbours
        while coV < length(coloredVertices)
            coV = coV + 1;

            % We pick the next 
            currentVertex = coloredVertices(coV); % current vertex
            currentColor = colors(currentVertex); % color of current vertex

            % We color all neighbours of this vertex
            for i = 1:size(pairs,1)
                newVertex = 0; % new vertex
                newColor = mod(currentColor,2)+1; % color of neighbouring vertices

                if pairs(i,1) == currentVertex
                    newVertex = pairs(i,2);
                elseif pairs(i,2) == currentVertex
                    newVertex = pairs(i,1);
                end

                if (newVertex ~= 0) && (newVertex ~= currentVertex)
                    if ismember(newVertex, coloredVertices)
                        % We already passed through this vertex, checking
                        % that coloring is consistent with previous choice
                        if colors(newVertex) ~= newColor
                            coloring = [];
                            return;
                        end
                    else
                        % We found a new vertex, let us color it
                        colors(newVertex) = newColor;
                        coColoredVertices = coColoredVertices + 1;
                        coloredVertices(coColoredVertices) = newVertex;
                    end
                end
            end
        end
    end
    ok = 1;
    coloring = [coloredVertices; colors(coloredVertices)];
end
