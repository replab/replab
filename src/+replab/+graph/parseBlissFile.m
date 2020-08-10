function [nVertices, edges, colors] = parseBlissFile(fileName)
% Parses a Bliss file
%
% Following the specifications from
% http://www.tcs.hut.fi/Software/bliss/fileformat.shtml
%
% Args:
%   fileName (charstring) : Name of the file to be loaded
%
% Returns:
% --------
%   nVertices (integer)
%   edges (integer (\*, 2))
%   colors (integer (\*, 1))

[fid, message] = fopen(fileName);

assert(fid >= 0, message);

problemLineFound = false;
nVertices = 0;
edges = zeros(0,2);
while true
    % Get the next line
    tline = fgetl(fid);

    if ~ischar(tline)
        break;
    end

    % Parse the line
    switch tline(1)
        case 'c'
            % This is a comment line
        case 'p'
            % The 'problem line
            assert(~problemLineFound, 'Multiple graph definitions in a unique file');
            problemLineFound = true;
            assert(isequal(tline(1:7), 'p edge '), 'Format mismatch');
            ve = regexp(tline(8:end), '\d+', 'match');
            assert(length(ve) == 2, 'Number of vertices and of edges expected');
            nVertices = str2double(ve{1});
            %edges = zeros(str2double(ve{2}),2);
            colors = zeros(1,nVertices);
        case 'n'
            % Color assignment
            assert(problemLineFound, 'Graph definition missing');
            vc = regexp(tline(3:end), '\d+', 'match');
            assert(length(vc) == 2, 'Vertex number and color expected');
            assert(str2double(vc{1}) <= nVertices);
            colors(str2double(vc{1})) = str2double(vc{2});
        case 'e'
            % Edge specification
            assert(problemLineFound, 'Graph definition missing');
            vs = regexp(tline(3:end), '\d+', 'match');
            assert(length(vs) == 2, 'Vertex number and color expected');
            assert(str2double(vs{1}) <= nVertices);
            assert(str2double(vs{2}) <= nVertices);
            edges = [edges; str2double(vs{1}) str2double(vs{2})];
    end
end

fclose(fid);

% Is there only one color?
if length(unique(colors)) == 1
    colors = colors(1);
end
