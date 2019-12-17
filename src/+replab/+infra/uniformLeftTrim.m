function lines = uniformLeftTrim(lines)
% Trims uniformly the left whitespace on the given lines
%
% We count the amount of leading whitespace for each line, and computes the minimum
% among all lines, then remove that number of whitespace characters from each line.
%
% We thus trim unnecessary leading whitespace while preserving the semantic information
% the remaining whitespace contains.
    if ~isempty(lines)
        l = [];
        for i = 1:length(lines)
            cl = lines{i};
            if ~isempty(strtrim(cl))
                token = regexp(cl, '^(\s*)', 'tokens', 'once');
                if isempty(l)
                    l = length(token);
                else
                    l = min(l, length(token));
                end
            end
        end
        if ~isempty(l)
            lines = cellfun(@(x) x(l+1:end), lines, 'uniform', 0);
        end
    end
end
