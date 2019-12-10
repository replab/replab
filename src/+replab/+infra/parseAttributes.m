function attributes = parseAttributes(string)
% Parses property/method block attributes into a structure
%
% An example of input is 'SetAccess = private, Access = protected'
    attrs = strsplit(string, ',');
    attributes = struct;
    for i = 1:length(attrs)
        tokens = regexp(attrs{i}, '\s*([A-Za-z]+)\s*=\s*([A-Za-z]+)\s*', 'tokens', 'once');
        if length(tokens) == 0
            error([attrs{i} ' is not a recognized property attribute']);
        end
        attributes.(tokens{1}) = tokens{2};
    end
end
