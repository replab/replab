function attributes = parseAttributes(string)
% Parses property/method block attributes into a structure
%
% An example of input is 'SetAccess = private, Static, Access = protected'
    attrs = strsplit(string, ',');
    attributes = struct;
    for i = 1:length(attrs)
        attr = strtrim(attrs{i});
        if contains(attr, '=')
            % this attribute is a key value pair, like 'Access = public'
            tokens = regexp(attr, '^([A-Za-z]+)\s*=\s*([A-Za-z]+)$', 'tokens', 'once');
            if isempty(tokens)
                error([attr ' is not a recognized property attribute']);
            end
            attributes.(tokens{1}) = tokens{2};
        else
            % this attribute is just a flag, like 'Static'
            tokens = regexp(attr, '^([A-Za-z]+)$', 'tokens', 'once');
            if isempty(tokens)
                error([attr ' is not a recognized property attribute']);
            end
            attributes.(tokens{1}) = true;
        end
    end
end
