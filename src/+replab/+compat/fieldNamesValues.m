function [names values] = fieldNamesValues(obj)
% Get the public field names and values of a struct or an object, taking care of an Octave bug
%
% Args:
%   obj: Object or struct to list the field names of
%
% Returns
% -------
%   names:
%     cell{:,1} of charstring: Field names
%   values:
%     cell{:,1} of values: Field values
    if replab.compat.isOctave
        % Octave has a bug in some versions: it prints a warning that the object is converted to a struct
        % and lists private/protected properties as well
        err = [];
        try
            candidates = fieldnames(obj);
        catch
            err = lasterror;
        end
        if ~isempty(err)
            rethrow(err);
        end
        warning(prev); % reset warnings
    else
        candidates = fieldnames(obj);
    end
    names = {};
    values = {};
    for i = 1:length(candidates)
        name = candidates{i};
        try
            value = obj.(name);
            names{end+1, 1} = name;
            values{end+1, 1} = value;
        catch
        end
    end
end
