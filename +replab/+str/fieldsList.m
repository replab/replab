function [names values] = fieldsList(obj)
% Returns a list of field names and values for the given object,
% returned as column vectors
    if isstruct(obj) || isobject(obj)
        % Octave has a bug in some versions: it prints a
        % warning that the object is converted to a struct
        % and lists private/protected properties as well
        prev = warning('off'); % turn off warnings
        candidates = fieldnames(obj);
        warning(prev);
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
    else
        names = {};
        values = {};
    end
    names = names(:);
    values = values(:);
    if isa(obj, 'replab.Str')
        hidden = obj.hiddenFields;
        [names I] = setdiff(names, hidden);
        values = values(I);
        [names1 values1] = obj.additionalFields;
        names1 = names1(:);
        values1 = values1(:);
        names = vertcat(names, names1);
        values = vertcat(values, values1);
    end
end
