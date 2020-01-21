function [names values] = fieldsList(obj)
% Returns a list of field names and values for the given object, returned as column vectors
    if isstruct(obj) || isobject(obj)
        [names values] = replab.compat.fieldNamesValues(obj);
    else
        names = {};
        values = {};
    end
    if isa(obj, 'replab.Str')
        hidden = obj.hiddenFields;
        [names I] = setdiff(names, hidden);
        values = values(I);
        [names1 values1] = obj.additionalFields;
        names = replab.str.horzcatForce(names, names1);
        values = replab.str.horzcatForce(values, values1);
    end
end
