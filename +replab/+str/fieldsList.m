function [names values] = fieldsList(obj)
% Returns a list of field names and values for the given object
    if isa(obj, 'replab.Str')
        try
            % We use the replab.Str.fieldsList method in
            % case the object implements that interface 
            [names values] = obj.fieldsList;
            return
        catch
        end
    end
    if isstruct(obj) || isobject(obj)
        % fallback using Matlab field enumeration
        names = fieldnames(obj);
        values = cellfun(@(f) obj.(f), names, 'uniform', 0);
    else
        names = {};
        values = {};
    end
end
