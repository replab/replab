function [s overLimit] = shortStr(obj, maxColumns)
% Default implementation for replab.shortStr
%   
% Thank you Matlab for being a regular language, we totally don't have to consider a thousand
% particular cases below.
    overLimit = false;
    if isscalar(obj) % prints scalars using the relevant method
        if isa(obj, 'vpi')
            s = strtrim(num2str(obj));
        elseif isnumeric(obj)
            s = num2str(obj);
        elseif isa(obj, 'sym')
            s = char(obj);
        else
            s = class(obj);
        end
    elseif ischar(obj) && length(obj) > 0 && size(obj, 1) == 1
        s = ['''' replab.str.escape(obj) ''''];
    elseif isstruct(obj) || isobject(obj)
        s = replab.str.tinyStr(obj);
        [names values] = replab.str.fieldsList(obj);
        s = [s ' with ' replab.str.pluralize(length(f), 'field') ': ' strjoin(names, ', ')];
    elseif numel(obj) == 0 % prints an empty object using the tiny description which has nice defaults
        s = replab.str.tinyStr(obj);
    elseif numel(obj) > maxColumns % prints too big to fit objects using the tiny description
        s = replab.str.tinyStr(obj);
    elseif isvector(obj)
        [lp rp] = replab.str.brackets(obj);
        elements = arrayfun(@(i) replab.str.cellStr(obj, maxColumns, i), 1:length(obj), 'uniform', 0);
        s = [lp strjoin(elements, ', ') rp];
    elseif ismatrix(obj)
        [lp rp] = replab.str.brackets(obj);
        rowFun = @(r) strjoin(arrayfun(@(c) replab.str.cellStr(obj, maxColumns, r, c), 1:size(obj, 2), 'uniform', 0), ', ');
        s = [lp strjoin(arrayfun(rowFun, 1:size(obj, 1), 'uniform', 0), '; ') rp];
    else
        s = tinyStr(obj);
    end
    if length(s) > maxColumns
        s = replab.str.tinyStr(obj);
        overLimit = true;
    else
        overLimit = false;
    end
end
