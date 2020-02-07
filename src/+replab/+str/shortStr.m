function s = shortStr(obj, maxColumns)
% Default implementation for `+replab.shortStr`
    overLimit = false;
    if isscalar(obj) % prints scalars using the relevant method
        if isa(obj, 'cell')
            s = ['{' replab.str.shortStr(obj{1}, maxColumns)];
            if length(s) < maxColumns
                s = [s '}'];
            end
        elseif isa(obj, 'vpi')
            s = num2str(obj);
            if size(s, 1) > 1
                % corrects for multiline vpi num2str
                s = s'; s = s(:)';
            end
            s = strtrim(s);
        elseif isa(obj, 'function_handle')
            s = func2str(obj);
        elseif islogical(obj)
            if obj
                s = 'true';
            else
                s = 'false';
            end
        elseif isnumeric(obj)
            s = num2str(obj);
        elseif isa(obj, 'sym')
            s = char(obj);
        elseif isa(obj, 'char')
            s = ['''' obj ''''];
        elseif isobject(obj)
            s = replab.headerStr(obj);
        elseif isstruct(obj)
            [names values] = replab.str.fieldsList(obj);
            s = ['struct with fields: ' strjoin(names, ', ')];
        else
            warning('Pretty printing not implemented')
            s = class(obj); % default print
        end
    elseif ischar(obj) && length(obj) > 0 && size(obj, 1) == 1
        s = ['''' replab.str.escape(obj) ''''];
    elseif numel(obj) == 0 % prints an empty object using the tiny description which has nice defaults
        s = replab.str.headerStr(obj);
    elseif numel(obj) > maxColumns % prints too big to fit objects using the tiny description
        s = replab.str.headerStr(obj);
    elseif isvector(obj)
        [lp rp] = replab.str.brackets(obj);
        elements = arrayfun(@(i) replab.str.cellStr(obj, maxColumns, i), 1:length(obj), 'uniform', 0);
        s = [lp strjoin(elements, ', ') rp];
        if size(obj, 1) > 1
            s = [s '.'''];
        end
    elseif ismatrix(obj)
        [lp rp] = replab.str.brackets(obj);
        rowFun = @(r) strjoin(arrayfun(@(c) replab.str.cellStr(obj, maxColumns, r, c), 1:size(obj, 2), 'uniform', 0), ', ');
        s = [lp strjoin(arrayfun(rowFun, 1:size(obj, 1), 'uniform', 0), '; ') rp];
    else % fallback
        s = replab.headerStr(obj);
    end
end
