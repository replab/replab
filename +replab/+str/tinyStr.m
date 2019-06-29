function s = tinyStr(obj)
% Returns a tiny string description -- class name and size (when relevant)
    if isa(obj, 'vpi')
        s = sprintf('vpi(~%.2e)', double(obj));
    elseif isequal(obj, {})
        s = '{}';
    elseif isequal(obj, [])
        s = '[]';
    elseif numel(obj) == 0
        s = [replab.str.sizeStr(size(obj)) ' empty ' class(obj) ' array'];
    elseif numel(obj) == 1
        s = class(obj);
    else
        s = [replab.str.sizeStr(size(obj)) ' ' class(obj)];
    end
end
