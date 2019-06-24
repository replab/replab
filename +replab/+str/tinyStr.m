function s = tinyStr(obj)
% Returns a tiny string description -- class name and size (when relevant)
    if iscell() 
    s = [replab.str.sizeStr(size(obj)) ' ' class(obj)];
end
