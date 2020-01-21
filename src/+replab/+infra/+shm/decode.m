function path = decode(id)
    if isequal(id, 'root__')
        path = {};
    else
        id = strrep(id, '__', '$');
        path = strsplit(id, '$');
    end
end
