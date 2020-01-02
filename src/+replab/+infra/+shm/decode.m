function path = decode(id)
    id = strrep(id, '__', '$');
    path = strsplit(id, '$');
end
