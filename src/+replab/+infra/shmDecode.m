function path = shmDecode(id)
    id = strrep(id, '__', '$');
    path = strsplit(id, '$');
end
