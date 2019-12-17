function s = shmUpdate(s, path, updateFun, ifEmpty)
% Updates the hash-map
    id = replab.infra.shmEncode(path);
    if isfield(s, id)
        s.(id) = updateFun(s.(id));
    else
        s.(id) = ifEmpty;
    end
end
