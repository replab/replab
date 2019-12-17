function id = shmEncode(parts)
    if isempty(parts)
        id = 'root__';
    else
        id = strjoin(parts, '__');
    end
end
