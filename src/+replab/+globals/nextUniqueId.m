function id = nextUniqueId
% Returns a unique integer
%
% Returns:
%   integer: Unique integer to be used as ID
    persistent lastId
    if isempty(lastId)
        lastId = 0;
    end
    id = lastId + 1;
    if lastId == id % wrap around when lastId = 2^53
        id = 0;
    end
    lastId = id;
end
