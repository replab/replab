function el = look_up(s, path, ifEmpty)
% Looks up an element in a struct-encoded String Hash Map (SHM)
%
% We do not name this function ``look_up`` because Octave emits warnings if
% package functions shadow internal functions.
%
% Args:
%   s (struct): Hash-map encoded as a Matlab struct
%   path (row cell vector of charstring): Path to lookup
%   ifEmpty (optional): Value to return if empty
%
% Raises:
%   An error if the path is not found and `ifempty` is not given
    id = replab.infra.shm.encode(path);
    if isfield(s, id)
        el = s.(id);
    else
        if nargin < 2
            error(sprintf('Path %s not found', strjoin(path, '.')));
        end
        el = ifEmpty;
    end
end
