function rmdirRec(path)
% Removes the folder with the given path including its contents without asking confirmation
%
% Args:
%   path (charstring): Path of the folder to remove
%
% Raises:
%   Throws errors if the folder doesn't exist, or if it can't be removed
    assert(exist(path) == 7, 'The path %s is not a folder', path);
    if replab.settings.isOctave
        confirm_recursive_rmdir (false, 'local');
    end
    rmdir(path, 's');
end
