function cleanDir(path, preserve, logFunction)
% Removes all files/subfolders from the given folder, except for the given list of preserved files
%
% Args:
%   path (charstring): Path to the folder to clean
%   preserve (cell(1,\*) of charstring): List of file/folder names to preserve
%   logFunction (function_handle, optional): Log function to call, single string argument
%
% Raises:
%   An error if the folder does not exist, or if any operation fails
    if nargin < 3 || isempty(logFunction)
        logFunction = @(s) [];
    end
    switch exist(path) % check that the base directory exists
      case 0
        error('Folder ''%s'' does not exist', folderPath);
      case 7
        % it's ok
      otherwise
        error('Path ''%s'' exists, but is of a different type than folder', folderPath);
    end
    children = dir(path);
    for i = 1:length(children)
        name = children(i).name;
        childpath = fullfile(path, name);
        if isequal(name, '.') || isequal(name, '..') || ismember(name, preserve)
            % do nothing
        elseif children(i).isdir
            logFunction(sprintf('Deleting subfolder %s\n', name));
            replab.compat.rmdirRec(childpath);
        else
            logFunction(sprintf('Deleting file %s\n', name));
            delete(childpath);
        end
    end
end
