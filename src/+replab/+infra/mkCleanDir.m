function mkCleanDir(baseDir, folderName, logFunction)
% Creates an empty directory, removing any previous directory
%
% Args:
%   baseDir (charstring): Base directory, which must exist
%   folderName (charstring): Folder to create, or to empty, in the base directory
%   logFunction (function_handle, optional): Log function to call, single string argument
%
% Raises:
%   An error if the base directory does not exist, or if any operation fails
    if nargin < 3 || isempty(logFunction)
        logFunction = @(s) [];
    end
    switch exist(baseDir) % check that the base directory exists
      case 0
        error('Base directory ''%s'' does not exist', baseDir);
      case 7
        % it's ok
      otherwise
        error('Path ''%s'' exists, but is of a different type', baseDir);
    end
    path = fullfile(baseDir, folderName);
    switch exist(path)
      case 0
        logFunction(sprintf('Folder ''%s'' does not exist in ''%s'', will be created', folderName, baseDir));
      case 7
        logFunction(sprintf('Folder ''%s'' exists in ''%s'', will be deleted', folderName, baseDir));
        replab.compat.rmdirRec(path);
      otherwise
        error('Element ''%s'' exists in ''%s'', but is not a folder', folderName, baseDir);
    end
    mkdir(path);
end
