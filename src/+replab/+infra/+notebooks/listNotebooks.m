function list = listNotebooks
% finds all jupyter notebooks in the sphinx folder
%
% The result is split in three parts: the sphinx path, any subfolder,
% finally the notebook name.
%
% Results:
%   cell array of charstring: list of notebook

    % Make sure that we are called from the right folder
    assert(isequal(pwd, replab.globals.replabPath));
    
    % find all matlab files
    scan = dir('sphinx/**/*.m');
    
    % keep only notebooks
    toKeep = cellfun(@(x) isempty(strfind(x, '/_src')), {scan.folder});
    scan = scan(toKeep);

    list = cell(length(scan), 3);
    for i = 1:length(scan)
        list{i,1} = [pwd, '/sphinx'];
        tmp = scan(i).folder;
        list{i,2} = tmp(length(list{i,1})+2:end);
        list{i,3} = scan(i).name;
    end 
end