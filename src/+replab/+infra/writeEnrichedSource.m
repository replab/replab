function writeEnrichedSource(codeBase, docSrcPath, packageElement)
% Writes the enriched source code for Sphinx documentation generation
%
% Args:
%   docSrcPath (charstring): Base folder for doctest generation, not including trailing path separator
%                            That folder must exist.
    if isa(packageElement, 'replab.infra.Class1')
        src = replab.infra.classWithTOC(codeBase, packageElement);
    else
        src = fileread(packageElement.fullFilename);
    end
    base = docSrcPath;
    for i = 1:length(packageElement.packageNameParts)
        part = ['+' packageElement.packageNameParts{i}];
        switch exist(fullfile(base, part))
          case 0
            assert(mkdir(base, part), sprintf('Could not create %s in directory %s', part, base));
          case 7
            % directory already exists
          otherwise
            error(sprintf('%s already exists in directory %s, but is not a folder', part, base));
        end
        base = fullfile(base, part);
    end
    filename = fullfile(base, [packageElement.name '.m']);
    fid = fopen(filename, 'w');
    fwrite(fid, src);
    fclose(fid);
end
