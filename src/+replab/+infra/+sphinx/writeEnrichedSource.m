function writeEnrichedSource(docSrcPath, el)
% Writes the enriched source code for Sphinx documentation generation
%
% Args:
%   docSrcPath (charstring): Base folder for doctest generation, not including trailing path separator
%                            That folder must exist.
%   el (`+replab.+infra.SourceElement`): Either a function or a class
    src = fileread(el.absoluteFilename);
    if isa(el, 'replab.infra.Class')
        lines = strsplit(src, '\n', 'CollapseDelimiters', false);
        % filter comments in property lines
        for i = el.propertyLines
            l = lines{i};
            pos = find(l == '%', 1);
            if ~isempty(pos)
                lines{i} = l(1:pos-1);
            end
        end
        pos = replab.infra.sphinx.findClassCommentEnd(lines);
        tocFN = fullfile(replab.globals.replabPath, 'src', '+replab', '+infra', '+sphinx', 'class.liquid');
        t = replab.lobster.Template.load(tocFN);
        content = t.render(struct('cl', el));
        src = strjoin(horzcat(lines(1:pos-1), content, lines(pos:end)), char(10));
    end
    pkgParts = cellfun(@(s) ['+' s], el.package.packagePath, 'uniform', 0);
    folder = fullfile(docSrcPath, pkgParts{:});
    % Make path
    switch exist(folder)
      case 0
        assert(mkdir(folder), sprintf('Could not create directory %s', folder));
      case 7
        % directory already exists
      otherwise
        error(sprintf('%s already exists but is not a folder', folder));
    end
    filename = fullfile(folder, [el.name '.m']);
    fid = fopen(filename, 'w');
    fwrite(fid, src);
    fclose(fid);
end
