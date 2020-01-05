function replab_generatesphinxsource
% Creates a copy of the source code with the Sphinx references fixed up
    folderName = 'tmp_sphinxsrc';
    rp = replab.settings.replabPath;
    srcRoot = fullfile(rp, 'src');
    docsrcRoot = fullfile(rp, folderName);

    %% Prepare test directory structure
    switch exist(docsrcRoot)
      case 7
        disp('Temp source directory exists, removing it');
        replab.infra.rmdirRec(docsrcRoot);
      case 0
        disp('Temp source directory does not exist yet');
        % Create subfolder if inexistent
        assert(mkdir(docsrcRoot), sprintf('Could not create directory %s', docsrcRoot));
      otherwise
        error('Unknown type')
    end

    disp('Crawling code base');
    cb = replab.infra.crawl(srcRoot);

    disp('Generating rich source code');
    af = cb.allFunctions;
    for i = 1:length(af)
        disp(af{i}.absoluteFilename)
        replab.infra.sphinx.writeEnrichedSource(docsrcRoot, af{i});
    end
    ac = cb.allClasses;
    for i = 1:length(ac)
        disp(ac{i}.absoluteFilename)
        replab.infra.sphinx.writeEnrichedSource(docsrcRoot, ac{i});
    end
end
