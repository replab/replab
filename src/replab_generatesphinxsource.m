function replab_generatesphinxsource
    folderName = 'tmp_sphinxsrc';
    [srcRoot, name, ~] = fileparts(mfilename('fullpath'));
    [root, ~] = fileparts(srcRoot);
    docsrcRoot = fullfile(root, folderName);

    %% Prepare test directory structure
    switch exist(docsrcRoot)
      case 7
        disp('./tmp_sphinxsrc directory exists, removing it');
        replab.infra.rmdirRec(docsrcRoot);
      case 0
        disp('./tmp_sphinxsrc directory does not exist yet');
      otherwise
        error('Unknown type')
    end

    % Create subfolder if inexistent
    [success, message, messageid] = mkdir(root, folderName);
    
    disp('Crawling code base');
    codeBase = replab.infra.OldCodeBase.crawl(fullfile(root, 'src'));
    
    disp('Writing tests');
    codeBase.writeEnrichedSource(docsrcRoot);
end
