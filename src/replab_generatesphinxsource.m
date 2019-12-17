function replab_generatesphinxsource
    folderName = 'tmp_sphinxsrc';
    [srcRoot, name, ~] = fileparts(mfilename('fullpath'));
    [root, ~] = fileparts(srcRoot);
    docsrcRoot = fullfile(root, folderName);

    %% Prepare test directory structure
    switch exist(docsrcRoot)
      case 7
        disp('docsrc directory exists, removing it');
        if replab.platformIsOctave
            confirm_recursive_rmdir (false, 'local');
        end
        rmdir(docsrcRoot, 's');
      case 0
        disp('Docsrc directory does not exist yet');
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
