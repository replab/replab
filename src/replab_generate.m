function result = replab_generate(what)
% Code generation function
%
% Depending on the value of the argument ``what``:
%
% - ``all`` regenerates all generated code/documentation
%
% - ``clear`` clears out all directories with autogenerated code/doc
%
% - ``sphinxsrc`` preprocesses the RepLAB source code to complement the
%   `Sphinx Matlab domain <https://github.com/sphinx-contrib/matlabdomain>`_ job.
%   It generates a table of contents for all source code files, generates method and
%   property occurences for inherited members, shifts the documentation of properties to
%   be explicit (to support multiline documentation and types).
%   It also copies the source files at the root of the source folder in a ``root`` subdirectory
%   so that the Sphinx Matlab domain has a "module" name for those.
%
% - ``sphinxbuild`` runs the Sphinx documentation generation.
%
% - ``sphinx`` runs all the Sphinx generation steps.
%
% - ``doctests`` extracts the doctests from the source code and writes them to the doctests code folder.
%
% - ``notebooks`` extracts the jupyter notebooks from the doc and write them to the notebooks code folder
%
% Args:
%   what ({'clear', 'sphinx*', 'sphinx', 'doctests', 'notebooks', 'all'}, optional): What to generate, default ``'all'``
%
% Results:
%     logical: True unless an error was detected

    if nargin < 1
        what = 'all';
    end
    
    result = true;
    
    % Make sure we are in the correct path
    initialPath = pwd;
    [pathStr, name, extension] = fileparts(which(mfilename));
    pathStr = strrep(pathStr, '\', '/');
    cd(pathStr)
    cd ..

    logFun = @(str) disp(str);
    valid = {'clear' 'sphinx' 'sphinxbuild' 'sphinxsrc' 'doctests' 'notebooks' 'all'};
    validStr = strjoin(cellfun(@(x) sprintf('''%s''', x), valid, 'uniform', 0), ', ');
    assert(ismember(what, valid), 'Argument must be one of: %s', validStr);

    rp = replab.globals.replabPath;
    srcRoot = fullfile(rp, 'src');

    if isequal(what, 'sphinxsrc') || isequal(what, 'sphinx') || isequal(what, 'doctests') || isequal(what, 'all')
        disp('Crawling code base');
        cb = replab.infra.crawl(srcRoot);
    end

    if isequal(what, 'sphinxsrc') || isequal(what, 'sphinx') || isequal(what, 'all') || isequal(what, 'clear')
        % Generate Sphinx preprocessed source files
        sphinxRoot = fullfile(rp, 'sphinx');
        sphinxSrcRoot = fullfile(sphinxRoot, '_src');
        replab.infra.mkCleanDir(sphinxRoot, '_src', logFun);
        if ~isequal(what, 'clear')
            logFun('Generating rich source code');
            els = cb.allSourceElements;
            pb = replab.infra.repl.ProgressBar(length(els));
            for i = 1:length(els)
                pb.step(i, els{i}.fullIdentifier);
                replab.infra.sphinx.writeEnrichedSource(sphinxSrcRoot, els{i});
            end
            pb.finish;
            % Copy root files from the root source folder to a subfolder named 'root'
            mkdir(sphinxSrcRoot, 'root');
            files = dir([sphinxSrcRoot filesep '*.m']);
            for i = 1:length(files)
                assert(~files(i).isdir, 'Files ending in .m cannot be directories');
                name = files(i).name;
                copyfile(fullfile(sphinxSrcRoot, name), fullfile(sphinxSrcRoot, 'root', name));
            end
            files = dir([rp filesep '*.m']);
            for i = 1:length(files)
                assert(~files(i).isdir, 'Files ending in .m cannot be directories');
                name = files(i).name;
                copyfile(fullfile(rp, name), fullfile(sphinxSrcRoot, 'root', name));
            end
            tocFile = fullfile(sphinxRoot, '_src', 'index.rst');
            replab.infra.sphinx.writeTOC(tocFile, cb);
        end
    end

    if isequal(what, 'sphinxbuild') || isequal(what, 'sphinx') || isequal(what, 'all')
        if ~exist(fullfile(rp, '_sphinx'))
            mkdir(rp, '_sphinx');
        end
        replab.infra.cleanDir(fullfile(rp, '_sphinx'), {'.git'});
        if ~exist(fullfile(rp, 'docs'))
            mkdir(rp, 'docs');
        end
        replab.infra.cleanDir(fullfile(rp, 'docs'), {'.git'});
        if ~isequal(what, 'clear')
            % Create a modifiable copy of the sphinx folder
            copyfile(fullfile(rp, 'sphinx/*'), fullfile(rp, '_sphinx'))

            % Load the conversion table to create API links in matlab files
            baseWeb = 'https://replab.github.io/replab';
            if exist([rp, '/objects.inv'], 'file')
                if unix(['python3 -m sphinx.ext.intersphinx ', rp, '/docs.objects.inv > ', rp, '/_sphinx/API_links.txt'])
                    warning('API conversion table not found, cross-links will not work in .m files');
                end
            else
                warning('No local objects.inv file found (to be copied manually into the sphinx folder), API links will be based on current online doc');
                if unix(['python3 -m sphinx.ext.intersphinx ', baseWeb, '/objects.inv > ', rp, '/_sphinx/API_links.txt'])
                    warning('API conversion table not found, cross-links will not work in .m files');
                end
            end

            % Select matlab files in the doc (excluding API) and compute
            % their links
            [status, fileList] = unix(['find ', rp, '/_sphinx -type f | grep -v ^"', rp, '/_sphinx/_src" | grep [.]m$']);
            fileList = regexp(fileList, '\n', 'split');
            fileList = fileList(1:end-1);
            if status == 0
                pb = replab.infra.repl.ProgressBar(length(fileList));
                for i = 1:length(fileList)
                    pb.step(i, fileList{i});
                    content = replab.infra.CodeTokens.fromFile(fileList{i});
                    lines = content.lines;
                    for j = find(content.tags == '%')
                        line = lines{j};
                        extents = regexp(line, '(`~?\+replab\.[\w,\.]+`|`~?root\.[\w,\.]+`)', 'tokenExtents');
                        extents{end+1} = length(line)+1;
                        if ~isempty(extents)
                            newLine = line(1:extents{1}(1)-1);
                            for k = 1:length(extents)-1
                                token = line(extents{k}(1)+1:extents{k}(2)-1);
                                silent = (token(1) == '~');
                                if silent
                                    tokenName = regexp(token, '\.*(\w+)$', 'tokens');
                                    tokenName = tokenName{1}{1};
                                    token = token(2:end);
                                else
                                    tokenName = token;
                                end
                                tokenName = tokenName(tokenName ~= '+');
                                [status, match] = unix(['cat ', rp, '/_sphinx/API_links.txt | grep "^[[:space:]]*', token, '\ "']);
                                if status == 0
                                    if length(regexp(match(1:end-1), '\n', 'split')) > 1
                                        warning(['Multiple references were found for ', token, ' in the API: ', match]);
                                    end
                                    link = regexp(match(1:end-1), '\ ([^\ ]+)$', 'tokens');
                                    link = [baseWeb, '/', link{1}{1}];
                                    newLine = [newLine, '[', tokenName, '](', link,')'];
                                else
                                    warning(['Reference ', token, ' in ', fileList{i}, ' was not found in the API.']);
                                    newLine = [newLine, token];
                                end
                                newLine = [newLine, line(extents{k}(2)+1:extents{k+1}(1)-1)];
                            end
                            lines{j} = newLine;
                        end
                    end
                    fid = fopen(fileList{i},'w');
                    for j = 1:length(lines)
                        fprintf(fid, '%s\n', lines{j});
                    end
                    fclose(fid);
                end
                pb.finish;
            end

            % Now launch sphinx in modified folder
            disp('Running Sphinx');
            lastPath = pwd;
            cd(rp);
            cmd = 'sphinx-build -b html _sphinx docs';
            disp(['Running ' cmd]);
            status = system(cmd);
            if status ~= 0
                result = false;
            end
            cd(lastPath);
        end
    end

    if isequal(what, 'doctests') || isequal(what, 'all') || isequal(what, 'clear')
        % Generate doctests
        testRoot = fullfile(rp, 'tests');
        doctestRoot = fullfile(rp, 'tests', 'doctests');
        replab.infra.mkCleanDir(testRoot, 'doctests');
        if ~isequal(what, 'clear')
            logFun('Generating doctests');
            els = cb.allSourceElements;
            pb = replab.infra.repl.ProgressBar(length(els));
            for i = 1:length(els)
                pb.step(i, els{i}.fullIdentifier);
                replab.infra.doctests.writeElementDocTests(doctestRoot, els{i});
            end
            pb.finish;
        end
    end

    if isequal(what, 'notebooks') || isequal(what, 'all') || isequal(what, 'clear')
        % Generate a copy of the jupyter notebooks
        testRoot = fullfile(rp, 'tests');
        notebooksRoot = fullfile(rp, 'tests', 'notebooks');
        replab.infra.mkCleanDir(testRoot, 'notebooks');
        if ~isequal(what, 'clear')
            logFun('Copying jupyter notebooks');
            els = replab.infra.notebooks.listNotebooks;
            pb = replab.infra.repl.ProgressBar(length(els));
            for i = 1:length(els)
                pb.step(i, els{i,3});
                replab.infra.notebooks.writeNotebook(notebooksRoot, els(i,:));
            end
            pb.finish;
        end
    end

    % return to the previous path
    cd(initialPath);
end
