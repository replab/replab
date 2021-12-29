function sphinx_doc(sphinxFolder, targetFolder, webBaseAddress, preferredInvFile)
% Preprocesses the Sphinx documentation folder
%
% Makes a clean copy of the Sphinx doc folder into the target folder and
% updates all links in the matlab jupyter files according to the API
% described in the file ``objects.inv``. References are recognized
% according to the usual (single) backtick syntax. In case no preferred
% ``objects.inv`` file is provided, or if it is not found, the inventory is
% downloaded from the web API.
%
% This function erases any file in ``targetFolder``.
%
% Args:
%   sphinxFolder (charstring): Location of the initial Sphinx folder
%   targetFolder (charstring): Where to store the modified Sphinx folder
%   webBaseAddress (charstring): base address of the API website, should
%     contain the file ``objects.inv``
%   preferredInvFile (optional, charstring): location of the inventory file
%     to be used instead of the online one (if applicable)

    logFun = @(str) disp(str);

    % Clean the target folder
    [path, directory] = fileparts(targetFolder);
    replab.infra.mkCleanDir(path, directory, logFun);

    % Create a copy of the Sphinx source folder
    copyfile(fullfile(sphinxFolder, '*'), targetFolder);

    % Load the conversion table to create API links in matlab files
    if (nargin >= 4) && exist(preferredInvFile, 'file')
        if unix(['python3 -m sphinx.ext.intersphinx ', preferredInvFile, ' > ', targetFolder, '/API_links.txt'])
            warning('API conversion table found but cound not be extracted, cross-links will not work in .m files');
        else
            logFun(['Using API referrences from ', preferredInvFile]);
        end
    else
        if (nargin >= 4)
            warning(['File ', preferredInvFile, ' not found (to be copied manually), API links will be based on current online inventory at ', webBaseAddress]);
        end
        if unix(['python3 -m sphinx.ext.intersphinx ', webBaseAddress, '/objects.inv > ', targetFolder, '/API_links.txt'])
            warning('API conversion table not found, cross-links will not work in .m files');
        end
    end

    % Select all matlab files and compute their links
    logFun('Updating API links');
    [status, fileList] = unix(['find ', targetFolder, ' -type f | grep [.]m$']);
    fileList = regexp(fileList, '\n', 'split');
    fileList = fileList(1:end-1);
    if status == 0
        pb = replab.infra.repl.ProgressBar(length(fileList));
        for i = 1:length(fileList)
            pb.step(i, fileList{i});
            content = replab.infra.CodeTokens.fromFile(fileList{i});
            lines = content.lines;
            lines = lines(1:find(cellfun(@(x) length(x), lines) ~= 0, 1, 'last'));
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
                        [status, match] = unix(['cat ', targetFolder, '/API_links.txt | grep "^[[:space:]]*', token, '\ "']);
                        if status == 0
                            if length(regexp(match(1:end-1), '\n', 'split')) > 1
                                warning(['Multiple references were found for ', token, ' in the API: ', match]);
                            end
                            link = regexp(match(1:end-1), '\ ([^\ ]+)$', 'tokens');
                            link = [webBaseAddress, '/', link{1}{1}];
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
end
