function replab_compiledoc
    % Creates pages for the documentation microsite from the matlab files 
    % listed in the folder 'srcdocs'.
    %
    % Script files with a name starting with 'MenuPositionXXX_' will appear
    % on the main menu of the documentation website at position 'XXX'.
    
    % Make sure we are in the right folder
    initialPath = pwd;
    [pathStr, name, extension] = fileparts(which(mfilename));
    cd(pathStr)
    cd docs_src

    % delete existing output and tmp folders
    status = rmdir('tmp','s');
    status = mkdir('tmp');

    % List all available scripts
    scripts = dir('*.m');
    scripts = {scripts.name};

    % Extract the position number of each script. A position of 0 -- or no
    % position -- means the item doesn't show up in the menu
    menuPositions = zeros(1,length(scripts));
    co = 0;
    finalScripts = scripts;
    for file = scripts
        co = co + 1;
        found = strfind(file{1}, 'MenuPosition');
        if ~isempty(found) && (found(1) == 1)
            endMatch = strfind(file{1}, '_');
            if ~isempty(endMatch)
                endMatch = endMatch(1);
                number = str2num(file{1}(13:endMatch-1));
                if ~isempty(number)
                    menuPositions(co) = floor(number);
                    
                    % The file name will change
                    finalScripts{co} = file{1}(endMatch+1:end);
                end
            end
        end
        % We copy the file to the right place
        copyfile(file{1}, ['tmp/', finalScripts{co}])
    end
    scripts = finalScripts;
    cd tmp;
    
    % Create the pages
    command = @(filename) publish(filename, 'stylesheet', '../stylesheet/mxdom2jekyll.xsl');
    cellfun(command, scripts, 'UniformOutput', 0);

    % We copy the files over and make the necessary adjustments
    copyfile('html/*','../../docs/docs/publish/');
    co = 0;
    for file = scripts
        co = co + 1;
        fidIn = fopen(['html/', file{1}(1:length(file{1})-2), '.html']);
        fidOut = fopen(['../../docs/docs/publish/', file{1}(1:length(file{1})-2), '.html'], 'w');
        title = '';
        while 1
            tline = fgets(fidIn);
            if ~ischar(tline)
                break
            end
            
            % If needed, we add the menu number
            if (menuPositions(co) ~= 0) && (length(tline) == 15) && isequal(tline(1:14), 'comments: true')
                fprintf(fidOut, '%s\n', ['position: ', num2str(menuPositions(co))]);
            end
            
            % add the page title
            % if we find a page title, we remember it
            if (length(tline) >= 7) && isequal(tline(1:7), 'title: ')
                title = erase(tline(8:end), char(10));
            end
            % if we find the beginning of the introduction, we add the title to it
            found = strfind(tline, '<!--introduction-->');
            if ~isempty(found) && ~isempty(title)
                found = found(1);
                tline = [tline(1:found+18), '<h1>', title, '</h1>', tline(found+19:end)];
            end
            
            % We remove any occurence of string like that:
            % </pre><p class="footer"><br><a href="http://www.mathworks.com/products/matlab/">Published with MATLAB&reg; R2016a</a><br></p></div>
            found = strfind(tline, '<p class="footer"><br><a href="http://www.mathworks.com/products/matlab/">Published with MATLAB');
            if ~isempty(found)
                found = found(1);
                endOfString = strfind(tline(found+95:end), '</a><br></p>');
                if ~isempty(endOfString)
                    endOfString = found + 94 +endOfString(1) + 11;
                    tline = [tline(1:found-1), tline(endOfString+1:end)];
                end
            end
            fprintf(fidOut, '%s', tline);
        end
        fclose(fidIn);
        fclose(fidOut);
    end

    % Delete temporary folder
    cd ..
    status = rmdir('tmp', 's');

    % go back to the initial path
    cd(initialPath);
end
