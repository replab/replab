function replab_compiledoc
    % Here we create all web pages for the documentation website

    % Make sure we are in the right folder
    initialPath = pwd;
    [pathStr, name, extension] = fileparts(which(mfilename));
    cd(pathStr)
    cd srcdocs

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

    % We remove the matlab trademark and copy files over
    copyfile('html/*','../../docs/docs/publish/');
    co = 0;
    for file = scripts
        co = co + 1;
        fidIn = fopen(['html/', file{1}(1:length(file{1})-2), '.html']);
        fidOut = fopen(['../../docs/docs/publish/', file{1}(1:length(file{1})-2), '.html'], 'w');
        while 1
            tline = fgets(fidIn);
            if ~ischar(tline)
                break
            end
            
            % If needed, we add the menu number
            if (menuPositions(co) ~= 0) && (length(tline) == 15) && isequal(tline(1:14), 'comments: true')
                fprintf(fidOut, '%s\n', ['position: ', num2str(menuPositions(co))]);
            end
            
            % We remove any occurence of string like that:
            % </pre><p class="footer"><br><a href="http://www.mathworks.com/products/matlab/">Published with MATLAB&reg; R2016a</a><br></p></div>
            found = strfind(tline, '</pre><p class="footer"><br><a href="http://www.mathworks.com/products/matlab/">Published with MATLAB');
            if ~isempty(found)
                foud = found(1);
                endOfString = strfind(tline(found+101:end), '</a><br></p></div>');
                if ~isempty(endOfString)
                    endOfString = found + 100 +endOfString(1) + 17;
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
