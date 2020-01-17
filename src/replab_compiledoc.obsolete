function replab_compiledoc
% Creates pages for the microsite from the .m files listed in the folder 'srcdocs'
%
% Script files with a name starting with 'MenuPositionXXX_' will
% appear on the main menu of the documentation website at position 'XXX'.
    
    % Make sure we are in the right folder
    initialPath = pwd;
    [pathStr, name, extension] = fileparts(which(mfilename));
    pathStr = strrep(pathStr, '\', '/');
    cd(pathStr)
    cd ../docs_src

    try
        % Check the presence of a SDP solver
        decentSDPSolverInPath = false;
        try
            x = sdpvar(2);
            F = [x >= 0, trace(x) == 1];
            [interfacedata,recoverdata,solver,diagnostic] = compileinterfacedata(F, [], [], [], sdpsettings, 0, 0);
            decentSDPSolverInPath = isempty(diagnostic);
            % If LMILAB was identified as the best solver to solve the
            % problem, this means that no good solver was found.
            if ~isempty(strfind(upper(solver.tag), 'LMILAB'))
                decentSDPSolverInPath = false;
            end
        catch
        end
        if ~decentSDPSolverInPath
            error('No SDP working SDP solver found. Did you run ''install_sdpt3''?');
        end

        % delete existing output and tmp folders
        status = rmdir('tmp','s');
        status = mkdir('tmp');

        % List all available scripts
        scripts = dir('*.m');
        scripts = {scripts.name};

        % Extract the position number of each script. A position of 0 -- or no
        % position -- means the item doesn't show up in the menu
        titles = cell(1,length(scripts));
        menuPositions = zeros(1,length(scripts));
        menuSubsections = cell(1,length(scripts));
        co = 0;
        finalScripts = scripts;
        for file = scripts
            co = co + 1;
            finalScripts{co} = file{1};
            
            % Extract page title
            split = regexp(finalScripts{co}, '(Title[^_]*_)', 'split');
            tokens = regexp(finalScripts{co}, '(Title[^_]*_)', 'tokens');
            if ~isempty(tokens)
                titles{co} = tokens{1}{1}(6:end-1);
            end
            finalScripts{co} = [split{:}];
            
            % Extract menu position
            split = regexp(finalScripts{co}, '(MenuPosition[^_]*_)', 'split');
            tokens = regexp(finalScripts{co}, '(MenuPosition[^_]*_)', 'tokens');
            if ~isempty(tokens)
                number = str2num(tokens{1}{1}(13:end-1));
                menuPositions(co) = floor(number);
            end
            finalScripts{co} = [split{:}];
            
            % Extract menu subsection
            split = regexp(finalScripts{co}, '(MenuSubsection[^_]*_)', 'split');
            tokens = regexp(finalScripts{co}, '(MenuSubsection[^_]*_)', 'tokens');
            if ~isempty(tokens)
                menuSubsections{co} = tokens{1}{1}(15:end-1);
            end
            finalScripts{co} = [split{:}];
            
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

                % If needed, we add the menu subsection
                if ~isempty(menuSubsections{co}) && (length(tline) == 15) && isequal(tline(1:14), 'comments: true')
                    fprintf(fidOut, '%s\n', ['subsection: ', menuSubsections{co}]);
                end

                % add the page title
                % if we find a page title, we remember it
                if (length(tline) >= 7) && isequal(tline(1:7), 'title: ')
                    title = erase(tline(8:end), char(10));
                    if ~isempty(titles{co})
                        % We replace the title with the desired one
                        tline = sprintf('title: %s\n', titles{co});
                    end
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
    catch
    end
    
    % go back to the initial path
    cd(initialPath);
end
