function result = runpublish
    % Here we create all web pages for the documentation website
    initialPath = pwd;
    [pathStr, name, extension] = fileparts(which(mfilename));
    cd(pathStr)

    % Make sure we are in the right folder
    initialPath = pwd;
    [pathStr, name, extension] = fileparts(which(mfilename));
    cd(pathStr)

    % delete existing output folder
    status = rmdir('html','s');

    % List all available scripts
    scripts = dir('*.m');
    scripts = {scripts.name};
    scripts = scripts(find(strcmp(scripts, [name, extension]) == 0));

    % Create the pages
    command = @(filename) publish(filename, 'stylesheet', 'stylesheet/mxdom2jekyll.xsl');
    cellfun(command, scripts, 'UniformOutput', 0);

    % We remove the matlab trademark and copy files over
    copyfile('html/*','../docs/docs/publish/');
    for file = scripts
        fidIn = fopen(['html/', file{1}(1:length(file{1})-2), '.html']);
        fidOut = fopen(['../docs/docs/publish/', file{1}(1:length(file{1})-2), '.html'], 'w');
        while 1
            tline = fgets(fidIn);
            if ~ischar(tline)
                break
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
    rmdir('html', 's');

    % go back to the initial path
    cd(initialPath);
end
