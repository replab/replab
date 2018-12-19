% Here we create all web pages for the documentation website

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


% go back to the initial path
cd(initialPath)
