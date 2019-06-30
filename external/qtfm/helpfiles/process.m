% Script to process XML help files into HTML. This script can be invoked
% from the Matlab Start button (via the function qftm_helpup).
%
% Copyright (c) 2008, 2016 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

Files = dir('xmlfiles/*.xml'); % Entries in the name field of this struct
                               % will be of the form 'functionname.xml'.

N = length(Files);

S = 'xmlfiles/qtfmfunction.xsl';

h = waitbar(0, 'Processing XML to HTML helpfiles ...');

for i = 1:N
    waitbar(i/N, h)
    
    F = Files(i).name;
    
    % end-4 strips the characters '.xml' from the end of the filename.
    xslt(['xmlfiles/', F], S, [F(1:end-4), '.html']);
end

close(h)

% $Id: process.m 1004 2017-11-15 17:14:09Z sangwine $
