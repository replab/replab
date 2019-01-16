function [major minor patch snapshot contents] = replab_version(contents)
% Parses a RepLAB version number
%
% If contents is provided, then parses that string; otherwise reads the version
% string from the replab_version.txt file.
    [pathStr, name, extension] = fileparts(which(mfilename));
    version_filename = fullfile(pathStr, 'replab_version.txt');
    contents = fileread(version_filename);
    contents = strrep(contents, sprintf('\n'), '');
    contents = strrep(contents, sprintf('\r'), '');    
    contents = strtrim(contents);
    snapshotString = '-SNAPSHOT';
    if length(contents) > length(snapshotString) && isequal(contents(end-length(snapshotString)+1:end), snapshotString)
        snapshot = true;
        mmp = contents(1:end-length(snapshotString));
    else
        snapshot = false;
        mmp = contents;
    end
    mmp = strsplit(mmp, '.');
    major = str2num(mmp{1});
    minor = str2num(mmp{2});
    patch = str2num(mmp{3});
end
