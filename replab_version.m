function [Major Minor Patch] = replab_version
% Returns RepLAB version read, from "replab_version.txt" 
    [pathStr, name, extension] = fileparts(which(mfilename));
    version_filename = fullfile(pathStr, 'replab_version.txt');
    contents = fileread(version_filename);
    mmp = strsplit(contents, '.');
    Major = str2num(mmp{1});
    Minor = str2num(mmp{2});
    Patch = str2num(mmp{3});
end
