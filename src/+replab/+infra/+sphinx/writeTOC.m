function writeTOC(filename, codeBase)
% Writes the list of all functions and classes in a ReST file
%
% Args:
%   filename (charstring): Filename to write
%   codeBase (`+replab.+infra.CodeBase`): Code base to enumerate
    tocFN = fullfile(replab.settings.replabPath, 'src', '+replab', '+infra', '+sphinx', 'toc.liquid');
    t = replab.lobster.Template.load(tocFN);
    content = t.render(struct('codeBase', codeBase));
    fid = fopen(filename, 'w');
    fwrite(fid, content);
    fclose(fid);
end
