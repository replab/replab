function dispStubs(pkg)
% Generate the API autodoc stubs
%
% Args:
%   pkg (`+replab.+infra.Package`): Package to enumerate elements in Sphinx documentation
    stubsFN = fullfile(replab.settings.replabPath, 'src', '+replab', '+infra', '+sphinx', 'stubs.liquid');
    t = replab.lobster.Template.load(stubsFN);
    content = t.render(struct('pkg', pkg));
    disp(content);
end
