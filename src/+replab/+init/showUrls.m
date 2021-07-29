function showUrls
    deps = {replab.init.cyclolab replab.init.MOcov replab.init.MOxUnit replab.init.sdpt3 replab.init.YALMIP};
    fprintf('\nDownload and place those ZIP files in the external/ folder:\n\n');
    for i = 1:length(deps)
        dep = deps{i};
        fprintf('%s\n', dep.zipUrl);
    end
    fprintf('\nThose files will have the names:\n\n');
    for i = 1:length(deps)
        dep = deps{i};
        fprintf('%s\n', dep.zipFilename);
    end
    fprintf('\n');
end
