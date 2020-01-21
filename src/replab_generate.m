function replab_generate(what)
% Code generation function
%
% With ``what = 'sphinxsrc'``, this preprocesses the RepLAB source code to complement the
% `Sphinx Matlab domain <https://github.com/sphinx-contrib/matlabdomain>`_ job.
%
% With ``what = 'sphinxjupyter'`` this processes the Jupyter notebooks in the Sphinx directory.
%
% With ``what = 'sphinxbuild'``, this runs the Sphinx documentation generation.
%
% With ``what = 'sphinx'``, this runs all the Sphinx generation steps.
%
% With ``what = 'doctests'``, this extracts the doctests from the source code and write
% them to the doctests code folder.
%
% Args:
%   what ({'sphinx*', 'sphinx', 'doctests', 'all'}, optional): What to generate, default ``'all'``

    if nargin < 1
        what = 'all';
    end

    logFun = @(str) disp(str);
    valid = {'sphinx' 'sphinxbuild' 'sphinxsrc' 'sphinxjupyter' 'doctests' 'all'};
    validStr = strjoin(cellfun(@(x) sprintf('''%s''', x), valid, 'uniform', 0), ', ');
    assert(ismember(what, valid), 'Argument must be one of: %s', validStr);

    rp = replab.settings.replabPath;
    srcRoot = fullfile(rp, 'src');

    disp('Crawling code base');
    cb = replab.infra.crawl(srcRoot);

    if isequal(what, 'sphinxjupyter') || isequal(what, 'sphinx') || isequal(what, 'all')
        toExplore = {{'sphinx'}};
        notebooks = {};
        sourceSuffix = '_source.ipynb';
        targetSuffix = '.ipynb';
        while length(toExplore) > 0
            subpath = toExplore{1};
            toExplore = toExplore(2:end);
            path = fullfile(rp, subpath{:});
            children = dir(path);
            for i = 1:length(children)
                name = children(i).name;
                if isequal(name, '.') || isequal(name, '..')
                    % do nothing
                elseif children(i).isdir
                    % new folder
                    toExplore{1,end+1} = horzcat(subpath, {name});
                else
                    if replab.compat.endsWith(name, sourceSuffix)
                        notebooks{1,end+1} = strjoin(horzcat(subpath, {name}), filesep);
                    end
                end
            end
        end
        lastPath = pwd;
        for i = 1:length(notebooks)
            source = notebooks{i};
            [~, sourceName, sourceExt] = fileparts(source);
            target = [sourceName sourceExt];
            target = [target(1:end-length(sourceSuffix)) targetSuffix];
            cmd = sprintf('jupyter nbconvert --to notebook --execute "%s" --output="%s"', source, target);
            disp(['Running ' cmd]);
            system(cmd);
        end
        cd(lastPath);
    end

    if isequal(what, 'sphinxsrc') || isequal(what, 'sphinx') || isequal(what, 'all')
        % Generate Sphinx preprocessed source files
        srcRoot = fullfile(rp, 'src');
        sphinxRoot = fullfile(rp, 'sphinx');
        sphinxSrcRoot = fullfile(sphinxRoot, '_src');
        replab.infra.mkCleanDir(sphinxRoot, '_src', logFun);
        logFun('Generating rich source code');
        els = cb.allSourceElements;
        pb = replab.infra.repl.ProgressBar(length(els));
        for i = 1:length(els)
            pb.step(i, els{i}.fullIdentifier);
            replab.infra.sphinx.writeEnrichedSource(sphinxSrcRoot, els{i});
        end
        pb.finish;
    end

    if isequal(what, 'sphinxbuild') || isequal(what, 'sphinx') || isequal(what, 'all')
        replab.infra.mkCleanDir(rp, 'docs', logFun);
        disp('Running Sphinx');
        lastPath = pwd;
        system('sphinx-build -b html sphinx docs');
        cd(lastPath);
    end

    if isequal(what, 'doctests') || isequal(what, 'all')
        % Generate doctests
        testRoot = fullfile(rp, 'tests');
        doctestRoot = fullfile(rp, 'tests', 'doctest');
        replab.infra.mkCleanDir(testRoot, 'doctest');
        logFun('Generating doctests');
        els = cb.allSourceElements;
        pb = replab.infra.repl.ProgressBar(length(els));
        for i = 1:length(els)
            pb.step(i, els{i}.fullIdentifier);
            replab.infra.doctests.writeElementDocTests(doctestRoot, els{i});
        end
        pb.finish;
    end

end
