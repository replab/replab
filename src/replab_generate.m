function replab_generate(what)
% Code generation function
%
% With ``what = 'sphinx'``, this preprocesses the RepLAB source code to complement the
% `Sphinx Matlab domain <https://github.com/sphinx-contrib/matlabdomain>`_ job.
%
% With ``what = 'doctests'``, this extracts the doctests from the source code and write
% them to the doctests code folder.
%
% Args:
%   what ({'sphinx', 'doctests', 'all'}, optional): What to generate, default ``'all'``

    if nargin < 1
        what = 'all';
    end

    logFun = @(str) disp(str);

    assert(ismember(what, {'sphinx' 'doctests' 'all'}), 'Argument must be one of: ''sphinx'', ''doctests'', ''all''');

    rp = replab.settings.replabPath;
    srcRoot = fullfile(rp, 'src');

    disp('Crawling code base');
    cb = replab.infra.crawl(srcRoot);

    if isequal(what, 'sphinx') || isequal(what, 'all')
        % Generate Sphinx preprocessed source files
        srcRoot = fullfile(rp, 'src');
        sphinxRoot = fullfile(rp, 'sphinxdocs');
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
