function test_suite = DocTestTest()
    disp(['Setting up tests in ', mfilename()]);
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
end

function test_docTestParse
    src = {'>>> 2 + 2', '  ans = 4'};
    ps = replab.infra.doctests.ParseState.fromDocTestBlock(src);
    dt = replab.infra.doctests.DocTest.parse(ps);
    assertEqual(dt.lineNumbers, 1);
    assertEqual(dt.commands, {{'2 + 2'}});
    assertEqual(dt.outputs, {{'ans = 4'}});
end

function test_flags
    flags = 'retries(3), matrixRange(1e-10)';
    s = replab.infra.doctests.parseFlags(flags);
    assertEqual(s, struct('retries', {{'3'}}, 'matrixRange', {{'1e-10'}}));
end

function test_parseTestsError
    lines = {'Blablabla' ...
             '' ...
             'Example:' ...
             '  >>> 2 + 2 % doctest: fun(' ...
             '      ans = 4'};
    assertError(@() replab.infra.doctests.parseTests(lines), '*', ...
                'Line 4: Argument list of fun should be closed by )');
end

% $$$ function test_docTestLineNumbers
% $$$     path = fullfile(replab.settings.replabPath, 'tests/infra/sample');
% $$$     c = replab.infra.CodeBase.crawl(path);
% $$$     cl = c.get('testpkg', 'DocSample');
% $$$     m = cl.get('staticMethod');
% $$$     dtcl = replab.infra.doctests.parseDoc(cl.doc);
% $$$     dtm = replab.infra.doctests.parseDoc(m.doc);
% $$$     assertEqual(length(dtcl), 1);
% $$$     assertEqual(dtcl{1}.relativeLineNumbers, [4]);
% $$$     assertEqual(length(dtm), 1);
% $$$     assertEqual(dtm{1}.relativeLineNumbers, [4 6]);
% $$$ end
