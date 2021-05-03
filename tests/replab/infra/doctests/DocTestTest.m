function test_suite = DocTestTest()
    disp(['Setting up tests in ', mfilename()]);
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
end

function test_docTestParseStatement
    str = {'>>> 2 + ...' '    2' '    4' '>>> 3 + 3' '    9'};
    dtt = replab.infra.doctests.DocTestTokens.lex(str);
    assertEqual(dtt.tags, 'Sooso$');
    [pos, dts] = replab.infra.doctests.DocTestStatement.parse(dtt, 1, @(ln) 1);
    assertEqual(pos, 4)
    assertEqual(dts.command, '2 + 2');
end

function test_docTestParse
    src = {'>>> 2 + 2', '  ans = 4'};
    dtt = replab.infra.doctests.DocTestTokens.lex(src);
    dt = replab.infra.doctests.DocTest.parse1(dtt, @(ln) 1);
    assert(dt.nStatements == 1);
    st = dt.statements{1};
    assertEqual(st.lineNumber, 1);
    assertEqual(st.command, '2 + 2');
    assertEqual(st.output, {'ans = 4'});
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
% $$$     path = fullfile(replab.globals.replabPath, 'tests/infra/sample');
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
