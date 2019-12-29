function test_suite = DocTestTest()
    disp(['Setting up tests in ', mfilename()]);
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
end

function test_docTestLineNumbers
    path = fullfile(replab.settings.replabPath, 'tests/infra/sample');
    c = replab.infra.CodeBase.crawl(path);
    cl = c.get('testpkg', 'DocSample');
    m = cl.get('staticMethod');
    dtcl = replab.infra.doctests.parseDoc(cl.doc);
    dtm = replab.infra.doctests.parseDoc(m.doc);
    assertEqual(length(dtcl), 1);
    assertEqual(dtcl{1}.relativeLineNumbers, [4]);
    assertEqual(length(dtm), 1);
    assertEqual(dtm{1}.relativeLineNumbers, [4 6]);
end
