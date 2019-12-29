function test_suite = CodeBaseTest()
    disp(['Setting up tests in ', mfilename()]);
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
end

function test_crawl
    path = fullfile(replab.settings.replabPath, 'tests/infra/sample');
    c = replab.infra.CodeBase.crawl(path);
    cl = c.get('testpkg').ownClasses;
    cl = cellfun(@(x) x.name, cl, 'uniform', 0);
    assertEqual(cl, {'DocSample' 'Group' 'Monoid'});
end
