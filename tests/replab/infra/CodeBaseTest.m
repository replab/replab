function tests = CodeBaseTest()
    disp(['Setting up tests in ', mfilename()]);
    tests = functiontests(localfunctions);
end

function test_crawl(param)
    path = fullfile(replab.settings.replabPath, 'tests/replab/infra/sample');
    c = replab.infra.crawl(path);
    cl = c.get('testpkg').ownClasses;
    cl = cellfun(@(x) x.name, cl, 'uniform', 0);
    assertEqual(cl, {'DocSample' 'Group' 'Monoid'});
end
