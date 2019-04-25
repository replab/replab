function test_suite = test()
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
    
    % MOxUnit starts monitoring source files for coverage only after the
    % test suite has been defined. However, tests in the folder 'tests'
    % typically use objects in the package +replab to define their tests.
    % This means that the tests are defined with some code which is not
    % monitored by MOxUnit.
    % 
    % This file avoids this behavior by defining a test suite which
    % consists of all tests in the 'tests' folder. This way, calling
    % MOxUnit on this file will do the following:
    % - define a test suite saying that we want to run MOxUnit on all tests
    %   in the 'tests' folder
    % - activate coverage of the source code in the replab package
    % - launch the test suite, which in turn defines the test suite we are
    %   interested in: tests in the 'tests' folder (code coverage is
    %   activated at this point
    % - run the tests
    % - finally, collect the coverage result (if requested)
    %
    % Here is a way to call this file to achieve this:
    % moxunit_runtests('test.m', '-verbose', '-with_coverage', '-cover', '.', '-cover_exclude', 'external', '-cover_json_file', 'coverage.json', '-cover_xml_file', 'coverage.xml', '-cover_html_dir', 'coverage_html')
    
    % At this point it is important to remove the current path so that the
    % default location for the source code will be the one with activated
    % coverage
    [pathStr, name, extension] = fileparts(which(mfilename));
    rmpath(pathStr);
    
    % We define a test which consists of all other tests
    f = @() moxunit_runtests('tests', '-recursive', '-verbose', '-junit_xml_file', 'testresults.xml');
    testCase = MOxUnitFunctionHandleTestCase('All tests with coverage', pwd, f);
    test_suite = addTest(test_suite, testCase);
    
    % We can restore the path
    addpath(pathStr);
end


