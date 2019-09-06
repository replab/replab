function test_suite = codeCoverageHelperFunction()
    % This file is a helper function for replab_runtests.m. To use it, call
    % 'replab_runtests(1)'.
    % 
    % The name of this function should not contain the string "test"
    % (otherwise, this will trigger an infinite loop of tests).
    
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
    % consists of all tests in the 'tests' folder (this one is not included
    % because its name doesn't contain the string 'test'). This way,
    % calling MOxUnit on this file will do the following:
    % - define a test suite saying that we want to run MOxUnit on all tests
    %   in the 'tests' folder
    % - activate coverage of the source code in the replab package
    % - launch the test suite, which in turn defines the test suite we are
    %   interested in: tests in the 'tests' folder (code coverage is
    %   activated at this point
    % - run the tests
    % - finally, collect the coverage result (if requested)
    %
    % This is what happens when calling replab_runtests with code coverage
    % enabled.
    
    % At this point it is important to remove the replab path so that the
    % default location for the source code will be the one with activated
    % coverage
    [pathStr, name, extension] = fileparts(which(mfilename));
    replabFolder = [pathStr(1:find(pathStr=='/',1,'last')-1), '/replab'];
    rmpath(replabFolder);
    
    % We define a test which consists of all other tests
    f = @() moxunit_runtests('tests', '-recursive', '-verbose', '-junit_xml_file', 'testresults.xml');
    testCase = MOxUnitFunctionHandleTestCase('All tests with coverage', replabFolder, f);
    test_suite = addTest(test_suite, testCase);
    
    % We can restore the path
    addpath(replabFolder);
end


