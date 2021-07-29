function test_suite = PartitionTest()
    disp(['Setting up tests in ', mfilename()]);
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;

    D = replab.domain.Partitions(10);
    laws = D.laws;
    test_suite = laws.addTestCases(test_suite);
end
