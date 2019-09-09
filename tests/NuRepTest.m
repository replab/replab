function test_suite = NuRepTest()
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
    youngs = {[1] [2] [1 1] [3] [2 1] [1 1 1] [4] [3 1] [2 2] [2 1 1] [1 1 1 1]};
    for i = 1:length(youngs)
        Y = youngs{i};
        nurep = replab.nu.specht(Y);
        test_suite = replab.nu.RepLaws(nurep).addTestCases(test_suite);
        nurep = replab.nu.youngSeminormal(Y);
        test_suite = replab.nu.RepLaws(nurep).addTestCases(test_suite);
        nurep = replab.nu.youngOrthogonal(Y);
        test_suite = replab.nu.RepLaws(nurep).addTestCases(test_suite);
    end
end
