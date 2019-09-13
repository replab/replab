function test_suite = NuRepTest()
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
    youngs = {[1] [2] [1 1] [3] [2 1] [1 1 1] [4] [3 1] [2 2] [2 1 1] [1 1 1 1]};
    for i = 1:length(youngs)
        Y = youngs{i};
        nurep = replab.sym.specht(Y);
        test_suite = replab.RepLaws(nurep).addTestCases(test_suite);
        nurep = replab.sym.youngSeminormal(Y);
        test_suite = replab.RepLaws(nurep).addTestCases(test_suite);
        nurep = replab.sym.youngOrthogonal(Y);
        test_suite = replab.RepLaws(nurep).addTestCases(test_suite);
    end
end
