function test_suite = DivisionAlgebraTest()
    disp(['Setting up tests in ', mfilename()]);
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
    D = replab.DivisionAlgebra('R->C');
    test_suite = D.laws.addTestCases(test_suite);
    D = replab.DivisionAlgebra('C->R');
    test_suite = D.laws.addTestCases(test_suite);
    D = replab.DivisionAlgebra('H->C');
    test_suite = D.laws.addTestCases(test_suite);
    D = replab.DivisionAlgebra('H->R:equivariant');
    test_suite = D.laws.addTestCases(test_suite);
    D = replab.DivisionAlgebra('H->R:rep');
    test_suite = D.laws.addTestCases(test_suite);
end
