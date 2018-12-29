function test_suite = DivisionAlgebraTest()
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
    R = replab.DivisionAlgebra.real;
    C = replab.DivisionAlgebra.complex;
    H = replab.DivisionAlgebra.quaternion;
    test_suite = replab.DivisionAlgebraLaws(R).addTestCases(test_suite);
    test_suite = replab.DivisionAlgebraLaws(C).addTestCases(test_suite);
    test_suite = replab.DivisionAlgebraLaws(H).addTestCases(test_suite);
end
