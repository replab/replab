function test_suite = DivisionAlgebraTest()
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
    n = 10;
    R = replab.rep.DivisionAlgebra.real;
    C = replab.rep.DivisionAlgebra.complex;
    H = replab.rep.DivisionAlgebra.quaternion;
    test_suite = R.lawsAddTestCases(test_suite);
    test_suite = C.lawsAddTestCases(test_suite);
    test_suite = H.lawsAddTestCases(test_suite);
end
