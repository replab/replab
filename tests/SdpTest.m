function test_suite = SdpTest()
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
end

function test_simple_structure
    % We do a sanity check with one group
    generators = {[2 3 4 5 1]};
    matrix = replab.SdpRep.fromGenerators(generators);
    matrix = matrix.fullMatrix;
    for i = 1:length(generators)
        difference = matrix - matrix(generators{i}, generators{i});
        vars = getvariables(difference);
        for j = 1:length(vars)
            coeffs = getbasematrix(difference, vars(i));
            assert(norm(coeffs(:)) <= replab.Settings.doubleEigTol);
        end
    end
end
