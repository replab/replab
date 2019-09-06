function test_suite = geTest()
    try
        yalmip('version');
        try
            test_functions = localfunctions();
        catch
        end
        initTestSuite;
    catch
        warning('Yalmip not found in the path, some tests will be skipped');
        test_suite=MOxUnitTestSuite();
    end
end

function test_cases
    % We do some sanity checks
    generators = {[2 3 4 5 1]};
    matrix = replab.CommutantVar.fromPermutations(generators);
    assert(length(matrix >= 2) == 3);
    assert(length(matrix >= eye(5)/2) == 3);
    assert(length(2 >= matrix) == 3);
    assert(length(eye(5)/2 >= matrix) == 3);
    
    generators = {[1 3 2]};
    matrix = replab.CommutantVar.fromPermutations(generators);
    R = rand(3);
    R = R + R([1 3 2], [1 3 2]);
    assert(length(matrix >= R+R') == 2);
    assert(length(R+R' >= matrix) == 2);

    R = sdpvar(3);
    R = R + R([1 3 2], [1 3 2]);
    assert(length(matrix >= R) == 2);
    assert(length(R >= matrix) == 2);
end

function test_inputs
    matrix = replab.CommutantVar.fromPermutations({[3 2 1]});
    shouldProduceAnError(@(x) matrix >= rand(5));
    shouldProduceAnError(@(x) matrix >= rand(3));
    shouldProduceAnError(@(x) rand(5) >= matrix);
    shouldProduceAnError(@(x) rand(3) >= matrix);
end
