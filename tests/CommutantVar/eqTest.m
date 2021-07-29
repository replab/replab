function test_suite = eqTest()
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
end

function test_general
    global matrix231 matrix23451 matrix23451H
    matrix = matrix23451;
    assert(length(matrix == 0) == 3);

    matrix = matrix231;
    R = rand(3);
    R = R + R([2 3 1], [2 3 1]) + R([3 1 2], [3 1 2]);
    assert(length(matrix == R+R') == 2);
    assert(length(R+R' == matrix) == 2);

    R = sdpvar(3);
    R = R + R([2 3 1], [2 3 1]) + R([3 1 2], [3 1 2]);
    assert(length(matrix == R) == 2);
    %assert(length(R == matrix) == 2); % We cannot check due to a bug in octave
end

function test_with_linear_constraint
    global matrix231 matrix23451 matrix23451H
    matrix = matrix23451H;
    assert(length(matrix == 0) == 4);

    if ReplabTestParameters.onlyFastTests
        return;
    end

    matrix = replab.CommutantVar.fromSdpMatrix(sdpvar(3,3,'hankel'), {[1 3 2]});
    R = rand(3);
    R = R + R([1 3 2], [1 3 2]);
    assert(length(matrix == R+R') == 3);
    assert(length(R+R' == matrix) == 3);

    R = sdpvar(3);
    R = R + R([1 3 2], [1 3 2]);
    assert(length(matrix == R) == 3);
    %assert(length(R == matrix) == 3); % We cannot check due to a bug in octave
end

function test_double_combination
    if ReplabTestParameters.onlyFastTests
        return;
    end

    matrix1 = replab.CommutantVar.fromPermutations({[2 3 4 5 1]}, 'symmetric', 'real');
    matrix2 = replab.CommutantVar.fromSdpMatrix(matrix1.fullMatrix, {[2 1 3 4 5]});

    obj = matrix2(1,5)-matrix2(1,4);
    optimize([matrix2 >= 0], obj, sdpsettings('verbose', 0));
    assert(abs(value(obj)) < replab.globals.doubleSdpTol);
    optimize([matrix2 >= 0], -obj, sdpsettings('verbose', 0));
    assert(abs(value(obj)) < replab.globals.doubleSdpTol);
end

function test_inputs
    global matrix231 matrix23451 matrix23451H
    matrix = matrix231;
    shouldProduceAnError(@(x) matrix == 1);
end
