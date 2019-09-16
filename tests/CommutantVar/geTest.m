function test_suite = geTest()
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
end

function test_general
    global matrix231 matrix23451 matrix23451H
    matrix = matrix23451;
    assert(length(matrix >= 2) == 3);
    assert(length(matrix >= eye(5)/2) == 3);
    assert(length(2 >= matrix) == 3);
    assert(length(eye(5)/2 >= matrix) == 3);
    
    matrix = matrix231;
    R = rand(3);
    R = R + R([2 3 1], [2 3 1]) + R([3 1 2], [3 1 2]);
    assert(length(matrix >= R+R') == 2);
    assert(length(R+R' >= matrix) == 2);

    R = sdpvar(3);
    R = R + R([2 3 1], [2 3 1]) + R([3 1 2], [3 1 2]);
    assert(length(matrix >= R) == 2);
    %assert(length(R >= matrix) == 2); % We cannot check due to a bug in octave
end

function test_with_linear_constraint
    global matrix231 matrix23451 matrix23451H
    matrix = matrix23451H;
    assert(length(matrix >= 2) == 4);
    assert(length(matrix >= eye(5)/2) == 4);
    assert(length(2 >= matrix) == 4);
    assert(length(eye(5)/2 >= matrix) == 4);
    
    if ReplabTestParameters.onlyFastTests
        return;
    end
    
    matrix = replab.CommutantVar.fromSdpMatrix(sdpvar(3,3,'hankel'), {[1 3 2]});
    R = rand(3);
    R = R + R([1 3 2], [1 3 2]);
    assert(length(matrix >= R+R') == 3);
    assert(length(R+R' >= matrix) == 3);

    R = sdpvar(3);
    R = R + R([1 3 2], [1 3 2]);
    assert(length(matrix >= R) == 3);
    %assert(length(R >= matrix) == 3); % We cannot check due to a bug in octave
end

function test_inputs
    global matrix231 matrix23451 matrix23451H
    matrix = matrix231;
    shouldProduceAnError(@(x) matrix >= rand(5));
    shouldProduceAnError(@(x) matrix >= rand(3));
    shouldProduceAnError(@(x) rand(5) >= matrix);
    shouldProduceAnError(@(x) rand(3) >= matrix);
end
