function test_suite = geTest()
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
end

function test_cases
    % We do some sanity checks
    matrix = replab.CommutantVar.fromPermutations({[2 3 4 5 1]});
    assert(length(matrix >= 2) == 3);
    assert(length(matrix >= eye(5)/2) == 3);
    assert(length(2 >= matrix) == 3);
    assert(length(eye(5)/2 >= matrix) == 3);
    
    matrix = replab.CommutantVar.fromPermutations({[1 3 2]});
    R = rand(3);
    R = R + R([1 3 2], [1 3 2]);
    assert(length(matrix >= R+R') == 2);
    assert(length(R+R' >= matrix) == 2);

    R = sdpvar(3);
    R = R + R([1 3 2], [1 3 2]);
    assert(length(matrix >= R) == 2);
    %assert(length(R >= matrix) == 2); % We cannot check due to a bug in octave
end

function with_linear_constraint
    matrix = replab.CommutantVar.fromSdpMatrix(sdpvar(5,5,'hankel'), {[2 3 4 5 1]});
    assert(length(matrix >= 2) == 4);
    assert(length(matrix >= eye(5)/2) == 4);
    assert(length(2 >= matrix) == 4);
    assert(length(eye(5)/2 >= matrix) == 4);
    
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
    matrix = replab.CommutantVar.fromPermutations({[3 2 1]});
    shouldProduceAnError(@(x) matrix >= rand(5));
    shouldProduceAnError(@(x) matrix >= rand(3));
    shouldProduceAnError(@(x) rand(5) >= matrix);
    shouldProduceAnError(@(x) rand(3) >= matrix);
end
