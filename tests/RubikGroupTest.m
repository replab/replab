function test_suite = RubikGroupTest()
    disp(['Setting up tests in ', mfilename()]);
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
end

function test_order
%for reference, see Rubik-3-facelet-48.png
%note that moves in this case correspond o anti-clockwise rotations
    if ReplabTestParameters.onlyFastTests
        return
    end
    F = cycle([1:48], [38,17,43,8], [39,20,42,5], [40,22,41,3], [9,11,16,14], [10,13,15,12]);
    R = cycle([1:48], [40,25,48,16], [37,28,45,13], [35,30,43,11], [17,19,24,22], [18,21,23,20]);
    L = cycle([1:48], [33,9,41,32], [36,12,44,29], [38,14,46,27], [1,3,8,6], [2,5,7,4]);
    B = cycle([1:48], [35,1,46,24], [34,4,47,21], [33,6,48,19], [26,29,31,28], [25,27,32,30]);
    D = cycle([1:48], [14,22,30,6], [15,23,31,7], [16,24,32,8], [42,45,47,44], [41,43,48,46]);
    U = cycle([1:48], [27, 19, 11, 3], [26, 18, 10, 2], [25, 17, 9, 1], [33, 35, 40, 38], [34, 37, 39, 36]);
    RC = replab.Permutations(48).subgroup({F, R, L, B, D, U});
    assert(RC.order == vpi('43252003274489856000'));
end
