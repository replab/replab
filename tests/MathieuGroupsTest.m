function test_suite = MathieuGroupsTest()
    disp(['Setting up tests in ', mfilename()]);
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
end

function test_M11
% Mathieu group M11
    if ReplabTestParameters.onlyFastTests
        return
    end
    S11 = replab.Permutations(11);
    g1 = S11.fromCycles([2,10], [4,11], [5,7], [8,9]);
    g2 = S11.fromCycles([1,4,3,8], [2,5,6,9]);
    M11 = S11.subgroup({g1, g2});
    assert(M11.order == 7920);
end

function test_M12
% Mathieu group M12
    S12 = replab.Permutations(12);
    g1 = S12.fromCycles([1,4], [3,10], [5,11], [6,12]);
    g2 = S12.fromCycles([1,8,9], [2,3,4], [5,12,11], [6,10, 7]);
    M12 = S12.subgroup({g1, g2});
    assert(M12.order == 95040);
end

function test_M22
% Mathieu group M22
    if ReplabTestParameters.onlyFastTests
        return
    end
    S22 = replab.Permutations(22);
    g1 = S22.fromCycles([1,13], [2,8], [3,16], [4,12], [6,22], [7,17], [9, 10], [11,14]);
    g2 = S22.fromCycles([1,22,3,21], [2,18,4,13], [5,12], [6,11,7,15], [8,14,20,10], [17, 19]);
    M22 = S22.subgroup({g1, g2});
    assert(M22.order == 443520);
end


function test_M23
% Mathieu group M23
    if ReplabTestParameters.onlyFastTests
        return
    end
    S23 = replab.Permutations(23);
    g1 = S23.fromCycles([1,2], [3,4], [7,8], [9,10], [13,14], [15,16], [19,20], [21,22]);
    g2 = S23.fromCycles([1,16,11,3], [2,9,21,12], [4,5,8,23], [6,22,14,18], [13,20], [15,17]);
    M23 = S23.subgroup({g1, g2});
    assert(M23.order == 10200960);
end

function test_M24
    if ReplabTestParameters.onlyFastTests
        return
    end
    S24 = replab.Permutations(24);
    g1 = S24.fromCycles([1,4], [2,7], [3,17], [5,13], [6,9], [8,15], [10, 19], [11, 18], [12, 21], [14, 16], [20, 24], [22, 23]);
    g2 = S24.fromCycles([1,4,6], [2,21,14], [3, 9, 15], [5, 18, 10], [13, 17, 16], [19, 24, 23]);
    M24 = S24.subgroup({g1, g2});
    assert(M24.order == 244823040);
end
