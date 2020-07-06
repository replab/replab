function test_suite = HigmanSimsGroupTest()
    disp(['Setting up tests in ', mfilename()]);
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
end

function test_order
    if ReplabTestParameters.onlyFastTests
        return
    end
    S100 = replab.S(100);
    g1 = S100.fromCycles([1, 60], [2, 72], [3, 81], [4, 43], [5, 11], [6, 87], [7, 34], [9, 63], [12, 46], [13, 28], [14, 71], [15, 42], [16, 97], [18, 57], [19, 52], [21, 32], [23, 47], [24, 54], [25, 83], [26, 78], [29, 89], [30, 39], [33, 61], [35, 56], [37, 67], [44, 76], [45, 88], [48, 59], [49, 86], [50, 74], [51, 66], [53, 99], [55, 75], [62, 73], [65, 79], [68, 82], [77, 92], [84, 90], [85, 98], [94, 100]);

    g2 = S100.fromCycles([1, 86, 13, 10, 47], [2, 53, 30, 8, 38], [3, 40, 48, 25, 17], [4, 29, 92, 88, 43], [5, 98, 66, 54, 65], [6, 27, 51, 73, 24], [7, 83, 16, 20, 28], [9, 23, 89, 95, 61], [11, 42, 46, 91, 32], [12, 14, 81, 55, 68], [15, 90, 31, 56, 37], [18, 69, 45, 84, 76], [19, 59, 79, 35, 93], [21, 22, 64, 39, 100], [26, 58, 96, 85, 77], [33, 52, 94, 75, 44], [34, 62, 87, 78, 50], [36, 82, 60, 74, 72], [41, 80, 70, 49, 67], [57, 63, 71, 99, 97]);

    HS = S100.subgroup({g1, g2});
    assert(HS.order == 44352000);
end
