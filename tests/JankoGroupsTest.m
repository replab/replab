function test_suite = JankoGroupsTest()
    disp(['Setting up tests in ', mfilename()]);
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
    %These groups have irreps with multiplicity 1, but high dimensions
end

function test_J1
% Janko group J1
%
% See https://en.wikipedia.org/wiki/Janko_group_J1
    if ReplabTestParameters.onlyFastTests
        return
    end
    S266 = replab.Permutations(266);
    g1 = S266.fromCycles([1, 262], [2, 107], [3, 21], [4, 213], [5, 191], [6, 22], [7, 133], [8, 234], [9, 232], [10, 151],  [11, 139], [12, 176], [13, 202], [14, 253], [15, 222], [17, 195], [18, 206], [19, 68], [20, 55], [23, 179], [24, 217], [25, 216], [26, 256], [27, 87], [28, 70], [29, 131], [30, 44], [31, 105], [32, 170], [33, 77], [34, 104], [35, 198], [36, 137], [37, 243], [38, 56], [39, 124], [40, 223], [41, 134], [43, 174], [46, 51], [47, 128], [48, 94], [49, 250], [50, 264], [52, 183], [53, 231], [54, 115], [57, 85], [58, 233], [59, 261], [60, 95], [61, 235], [62, 177], [63, 249], [64, 91], [65, 247], [66, 155], [69, 219], [71, 237], [72, 211], [73, 84], [74, 192], [75, 130], [76, 251], [79, 260], [80, 112], [81, 193], [82, 156], [83, 242], [86, 238], [88, 143], [89, 168], [90, 148], [92, 119], [93, 212], [96, 150], [97, 199], [98, 140], [99, 189], [100, 180], [101, 147], [102, 111], [103, 159], [106, 162], [108, 194], [109, 166], [110, 200], [113, 120], [114, 141], [116, 182], [117, 181], [118, 225], [121, 254], [122, 125], [123, 146], [126, 208], [127, 221], [129, 210], [132, 255], [136, 175], [138, 207], [142, 240], [144, 172], [145, 185], [149, 224], [152, 169], [153, 241], [154, 190], [157, 214], [158, 161], [160, 236], [163, 239], [164, 229], [165, 230], [167, 188], [171, 258], [173, 186], [178, 245], [184, 205], [187, 228], [197, 203], [201, 252], [209, 248], [215, 259], [218, 246], [220, 227], [257, 263], [265, 266]);

    g2 = S266.fromCycles([1, 146, 21], [2, 132, 82], [4, 156, 166], [5, 242, 253], [6, 107, 28], [7, 125, 76], [8, 245, 130], [9, 174, 42], [10, 241, 244], [11, 264, 63], [12, 248, 234], [13, 36, 44], [14, 116, 128], [15, 47, 25], [16, 178, 112], [17, 170, 110], [18, 197, 74], [19, 233, 180], [20, 121, 96], [22, 228, 155], [23, 48, 173], [24, 201, 187], [26, 136, 190], [27, 212, 94], [29, 175, 52], [30, 77, 32], [31, 237, 34], [33, 226, 90], [35, 129, 54], [37, 161, 114], [38, 232, 87], [39, 219, 192], [40, 78, 159], [41, 139, 71], [43, 211, 251], [45, 222, 240], [46, 97, 135], [49, 70, 131], [50, 153, 200], [51, 186, 209], [53, 203, 216], [55, 169, 64], [56, 140, 230], [57, 260, 118], [58, 91, 243], [59, 199, 227], [60, 108, 164], [61, 208, 101], [62, 206, 106], [65, 103, 66], [67, 95, 205], [68, 73, 225], [69, 151, 113], [72, 221, 152], [75, 143, 202], [79, 217, 254], [80, 93, 122], [81, 181, 252], [83, 258, 126], [84, 163, 177], [85, 154, 213], [86, 182, 196], [88, 133, 215], [89, 117, 247], [92, 191, 160], [99, 229, 263], [100, 138, 188], [102, 194, 157], [105, 149, 184], [109, 123, 193], [111, 137, 183], [115, 238, 235], [119, 167, 147], [120, 134, 189], [124, 185, 265], [127, 218, 261], [141, 231, 210], [142, 239, 236], [144, 224, 249], [145, 158, 220], [148, 214, 172], [150, 250, 259], [162, 257, 256], [165, 179, 246], [176, 195, 266], [198, 204, 207], [223, 262, 255]);

    J1 = S266.subgroup({g1, g2});
    assert(J1.order == 175560);
end

function test_J2
% Janko group J2
%
% See https://en.wikipedia.org/wiki/Janko_group_J2
    if ReplabTestParameters.onlyFastTests
        return
    end
    S100 = replab.Permutations(100);
    g1 = S100.fromCycles([1, 84], [2, 20], [3, 48], [4, 56], [5, 82], [6, 67], [7, 55], [8, 41], [9, 35], [10, 40], [11, 78], [12, 100], [13, 49], [14, 37], [15, 94], [16, 76], [17, 19], [18, 44], [21, 34], [22, 85], [23, 92], [24, 57], [25, 75], [26, 28], [27, 64], [29, 90], [30, 97], [31, 38], [32, 68], [33, 69], [36, 53], [39, 61], [42, 73], [43, 91], [45, 86], [46, 81], [47, 89], [50, 93], [51, 96], [52, 72], [54, 74], [58, 99], [59, 95], [60, 63], [62, 83], [65, 70], [66, 88], [71, 87], [77, 98], [79, 80]);  
    g2 = S100.fromCycles([1, 80, 22], [2, 9, 11], [3, 53, 87], [4, 23, 78], [5, 51, 18], [6, 37, 24], [8, 27, 60], [10, 62, 47], [12, 65, 31], [13, 64, 19], [14, 61, 52], [15, 98, 25], [16, 73, 32], [17, 39, 33], [20, 97, 58], [21, 96, 67], [26, 93, 99], [28, 57, 35], [29, 71, 55], [30, 69, 45], [34, 86, 82], [38, 59, 94], [40, 43, 91], [42, 68, 44], [46, 85, 89], [48, 76, 90], [49, 92, 77], [50, 66, 88], [54, 95, 56], [63, 74, 72], [70, 81, 75], [79, 100, 83]);
    J2 = S100.subgroup({g1, g2});
    assert(J2.order == 604800);
end
