function test_suite = HashIntMatrixTest()
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
end

function test_hashIntMatrix()
    A1 = randi(30, 5, 30);
    A2 = randi(30, 5, 20);
    H = replab.prv.HashIntMatrix(A1);
    H1 = H.append(A2);
    H2 = replab.prv.HashIntMatrix([A1 A2]);
    H3 = H1.lexColSorted;
    i = randi(30);
    assertEqual(H1.find(A1(:, i)), i);
    assertEqual(H2.find(A1(:, i)), i);
    assertGreaterThan(H3.find(A1(:, i)), 0);
    i = randn(30);
end
