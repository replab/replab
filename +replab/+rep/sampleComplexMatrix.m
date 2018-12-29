function M = sampleComplexMatrix(nRows, nCols)
    M = (randn(nRows, nCols) + randn(nRows, nCols) * 1i)/sqrt(2);
end
