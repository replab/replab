function D = intAsDoubleMatrix(nRows, nCols, rangeMin, rangeMax)
    rangeMin = double(rangeMin);
    rangeMax = double(rangeMax);
    if rangeMax < rangeMin
        sampleFun = @() replab.domain.inexistent('Cannot sample from empty range');
    else
        sampleFun = @() randi([rangeMin rangeMax], nRows, nCols);
    end
    desc = sprintf('%d x %d integer (double) matrix with coefficients between %d and %d', nRows, nCols, rangeMin, rangeMax);
    D = replab.DomainFun(desc, @isequal, sampleFun);
end
