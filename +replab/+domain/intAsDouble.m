function D = intAsDouble(rangeMin, rangeMax)
    rangeMin = double(rangeMin);
    rangeMax = double(rangeMax);
    if rangeMax < rangeMin
        sampleFun = @() replab.domain.inexistent('Cannot sample from empty range');
    else
        sampleFun = @() randi([rangeMin rangeMax]);
    end
    desc = sprintf('Integers (double) between %d and %d', rangeMin, rangeMax);
    D = replab.Domain.lambda(desc, @isequal, sampleFun);
end
