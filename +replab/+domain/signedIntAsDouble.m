function D = signedIntAsDouble(rangeMin, rangeMax)
    assert(rangeMin >= 0);
    assert(rangeMax >= 0);
    rangeMin = double(rangeMin);
    rangeMax = double(rangeMax);
    if rangeMax < rangeMin
        sampleFun = @() replab.domain.inexistent('Cannot sample from empty range');
    else
        sampleFun = @() randi([rangeMin rangeMax])*(randi([0 1])*2-1);
    end
    desc = sprintf('Signed integers (double) in {-%d,...,-%d,%d,...,%d}', ...
                   rangeMax, rangeMin, rangeMin, rangeMax);
    D = replab.DomainFun(desc, @isequal, sampleFun);
end
