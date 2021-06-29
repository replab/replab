function ub = condUpperBound(A, Ainv)
    if isa(A, 'replab.cyclotomic')
        A = double(A);
    end
    ub = replab.numerical.norm2UpperBound(A) * replab.numerical.norm2UpperBound(Ainv);
end
