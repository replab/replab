function ub = norm2UpperBound(A)
    if isa(A, 'replab.cyclotomic')
        A = double(A);
    end
    ub1 = norm(A, 'fro'); % ||A||_2 <= ||A||_fro
    ub2 = sqrt(norm(A, 1)*norm(A, inf)); % ||A||_2 <= sqrt(||A||_1 ||A||_inf)
    ub = min(ub1, ub2);
    if isa(ub, 'intval')
        ub = sup(ub);
    end
end
