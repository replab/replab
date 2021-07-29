function res = starCyclotomic(cyc)
% Returns the unique nontrivial Galois conjugate of cyclotomics or zero
%
% Args:
%   cyc (`+replab.cyclotomic`): Scalar cyclotomic number
%
% Returns:
%   `+replab.cyclotomic` or ``[]``: Unique Galois conjugate or ``y[]``
    res = []; % fail
    assert(numel(cyc) == 1);
    javaCyc = cyc.data;
    javaCyc = javaCyc(1);
    n = javaMethod('order', javaCyc);
    if n == 1
        return
    end
    gens = replab.numerical.integer.generatorsPrimeResidues(n);
    gens = [gens{:}];
    cand = [];
    for exp = gens
        img = galois(cyc, exp);
        if img ~= cyc
            if isempty(cand)
                cand = img;
            elseif cand ~= img
                return
            end
        end
    end
    for exp = gens
        img = galois(cand, exp);
        if img ~= cyc && img ~= cand
            return
        end
    end
    res = cand;
end
