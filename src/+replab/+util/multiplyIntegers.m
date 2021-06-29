function o = multiplyIntegers(factors)
    if isempty(factors)
        o = vpi(1);
        return
    end
    o = factors(1);
    i = 2;
    n = length(factors);
    while i <= n && log2(o) + log2(factors(i)) < 53
        o = o * factors(i);
        i = i + 1;
    end
    if i > n
        o = vpi(o);
        return
    end
    if exist('java.math.BigInteger')
        o = java.math.BigInteger.valueOf(o);
        while i <= n
            o = o.multiply(java.math.BigInteger.valueOf(factors(i)));
            i = i + 1;
        end
        o = vpi(char(o.toString));
    else
        o = vpi(o);
        while i <= n
            o = o * factors(i);
            i = i + 1;
        end
    end
end
