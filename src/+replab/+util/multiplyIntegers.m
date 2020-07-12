function o = multiplyIntegers(factors)
    if exist('java.math.BigInteger')
        o = java.math.BigInteger.valueOf(factors(1));
        for i = 2:length(factors)
            o = o.multiply(java.math.BigInteger.valueOf(factors(i)));
        end
        o = vpi(char(o.toString));
    else
        o = factors(1);
        i = 2;
        n = length(factors);
        while i <= length(n) && log2(o) + log2(factors(i)) < 53
            o = o * factors(i);
            i = i + 1;
        end
        o = vpi(o);
        while i <= n
            o = o * factors(i);
            i = i + 1;
        end
    end
end
