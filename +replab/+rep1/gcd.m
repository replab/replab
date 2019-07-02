function g = gcd(vec)
    assert(length(vec) > 0);
    g = vec(1);
    for i = 2:length(vec)
        g = gcd(g, vec(i));
    end
end
