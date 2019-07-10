function l = veclcm(vec)
    assert(length(vec) > 0);
    l = vec(1);
    for i = 2:length(vec)
        l = lcm(l, vec(i));
    end
end
