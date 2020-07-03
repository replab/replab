function ordering = baseOrdering(degree, base)
    base_len = length(base);
    ordering = zeros(1, degree);
    ordering(base) = 1:base_len;
    rest = setdiff(1:degree, base);
    ordering(rest) = base_len + (1:length(rest));
end
