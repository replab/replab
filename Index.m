function index = Index(perm)
% Note: perm must be a permutation group element i.e. [1, 3, 2]

eltlen = length(perm);

% concatenate numbers to from unique identifier of permutation group
index = 0;
j = 1;
for i = 0:(eltlen-1)
    index = index + perm(eltlen - i) * 10^i * j;
    j = j * 10^(fix(perm(eltlen - i)/10));
end
end