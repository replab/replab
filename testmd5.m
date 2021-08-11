m1 = replab.util.md5('Test');
m2 = replab.util.md5_plain('Test');
assert(isequal(m1, m2))