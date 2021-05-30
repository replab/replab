

function g = primitiveRootMod(p)

if ~isprime(p), error('p is not a prime'); end;

primef = unique(factor(p-1));
g = 2; i = 1;
while i <= length(primef)
    if powmod(g, (p-1)/primef(i), p) == 1
        g = g + 1; i = 0;
    end;
    i = i + 1;
end;
