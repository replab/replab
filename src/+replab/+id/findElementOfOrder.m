function el = findElementOfOrder(G, elementOrder, prop)
% Finds an element of a given order satisfying the given property
    assert(isa(elementOrder, 'double'));
    chain = G.niceGroup.chain.mutableCopy;
    n = chain.n;
    chain.baseChange(1:chain.n);
    chain.makeImmutable;
    cycleLengths = arrayfun(@(i) mod(elementOrder, i) == 0, 1:n+1);
    tests = cell(1, n);
    for i = 1:n
        tests{i} = @(x, data) test(i, cycleLengths, x, data);
    end
    if nargin < 3
        fullProp = @(g) replab.Permutation.order(g) == elementOrder;
    else
        fullProp = @(g) fullProperty(g, elementOrder, prop);
    end
    el = replab.bsgs.backtrackSearch(chain, fullProp, tests, []);
end

function ok = fullProperty(g, elementOrder, prop)
    ok = replab.Permutation.order(g) == elementOrder && prop(g);
end

function [ok outdata] = test(l, cycleLengths, g, indata)
% indata is either [cycleStart count] or []
    if isempty(indata)
        if g(l) == l
            ok = true;
            outdata = [];
        else
            ok = true;
            outdata = [l 1];
        end
    else
        if g(l) == indata(1)
            cycleLength = indata(2) + 1;
            if cycleLengths(cycleLength)
                ok = true;
                outdata = [];
            else
                ok = false;
                outdata = [];
            end
        else
            ok = true;
            outdata = indata;
            outdata = indata + [0 1]; % increment the cycle length
        end
    end
end
