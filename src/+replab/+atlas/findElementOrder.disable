function el = findElementOfOrder(group, prop)
% DOES NOT WORK FOR NOW
% Finds an element of a given order satisfying the given property
    assert(isa(elementOrder, 'double'));
    n = chain.n;
    L = chain.length;
    stackG = zeros(n, L+1); % current group element
    stackG(:,1) = 1:n;
    stackS = cell(1, L+1); % subgroup chain
    stackS{1} = chain;
    stackB = cell(1, L+1);
    stackB{1} = sort(stackS{1}.Delta{1});
    stackO = zeros(2, L+1); % first line is where to go, second line is lcm so far
    stackO(:, 1) = [0; 1];
    cycleLengths = arrayfun(@(i) mod(elementOrder, i) == 0, 1:n+1);
    l = 1;
    while l >= 1
        g = stackG(:,l)';
        while 1
            B = stackB{l};
            if isempty(B)
                l = l - 1;
                break
            end
            b = B(1);
            stackB{l} = B(2:end);
            if stackO(1,l) == 0
                stackO(:,l+1) = [stackS{l}.B(1); stackO(2,l)];
            elseif stackO(1,l) == b
                cycleLength = sum(stackO(1:l) == b) + 1;
                if ~cycleLengths(cycleLength)
                    continue
                end
                stackO(:,l+1) = [0; lcm(stackO(2,l), cycleLength)];
            else
                stackO(:,l+1) = stackO(:,l);
            end
            gnext = g(stackS{l}.u(1, b));
            if l == L
                if stackO(2,l+1) == elementOrder && prop(gnext)
                    gnext
                end
            else
                stackS{l+1} = stackS{l}.stabilizer(b);
                stackG(:,l+1) = gnext;
                stackB{l+1} = sort(stackS{l+1}.Delta{1});
                l = l + 1;
                break
            end
        end
    end
    el = [];
end

function [ok outdata] = finaltest(l, elementOrder, g, indata)
    outdata = indata(1);
    if length(indata) == 1
        assert(g(l) == l);
        ok = indata(1) == elementOrder;
    else
        assert(g(l) == indata(2));
        cycleLength = indata(3) + 1;
        ok = lcm(indata(1), cycleLength) == elementOrder;
    end
end

function [ok outdata] = test(l, cycleLengths, g, indata)
% indata is either [lcmOfCycleLengths cycleStart count] or [lcm]
    currentLcm = indata(1);
    if length(indata) == 1
        if g(l) == l
            ok = true;
            outdata = currentLcm;
        else
            ok = true;
            outdata = [currentLcm l 1];
        end
    else
        cycleStart = indata(2);
        if g(l) == cycleStart
            cycleLength = indata(3) + 1;
            if cycleLengths(cycleLength)
                ok = true;
                outdata = [lcm(currentLcm, cycleLength)];
            else
                ok = false;
                outdata = [];
            end
        else
            ok = true;
            outdata = indata + [0 0 1]; % increment the cycle length
        end
    end
end
