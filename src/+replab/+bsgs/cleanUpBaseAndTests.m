function [group1, tests1, startData1] = cleanUpBaseAndTests(group, tests, startData)
% Removes redundant base points from ``group`` and fixes the rest of the data if necessary
    for i = length(tests)+1:group.length
        % handle tests: if tests do not cover base fully, we complete by trivial tests
        tests{1,i} = @(g, indata) deal(true, []);
    end
    if any(group.orbitSizes == 1)
        base = group.B;
        group1 = group.mutableCopy;
        group1.removeRedundantBasePoints;
        group1.makeImmutable;
        if all(group.orbitSizes == 1)
            tests1 = {};
            startData1 = [];
            return
        end
        base1 = group1.base;
        start = find(base == base1(1));
        startData1 = startData;
        identity = 1:group.n;
        for j = 1:start-1
            [ok, startData1] = tests{j}(identity, startData1);
            assert(ok == 1);
        end
        current = start;
        for i = 1:length(base1)
            if i < length(base1)
                next = find(base == base1(i+1));
            else
                next = length(base) + 1;
            end
            if current + 1 == next
                tests1{i} = tests{current};
            else
                tests1{i} = @(g, inData) replab.bsgs.concatTests(tests(current:next-1), g, inData);
            end
        end
    else
        group1 = group;
        tests1 = tests;
        startData1 = startData;
    end
end
