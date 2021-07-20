function H = inductiveStep(H, generators, identity)
% Inductive step of Dimino's algorithm
    subgroupSize = H.nElements;
    cosetReps = {identity};
    i = 1;
    while i <= length(cosetReps)
        g = cosetReps{i};
        for j = 1:length(generators)
            s = generators{j};
            newEl = g * s;
            if H.find(newEl) == 0
                cosetReps{1,end+1} = newEl;
                for k = 1:subgroupSize
                    H.insert(H.at(k) * newEl);
                end
            end
        end
        i = i + 1;
    end
end
