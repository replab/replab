function c = javaArrayToCell(ja)
    if replab.compat.isOctave
        c = cell(1, length(ja));
        for i = 1:length(ja)
            c{i} = ja(i);
        end
    else
        c = cell(ja);
    end
end
