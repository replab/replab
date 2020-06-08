function elmtStr(elements)
% prints out up to 10 of a cell array of elements and the array length

nelmts = 12;

if length(elements) <= nelmts
    for i = 1:length(elements)
        disp(replab.str.shortStr(elements{i}, []))
    end

else
    for i = 1:nelmts
        disp(replab.str.shortStr(elements{i}, []))
    end
    disp(['.. ', num2str(length(elements) - nelmts), ' elements omitted'])
end

end