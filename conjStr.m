function conjStr(conjugacyclass)
% prints out up to 12 of elements of up to 3 conjugacy classes

nclasses = 3;
maxCol = 100;

if length(conjugacyclass) <= nclasses
    for i = 1:length(conjugacyclass)
        disp(['Conjugacy class representative ', ... 
            replab.shortStr(conjugacyclass{i}{1}, maxCol)])
        elmtStr(conjugacyclass{i})
    end
else
    for i = 1:nclasses
        disp(['Conjugacy class representative ', ... 
            replab.shortStr(conjugacyclass{i}{1}, maxCol)])
        elmtStr(conjugacyclass{i})
    end
    disp(['.. ', num2str(length(conjugacyclass{i}) - nclasses), ...
        ' conjugacy classes omitted'])
end

end