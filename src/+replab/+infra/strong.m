function res = strong(str)
% Returns the text surrounded by strong HTML tags if the console output supports that
    if replab.globals.consoleUseHTML
        res = ['<strong>' str '</strong>'];
    else
        res = str;
    end
end
