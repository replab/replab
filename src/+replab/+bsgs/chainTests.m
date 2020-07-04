function [ok, data] = chainTests(tests, g, data)
    for i = 1:length(tests)
        [ok, data] = tests{i}(g, data);
        if ~ok
            return
        end
    end
end
