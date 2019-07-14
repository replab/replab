classdef EnumeratorFun < replab.Enumerator
    properties (SetAccess = protected)
        atFun; % Handle that implements Enumerator.at
        findFun; % Handle that implements Enumerator.find
    end
    methods
        function self = EnumeratorFun(size, atFun, findFun)
            self@replab.Enumerator(size);
            self.atFun = atFun;
            self.findFun = findFun;            
        end
        function names = hiddenFields(self)
            names = hiddenFields@replab.Enumerator(self);
            names{end+1, 1} = 'atFun';
            names{end+1, 1} = 'sampleFun';
            names = unique(names);
        end
    end
end
