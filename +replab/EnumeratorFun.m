classdef EnumeratorFun < replab.Enumerator
    properties (SetAccess = protected)
        atFun; % Handle that implements Enumerator.at
        findFun; % Handle that implements Enumerator.find
    end
    methods
        function self = EnumeratorFun(size, atFun, findFun)
            self = self@replab.Enumerator(size);
            self.atFun = atFun;
            self.findFun = findFun;            
        end
    end
end
