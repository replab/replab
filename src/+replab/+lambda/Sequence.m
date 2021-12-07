classdef Sequence < replab.Sequence
% An implementation of a sequence defined by image functions

    properties (SetAccess = protected)
        atFun % (function_handle): Handle that implements `+replab.Sequence.at`
        findFun % (function_handle): Handle that implements `+replab.Sequence.find`
    end

    methods

        function self = Sequence(nElements, atFun, findFun)
            self@replab.Sequence(nElements);
            self.atFun = atFun;
            self.findFun = findFun;
        end

        function names = hiddenFields(self)
            names = replab.str.uniqueNames( ...
                hiddenFields@replab.Sequence(self), ...
                {'atFun' 'sampleFun'} ...
                );
        end

        function obj = at(self, ind)
            obj = self.atFun(vpi(ind));
        end

        function ind = find(self, obj)
            ind = self.findFun(obj);
        end

    end

end
