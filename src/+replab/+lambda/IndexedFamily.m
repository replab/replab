classdef IndexedFamily < replab.IndexedFamily
% An implementation of an indexed family defined by image functions

    properties (SetAccess = protected)
        atFun % (function_handle): Handle that implements `+replab.IndexedFamily.at`
        findFun % (function_handle): Handle that implements `+replab.IndexedFamily.find`
    end

    methods

        function self = IndexedFamily(size, atFun, findFun)
            if isa(size, 'vpi')
                self.size = size;
            else
                self.size = vpi(size);
            end
            self.atFun = atFun;
            self.findFun = findFun;
        end

        function names = hiddenFields(self)
            names = replab.str.uniqueNames( ...
                hiddenFields@replab.IndexedFamily(self), ...
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
