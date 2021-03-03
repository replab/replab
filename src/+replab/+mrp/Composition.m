classdef Composition < replab.Morphism

    properties (SetAccess = protected)
        first % (`+replab.Morphism`): First morphism
        second % (`+replab.Morphism`): Second morphism
        imageElementFun % (function_handle or ``[]``): Optional function_handle that implements `.imageElement`
    end


    methods

        function self = Composition(second, first, imageElementFun)
            self.target = second.target;
            self.source = first.source;
            self.first = first;
            self.second = second;
            if ~isempty(first.torusMap) && ~isempty(second.torusMap)
                self.torusMap = second.torusMap * first.torusMap;
            end
            self.imageElementFun = imageElementFun;
        end

    end

    methods % Implementations

        function t = imageElement(self, s)
            if ~isempty(self.imageElementFun)
                f = self.imageElementFun;
                t = f(s);
                return
            end
            t = self.second.imageElement(self.first.imageElement(s));
        end

    end

end
