classdef Context < replab.Str

    properties (SetAccess = protected)
        id % (charstring): Unique ID of the context
        closed % (logical): Whether the context has been closed
    end

    properties (Access = protected)
        instances_ % (cell(1,*) of `+replab.Equivariant`): Equivariant instances with cached samples
    end

    methods (Access = protected)

        function self = Context(id)
            self.id = id;
            self.closed = false;
        end

    end

    methods

        function register(self, e)
        % Registers an equivariant instance to clean the cache of
        %
        % Args:
        %   e (`replab.Equivariant`): Equivariant instance to register
            assert(~self.closed);
            self.instances_{1,end+1} = e;
        end

        function close(self)
        % Cleans the caches of this context
            assert(~self.closed);
            for i = 1:length(self.instances_)
                self.instances_{i}.clearCache(self);
            end
            self.closed = true;
        end

    end

    methods (Static)

        function id = indexToId(i)
        % Converts a 1-based index below 2^32 to a 7 character alphabetical string
            assert(1 <= i && i < 2^32);
            i = i - 1;
            id = repmat(' ', 1, 7);
            for p = 7:-1:1
                m = mod(i, 26);
                i = (i - m)/26;
                id(p) = char('a' + m);
            end
        end

        function c = make
            persistent index
            if isempty(index)
                index = 1;
            end
            id = replab.Context.indexToId(index);
            c = replab.Context(id);
            index = index + 1;
        end

    end

end
