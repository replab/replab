classdef Obj < replab.Str
% Base class that provides sane pretty printing and instance equality tests

    methods (Static)

        function id = newId
            persistent lastId
            if isempty(lastId)
                lastId = uint64(0);
            end
            id = lastId + uint64(1);
            lastId = id;
        end

    end

    properties (Access = protected)
        cache_ % (struct): Contains the computed properties of this object
        id_ % (uint64): Unique object ID
    end

    methods

        function l = inCache(self, name)
        % Returns whether the value of the given property has already been computed
        %
        % Args:
        %   name (charstring): Name of the property
        %
        % Returns:
        %   logical: Whether the property value is in the cache
            if isempty(self.cache_)
                self.cache_ = struct;
            end
            l = isfield(self.cache_, name);
        end

        function cache(self, name, value, handleExisting)
        % Sets the value of the designated property in the cache
        %
        % If the property value is already in the cache, the behavior depends on the argument ``handleExisting``:
        %
        % - ``overwrite``: The existing value is overwritten
        % - ``ignore``: The existing value is left unchanged
        % - ``isequal``: The existing value and the given value are compared for equality using ``isequal``
        % - ``=``: The existing value and the given value are compared for equality using ``==``
        %
        % Args:
        %   name (charstring): Name of the property
        %   value: Value of the property
        %   handleExisting ({'overwrite', 'ignore', 'isequal', '=='}, optional): What to do if the value is already known, default ``'ignore'``
            if isempty(self.cache_)
                self.cache_ = struct;
            end
            if self.inCache(name)
                if nargin <4
                    handleExisting = 'ignore';
                end
                switch handleExisting
                  case 'overwrite'
                    self.cache_.(name) = value;
                  case 'ignore'
                  case '=='
                    assert(self.cache_.(name) == value);
                  case 'isequal'
                    assert(isequal(self.cache_.(name), value));
                  otherwise
                    error('Invalid argument handleExisting')
                end
            else
                self.cache_.(name) = value;
            end
        end

        function res = cachedOrEmpty(self, name)
        % Returns the cached property if it exists, or ``[]`` if it is unknown yet
        %
        % See `.cached`.
        % Args:
        %   name (charstring): Name of the property
        %
        % Returns:
        %   The property value if known, otherwise ``[]``
            if isempty(self.cache_)
                self.cache_ = struct;
            end
            if isfield(self.cache_, name)
                res = self.cache_.(name);
            else
                res = [];
            end
        end

        function res = cached(self, name, fun)
        % Returns the cached property if it exists, computing it if necessary
        %
        % If the property is not known, it will call the given function handle and
        % cache the result.
        %
        % Args:
        %   name (charstring): Name of the property
        %   fun (function_handle): Function that computes the property
        %
        % Returns:
        %   The property value
            if isempty(self.cache_)
                self.cache_ = struct;
            end
            if isfield(self.cache_, name)
                res = self.cache_.(name);
            else
                res = fun();
                self.cache_.(name) = res;
            end
        end

        function i = id(self)
            if isempty(self.id_)
                self.id_ = replab.Obj.newId;
            end
            i = self.id_;
        end

        function res = ne(self, rhs)
        % Non-equality test
        %
        % Workaround bug of == not implemented for handles in Octave
        %
        % Args:
        %   self (object): first object
        %   rhs (object): second object to compare to
        %
        % Returns:
        %   logical: true iff self ~= rhs
            res = ~self.rq(rhs);
        end

        function res = eq(self, rhs)
        % Equality test
        %
        % Workaround bug of == not implemented for handles in Octave
        %
        % Args:
        %   self (object): first object
        %   rhs (object): second object to compare to
        %
        % Returns:
        %   logical: true iff self == rhs
            if ~isa(rhs, 'replab.Obj')
                res = false;
                return
            end
            l = arrayfun(@(x) x.id, self);
            r = arrayfun(@(x) x.id, rhs);
            res = (l == r);
        end

    end

end
