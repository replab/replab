classdef Obj < replab.Str
% Base class that provides sane pretty printing and instance equality tests

    properties (Access = protected)
        cache_ % (struct): Contains the computed properties of this object
        id_ % (integer): Unique object ID
    end

    methods % Laws

        function check(self)
        % Checks the consistency of this object
        %
        % The default implementation checks the declared laws.
            self.laws.check;
        end

        function l = laws(self)
        % Returns the laws that this object obeys
        %
        % Returns:
        %   `+replab.Laws`: The `.Laws` instance that is relevant for this object.
            l = replab.Laws; % Empty instance by default
        end

    end

    methods % Property cache

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
        % - ``error``: Raises an error
        %
        % Note that function handles cannot be stored in the cache (if that is necessary, they can be wrapped in
        % a scalar cell array); function handles are used to provide lazy evaluation of properties, and will be
        % called once when the property is requested.
        %
        % Args:
        %   name (charstring): Name of the property
        %   value: Value of the property, or a function handle able to compute the value
        %   handleExisting ({'overwrite', 'ignore', 'isequal', '=='}, optional): What to do if the value is already known, default ``'ignore'``
            if isempty(self.cache_)
                self.cache_ = struct;
            end
            if self.inCache(name)
                if nargin <4
                    handleExisting = 'ignore';
                end
                switch handleExisting
                  case 'error'
                    error('Value for %s already present in cache', name);
                  case 'overwrite'
                    self.cache_.(name) = value;
                  case 'ignore'
                  case '=='
                    assert(self.cachedOrEmpty(name) == value);
                  case 'isequal'
                    assert(isequal(self.cachedOrEmpty(name), value));
                  otherwise
                    error('Invalid argument handleExisting')
                end
            else
                self.cache_.(name) = value;
            end
        end

        function res = cachedOrDefault(self, name, defaultValue)
        % Returns the cached property if it exists, or the provided default value if it is unknown yet
        %
        % See `.cached` and `.cachedOrEmpty`.
        %
        % Args:
        %   name (charstring): Name of the property
        %   defaultValue: Value returned in case the property is unknown
        %
        % Returns
        %   The property value if known, otherwise ``defaultValue``
            if self.inCache(name)
                res = self.cache_.(name);
                if isa(res, 'function_handle')
                    res = res();
                    self.cache_.(name) = res;
                end
            else
                res = defaultValue;
            end
        end

        function res = cachedOrEmpty(self, name)
        % Returns the cached property if it exists, or ``[]`` if it is unknown yet
        %
        % See `.cached` and `.cachedOrDefault`.
        %
        % Args:
        %   name (charstring): Name of the property
        %
        % Returns:
        %   The property value if known, otherwise ``[]``
            res = self.cachedOrDefault(name, []);
        end

        function res = cached(self, name, fun)
        % Returns the cached property if it exists, computing it if necessary
        %
        % If the property is not known, it will call the given function handle and cache the result.
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
                if isa(res, 'function_handle')
                    res = res();
                    self.cache_.(name) = res;
                end
            else
                res = fun();
                self.cache_.(name) = res;
            end
        end

    end

    methods % Unique ID

        function i = id(self)
        % Returns the unique ID of this object (deprecated)
            if isempty(self.id_)
                self.id_ = replab.globals.nextUniqueId;
            end
            i = self.id_;
        end

        function res = ne(self, rhs)
        % Non-equality test (deprecated)
        %
        % Workaround bug of == not implemented for handles in Octave
        %
        % Args:
        %   self (object): first object
        %   rhs (object): second object to compare to
        %
        % Returns:
        %   logical: true iff self ~= rhs
            res = ~self.eq(rhs);
        end

        function res = eq(self, rhs)
        % Equality test (deprecated)
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
