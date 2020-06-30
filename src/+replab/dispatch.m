function varargout = dispatch(cmd, name, varargin)
% Supports multiple and flexible dispatch through a registry of functions
%
% This function has an internal registry of implementations which can be enriched,
% probed, used, depending on the passed command ``cmd`` argument.
%
% - For the ``exists`` command: the next argument is ``name``, and the function
%   returns whether the named function has been registered.
%
% - For the ``register`` command:
%   The next four arguments are ``name``, ``description``, ``priority`` and ``handle``,
%   in that order. The argument ``name`` corresponds to the name of the implemented function
%   The argument ``description`` is a short description of the particular implementation of
%   that function, while ``priority`` is an integer with bigger numbers take higher priority, and
%   finally ``handle`` is a function handle implementing the function.
%   The implemented function can return an instance of `+replab.DispatchNext` on call as the
%   first output argument,
%
% - For the ``get`` command: no additional arguments are required, as the whole
%   dispatch registry is returned.
%
% - For the ``call`` command, the ``name`` argument is the function being dispatched,
%   followed by that function arguments.
%
% Args:
%   cmd ({'exists', 'register', 'get', 'call'}): Action to take
%   name (charstring, optional): Function name to act on.
%                                Corresponds to the fully qualified function identifier
%                                such as ``replab.makeEquivariant``
%   varargin (optional): Additional arguments

    persistent registry
    persistent level

    if isempty(level) || ~replab.globals.verboseDispatch
        level = 0;
    end

    if isempty(registry)
        registry = struct;
    end

    if isequal(cmd, 'get')
        varargout = {registry};
        return
    end

    assert(isa(name, 'char'));
    ident = strrep(name, '.', '_');

    switch cmd
      case 'exists'
        varargout = {isfield(registry, ident)};
      case 'register'
        description = varargin{1};
        assert(isa(description, 'char'));
        priority = varargin{2};
        assert(isa(priority, 'double'));
        handle = varargin{3};
        assert(isa(handle, 'function_handle'));
        if ~isfield(registry, ident)
            s = struct('description', {description}, 'priority', {priority}, 'handle', {handle});
            registry.(ident) = s;
        else
            s = registry.(ident);
            ind = length(s) + 1;
            s(ind).description = description;
            s(ind).priority = priority;
            s(ind).handle = handle;
            % sort in decreasing order
            p = replab.Permutations.sorting([s.priority], @(x, y) x < y);
            s = s(p);
            registry.(ident) = s;
        end
      case 'call'
        if replab.globals.verboseDispatch
            fprintf('%sDispatching %s\n', repmat(' ', 1, level), name);
            level = level + 1;
        end
        ident = strrep(name, '.', '_');
        s = registry.(ident);
        n = min(1, nargout);
        res = cell(1, n);
        for i = 1:length(s)
            if replab.globals.verboseDispatch
                fprintf('%sTrying %s\n', repmat(' ', 1, level), s(i).description);
            end
            h = s(i).handle;
            [res{1:n}] = h(varargin{:});
            if ~isa(res{1}, 'replab.DispatchNext')
                if replab.globals.verboseDispatch
                    fprintf('%s Success\n', repmat(' ', 1, level));
                    level = level - 1;
                end
                varargout = res(1:nargout);
                return
            end
            if replab.globals.verboseDispatch
                if isempty(res{1}.message)
                    fprintf('%s Failed\n', repmat(' ', 1, level));
                else
                    fprintf('%s Failed: %s\n', repmat(' ', 1, level), res{1}.message);
                end
            end
        end
        if replab.globals.verboseDispatch
            level = level - 1;
        end
        error('No registered implementation worked.')
    end
end
