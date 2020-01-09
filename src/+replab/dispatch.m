function varargout = dispatch(cmd, name, varargin)
% Supports multiple and flexible dispatch through a registry of functions
%
% For the 'register' command, the next four arguments are 'name', ``description``,
% ``priority`` and ``handle`` in that order. The argument ``name`` corresponds to the
% name of the function implementing multiple dispatch; ``description`` is a short
% description of the particular implementation of that funciton, while
% ``priority`` is an integer with bigger numbers take higher priority, and
% finally ``handle`` is a function handle implementing the function.
%
% For the ``get`` command, no additional arguments are required and the whole
% dispatch registry is returned.
%
% For the 'call' command, the ``name`` argument is the function being dispatched,
% followed by that function arguments.
%
% Args:
%   cmd ({'register', 'get', 'call'}): Action to take
%   name (char, optional): Function name to act on, corresponds to the fully qualified function identifier
%                          such as ``replab.makeEquivariant``
%   varargin (optional): Additional arguments

    persistent registry;
    if isempty(registry)
        registry = struct;
        replab.dispatchDefaults;
    end
    switch cmd
      case 'register'
        assert(isa(name, 'char'));
        ident = strrep(name, '.', '_');
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
      case 'get'
        varargout = {registry};
      case 'call'
        ident = strrep(name, '.', '_');
        s = registry.(ident);
        for i = 1:length(s)
            h = s(i).handle;
            try
                [varargout{1:nargout}] = h(varargin{:});
                return
            catch ME
                if ~strcmp(ME.identifier, 'replab:dispatch:tryNext')
                    rethrow(ME)
                end
            end
        end
        error('No registered implementation worked.')
    end
end
