classdef Dispatch < handle
% A class handling functions/methods with multiple implementations
    
    properties (SetAccess = protected)
        names % row cell array of char: Function names
        handles % row cell array of function_handle: Function handles
        priorities % row integer vector: Implementation priority
    end
    
    methods
        
        function self = Dispatch
            self.names = {};
            self.handles = {};
            self.priorities = [];
        end
    
        function register(self, name, handle, priority)
            newPriorities = [self.priorities priority];
            % Decreasing order
            p = replab.Permutations.sorting(newPriorities, @(x, y) x < y);
            newPriorities = newPriorities(p);
            newNames = {self.names{:} name};
            newNames = newNames(p);
            newHandles = {self.handles{:} handle};
            newHandles = newHandles(p);
            self.names = newNames;
            self.handles = newHandles;
            self.priorities = newPriorities;
        end
        
        function varargout = call(self, varargin)
            for i = 1:length(self.handles)
                h = self.handles{i};
                try
                    
                    [varargout{1:nargout}] = h(varargin{:});
                    return
                catch ME
                    if ~strcmp(ME.identifier, 'replab:Dispatch:tryNext')
                        rethrow(ME)
                    end
                end
            end
            error('No registered implementation worked.');
        end
        
    end
    
end
