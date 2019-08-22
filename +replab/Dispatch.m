classdef Dispatch < handle
    
    properties (SetAccess = protected)
        names
        handles
        priorities
    end
    
    methods
        
        function self = Dispatch
            self.names = {};
            self.handles = {};
            self.priorities = [];
        end
    
        function register(self, name, handle, priority)
            newPriorities = [self.priorities priority];
            % Descreasing order
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
    
    methods (Static)
        
        function error(message)
            if nargin < 1
                message = '';
            end
            error('replab:Dispatch:tryNext', message);
        end
        
        function assert(condition, message)
            if nargin < 2
                message = '';
            end
            if ~condition
                error('replab:Dispatch:tryNext', message);
            end
        end
        
    end
end
