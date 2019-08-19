function varargout = subsref(self, varargin)
% Extract one part of the SDP matrix
    switch varargin{1}(1).type
        case '()'
            [varargout{1:nargout}] = subsref(self.fullMatrix, varargin{1});
        case '.'
            % This was actually a function call
            [varargout{1:nargout}] = builtin('subsref',self,varargin{1});
        otherwise
            error('Not a valid indexing expression');
    end
end
