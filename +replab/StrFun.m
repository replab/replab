classdef StrFun < replab.Str
% Base class that enables the use of lambda functions to print the object 

    properties (SetAccess = protected)
        headerStrFun; % function handle with the conventions of replab.Str.headerStr
        shortStrFun;  % function handle with the conventions of replab.Str.shortStr
        longStrFun;   % function handle with the conventions of replab.Str.longStr
    end
    
    methods
        
        function self = StrFun(headerStrFun, shortStrFun, longStrFun)
        % Constructs a StrFun instance with the given (optional) printing functions
        %
        % Each parameter is optional, can be omitted by providing [] instead,
        % or can be a constant character string; adaptations are made in this
        % constructor
            if nargin < 1
                headerStrFun = [];
            elseif isa(headerStrFun, 'char')
                headerStrFun = @(self) headerStrFun;
            end
            if nargin < 2
                shortStrFun = [];
            elseif isa(shortStrFun, 'char')
                shortStrFun = @(self, mc) shortStrFun;
            end
            if nargin < 3
                longStrFun = [];
            elseif isa(longStrFun, 'char')
                longStrFun = @(self, mr, mc) longStrFun;
            end
            self.headerStrFun = headerStrFun;
            self.shortStrFun = shortStrFun;
            self.longStrFun = longStrFun;
        end
        
        function names = hiddenFields(self)
            names = hiddenFields@replab.Str(self);
            names{1, end+1} = 'headerStrFun';
            names{1, end+1} = 'shortStrFun';
            names{1, end+1} = 'longStrFun';
        end
            
    end
    
end
