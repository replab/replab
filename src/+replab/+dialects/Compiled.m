classdef Compiled < handle
% Describes compiled code
%
% This class provides an interface to code written in another language than 
% Matlab/Octave that needs to be compiled.
%
% To make the function coded in 'aFunction_mex.cpp' available through this 
% interface, a matlab funciton called 'aFunction.m' should be available in 
% the same directory as 'aFunction_mex.cpp'. This file should contain a 
% function, which instantiates the class `+replab.+dialect.Compiled` with 
% the extension '.cpp'.

    properties (SetAccess = protected)
        path % (charstring) directory of the source code
        fileName % (charstring) name of the source file
        extension % (charstring) extension of the source code
        package % (charstring) name of the package containing the function
        triedBefore % (boolean) whether the compiled binary was called before in the current session
        isWorking % (boolean) whether the compiled binary works
    end

    methods

        function self = Compiled(extension)
        % Constructor
        %
        % Constructs an interfact to a code written in another language
        % than Matlab/Octave.
        %
        % Args:
        %   extension (charstring) : extension of the source code file
            [ST, I] = dbstack('-completenames');
            
            self.path = fileparts(ST(2).file);
            self.fileName = ST(2).name;
            self.extension = extension;
            
            % we extract the package that contains the source code
            pathParts = regexp(self.path, '/', 'split');
            startsWithPlus = cellfun(@(x) ((length(x) > 0) && isequal(x(1),'+')), pathParts);
            self.package = '';
            for i = find(startsWithPlus,1) : length(startsWithPlus)
                self.package = [self.package, pathParts{i}(2:end), '.'];
            end
            if length(self.package) >= 1
                self.package = self.package(1:end-1);
            end

            self.triedBefore = false;
            self.isWorking = false;
        end
        
        function results = call(self, varargin)
        % Calls the compiled function
        %
        % This function checks if the corresponding function is
        % operational, trying to compile it if needed. It then calls the
        % function with the provided arguments and returns the (single) 
        % output of that function. If the call to the function fails, the
        % returned value is `+replab.DisplatchNext`.
        %
        % Args:
        %   varargin : all the arguments that need to be passed to the
        %     compiled function
        
            % We try to use the compiled code
            if (~self.triedBefore) || self.isWorking
                firstPartWorks = true;
                if ~self.triedBefore
                    % We try to use the compiled option
                    % Make sure we are in the current path
                    initialPath = pwd;
                    try
                        cd(self.path);

                        needToCompile = true;
                        if (exist([self.fileName, '_mex.', mexext], 'file') == 2) || (exist([self.fileName, '_mex.', mexext], 'file') == 3)
                            % If the program was already compiled, we check that the
                            % compilation is up to date
                            if exist([self.fileName, '_timestamp_', mexext, '.mat'], 'file')
                                savedData = load([self.fileName, '_timestamp_', mexext, '.mat']);
                                fileProperties = dir([self.fileName, '_mex.', self.extension]);
                                if isequal(savedData.timestamp, fileProperties.date)
                                    needToCompile = false;
                                end
                            end
                        end

                        if needToCompile
                            % If the program was never compiled, or the source was
                            % modified, we try to compile it
                            mex('-largeArrayDims', [self.fileName, '_mex.', self.extension]);

                            % We save the timestamp corresponding to this
                            % compilation
                            fileProperties = dir([self.fileName, '_mex.', self.extension]);
                            timestamp = fileProperties.date;
                            save([self.fileName, '_timestamp_', mexext, '.mat'], 'timestamp')
                        end
                    catch
                        firstPartWorks = false;
                    end
                    % return to the previous path
                    cd(initialPath);
                    self.triedBefore = true;
                end
                if firstPartWorks
                    % The preparation worked, we try call the optimized method
                    self.isWorking = true;
                    try
                        command = {self.package, '.', self.fileName, '_mex(varargin{:});'};
                        command = cat(2, command{:});
                        results = eval(command);
                    catch
                        self.isWorking = false;
                    end
                else
                    self.isWorking = false;
                end
            end

            if ~self.isWorking
                % Inform that the method did not succeed
                results = replab.DispatchNext;
            end
        end
        
    end
    
end

