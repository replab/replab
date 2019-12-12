classdef CodeBase < replab.Str
    
    properties
        rootDirectoryName % charstring: Path of the root directory
        packages % struct-based hash map
    end
    
    methods

        function self = CodeBase(packages)
            self.packages = packages;
        end

        function p = package(self, varargin)
            p = self.lookupPackage(varargin);
        end
        
        function names = hiddenFields(self)
            names = hiddenFields@replab.Str(self);
            names{1, end+1} = 'packages';
        end

        function [names values] = additionalFields(self)
            [names values] = additionalFields@replab.Str(self);
            fn = fieldnames(self.packages);
            for i = 1:length(fn)
                if isequal(fn{i}, 'root_')
                    nameParts = {};
                else
                    nameParts = strsplit(fn{i}, '_');
                end
                args = strjoin(cellfun(@(x) ['''' x ''''], nameParts, 'uniform', 0), ', ');
                names{1, end+1} = sprintf('package(%s)', args);
                values{1, end+1} = self.packages.(fn{i});
            end
        end
        
        function subpackageNames = subPackagesNames(self, packageNameParts)
            fn = fieldnames(self.packages);
            subpackageNames = {};
            if length(packageNameParts) == 0
                for i = 1:length(fn)
                    if all(fn{i} ~= '_')
                        subpackageNames{1, end+1} = fn{i};
                    end
                end
            else
                pn = strjoin(packageNameParts, '_');
                for i = 1:length(fn)
                    fni = fn{i};
                    if length(fni) > length(pn) && isequal(fni(1:length(pn)), pn)
                        rest = fni(length(pn)+2:end);
                        if all(rest ~= '_')
                            subpackageNames{1, end+1} = rest;
                        end
                    end
                end
            end
        end
        
        function p = lookupPackage(self, nameParts)
        % Looks for a package from its name parts
        %
        % Args:
        %   nameParts (cell row vector of charstring): Parts of the package name
        %
        % Returns:
        %   :class:`+replab.+infra.Package`: The corresponding package, or ``[]`` if not found
            fname = replab.infra.CodeBase.fieldName(nameParts);
            if isfield(self.packages, fname)
                p = getfield(self.packages, fname);
            else
                p = [];
            end
        end
        
        function [package packageNameParts restNameParts] = lookupGreedy(self, nameParts)
        % Greedily looks up for a package from its name parts, disregarding the suffix that does not match 
        %
        % Returns `package` which matches `packageNameParts`, and `restNameParts` such that
        % ``nameParts = horzcat(packageNameParts, restNameParts)``.
        %
        % Args:
        %   nameParts (cell row vector of charstring): Parts of a fully qualified object name
        %
        % Returns
        % -------
        %   package: 
        %     :class:`+replab.+infra.Package`: The package looked up for
        %   packageNameParts:
        %     cell row vector of charstring: The maximal prefix of `nameParts` that matches a package name
        %   restNameParts:
        %     cell row vector of charstring: The tail of `nameParts` that could not be matched
            n = length(nameParts);
            for i = n:-1:0
                packageNameParts = nameParts(1:i);
                restNameParts = nameParts(i+1:n);
                package = self.lookupPackage(packageNameParts);
                if ~isempty(package)
                    return
                end
            end
            error('Should not happen: empty name parts match the root package');
        end
        
    end
   
    methods (Static)
        
        function fn = fieldName(nameParts)
            if isempty(nameParts)
                fn = 'root_';
            else
                fn = strjoin(nameParts, '_');
            end
        end
        
        function c = crawl(rootDirectoryName)
            packages = struct;
            % toExplore represents a stack of subpaths to explore
            % toExplore is a row cell array, each element inside
            % is a cell array of char strings, which represent a
            % sequence of subfolders of pathStr
            toExplore = {{}};
            while length(toExplore) > 0
                % the current path explored
                subpath = toExplore{1};
                packageNameParts = cellfun(@(x) x(2:end), subpath, 'uniform', 0);
                toExplore = toExplore(2:end);
                % the path to explore
                path = fullfile(rootDirectoryName, subpath{:});
                children = dir(path);
                members = {};
                for i = 1:length(children)
                    name = children(i).name;
                    if isequal(name, '.') || isequal(name, '..')
                        % do nothing
                    elseif children(i).isdir
                        % folder
                        assert(name(1) == '+', 'We only support crawling subpackages');
                        newsubpath = horzcat(subpath, {name});
                        toExplore{1,end+1} = newsubpath;
                    elseif isequal(name(end-1:end), '.m')
                        % is not a folder and has a Matlab file extension
                        filename = fullfile(rootDirectoryName, subpath{:}, name);
                        parseState = replab.infra.CodeParseState.fromFile(filename);
                        switch parseState.peek
                          case 'CLASSDEF'
                            member = replab.infra.Class.fromParseState(parseState, packageNameParts);
                          case 'FUNCTION'
                            member = replab.infra.Function.fromParseState(parseState, packageNameParts);
                          otherwise
                            error(['Unrecognized first tag ' parseState.peek]);
                        end
                        members{1,end+1} = member;
                    end
                end
                package = replab.infra.Package(packageNameParts, members);
                fname = replab.infra.CodeBase.fieldName(packageNameParts);
                packages.(fname) = package;
            end
            c = replab.infra.CodeBase(packages);
        end
        
    end
    
end
