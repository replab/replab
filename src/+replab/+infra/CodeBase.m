classdef CodeBase < replab.Str
% Describes a code base
%
% Terminology:
%
% - ``element``: Generic term for a subobject (class in a package, method in a class, etc)
% - ``member``: Methods/properties of a class, seen across the inheritance tree
% - ``parent``, ``children``: Package had subpackage children, 
%                             with functions/classes as children, classes have methods/properties children
% - ``superclass``, ``subclass``: Terminology used in inheritance hierarchies
% 
% 
    properties
        rootFolder % charstring: Path to the root folder
        packages % struct-based hash map
    end
    
    properties% (Access = protected)
        subpackages_ % struct-based hash map: Maps package paths to a cell array of their subpackages
        subclasses_ % struct-based hash map: Maps class paths to a cell array of their subclasses
    end
    
    methods

        function self = CodeBase(rootFolder, packageData)
            self.rootFolder = rootFolder;
            packages = struct;
            for i = 1:length(packageData)
                pkg = replab.infra.Package(self, packageData{i});
                id = replab.infra.shmEncode(pkg.path);
                packages.(id) = pkg;
            end
            self.packages = packages;
        end
        
        function p = package(self, varargin)
            id = replab.infra.shmEncode(varargin);
            if isfield(self.packages, id)
                p = self.packages.(id);
            else
                p = [];
            end
        end
        
        function c = subpackages(self, package)
        % Returns the direct subpackages of the given package
        %
        % Args:
        %   package (`.Package`): Package
        %
        % Returns:
        %   row cell vector of `.Package`: Subpackages of the given package
            if isempty(self.subpackages_)
                s = struct;
                pkgNames = fieldnames(self.packages);
                for i = 1:length(pkgNames)
                    pkg = self.packages.(pkgNames{i});
                    pp = pkg.path;
                    if ~isempty(pp)
                        parentPath = pp(1:end-1);
                        s = replab.infra.shmUpdate(s, parentPath, @(c) horzcat(c, {pkg}), {pkg});
                    end
                end
                self.subpackages_ = s;
            end
            c = replab.infra.shmLookup(self.subpackages_, package.path, {});
        end
        
        function p = allPackages(self)
        % Returns all packages present in this code base
            pkgNames = fieldnames(self.packages);
            p = cellfun(@(name) self.packages.(name), pkgNames, 'uniform', 0);
        end
        
        function f = allFunctions(self)
        % Returns all functions present in this code base
            p = self.allPackages;
            f = {};
            for i = 1:length(p)
                f = horzcat(f, p{i}.ownFunctions);
            end
        end
        
        function c = allClasses(self)
        % Returns all classes present in this code base
            p = self.allPackages;
            c = {};
            for i = 1:length(p)
                c = horzcat(c, p{i}.ownClasses);
            end
        end
        
        function c = subclasses(self, cls)
        % Returns the direct subclasses of the given class
        %
        % Args:
        %   cls (`.Class`): Class
        %
        % Returns:
        %   row cell vector of `.Class`: Subclasses of the given class
            if isempty(self.subclasses_)
                s = struct;
                classes = self.allClasses;
                for i = 1:length(classes)
                    cl = classes{i};
                    sup = cl.ownSuperclasses;
                    for j = 1:length(sup)
                        s = replab.infra.shmUpdate(s, sup{j}.path, @(c) horzcat(c, {cl}), {cl});
                    end
                end
                self.subclasses_ = s;
            end
            c = replab.infra.shmLookup(self.subclasses_, cls.path, {});
        end
        
        function e = getIdentifier(self, id)
            path = strsplit(id, '.');
            e = self.get(path{:});
        end
        
        function e = get(self, varargin)
            e = self.root.get(varargin{:});
        end
        
        function p = root(self)
            p = self.packages.(replab.infra.shmEncode({}));
        end
        
% $$$        function subpackageNames = subPackagesNames(self, packageNameParts)
% $$$             fn = fieldnames(self.packages);
% $$$             subpackageNames = {};
% $$$             if length(packageNameParts) == 0
% $$$                 for i = 1:length(fn)
% $$$                     if all(fn{i} ~= '_')
% $$$                         subpackageNames{1, end+1} = fn{i};
% $$$                     end
% $$$                 end
% $$$             else
% $$$                 pn = strjoin(packageNameParts, '_');
% $$$                 for i = 1:length(fn)
% $$$                     fni = fn{i};
% $$$                     if length(fni) > length(pn) && isequal(fni(1:length(pn)), pn)
% $$$                         rest = fni(length(pn)+2:end);
% $$$                         if all(rest ~= '_')
% $$$                             subpackageNames{1, end+1} = rest;
% $$$                         end
% $$$                     end
% $$$                 end
% $$$             end
% $$$         end
        
% $$$         function p = lookupPackage(self, nameParts)
% $$$         % Looks for a package from its name parts
% $$$         %
% $$$         % Args:
% $$$         %   nameParts (cell row vector of charstring): Parts of the package name
% $$$         %
% $$$         % Returns:
% $$$         %   :class:`+replab.+infra.Package`: The corresponding package, or ``[]`` if not found
% $$$             fname = replab.infra.OldCodeBase.fieldName(nameParts);
% $$$             if isfield(self.packages, fname)
% $$$                 p = getfield(self.packages, fname);
% $$$             else
% $$$                 p = [];
% $$$             end
% $$$         end
        
% $$$         function [package packageNameParts restNameParts] = lookupGreedy(self, nameParts)
% $$$         % Greedily looks up for a package from its name parts, disregarding the suffix that does not match 
% $$$         %
% $$$         % Returns `package` which matches `packageNameParts`, and `restNameParts` such that
% $$$         % ``nameParts = horzcat(packageNameParts, restNameParts)``.
% $$$         %
% $$$         % Args:
% $$$         %   nameParts (cell row vector of charstring): Parts of a fully qualified object name
% $$$         %
% $$$         % Returns
% $$$         % -------
% $$$         %   package: 
% $$$         %     :class:`+replab.+infra.Package`: The package looked up for
% $$$         %   packageNameParts:
% $$$         %     cell row vector of charstring: The maximal prefix of `nameParts` that matches a package name
% $$$         %   restNameParts:
% $$$         %     cell row vector of charstring: The tail of `nameParts` that could not be matched
% $$$             n = length(nameParts);
% $$$             for i = n:-1:0
% $$$                 packageNameParts = nameParts(1:i);
% $$$                 restNameParts = nameParts(i+1:n);
% $$$                 package = self.lookupPackage(packageNameParts);
% $$$                 if ~isempty(package)
% $$$                     return
% $$$                 end
% $$$             end
% $$$             error('Should not happen: empty name parts match the root package');
% $$$         end
        
% $$$         function writeDocTests(self, doctestPath)
% $$$         % Writes the doc tests of the whole code base in the specified folder
% $$$         %
% $$$         % Args:
% $$$         %   doctestPath (charstring): Path of the doctests generated code, must exist
% $$$             names = fieldnames(self.packages);
% $$$             for i = 1:length(names)
% $$$                 package = self.packages.(names{i});
% $$$                 fprintf('Writing tests for package %s:\n', package.fullName);
% $$$                 memberNames = fieldnames(package.members);
% $$$                 for j = 1:length(memberNames)
% $$$                     fprintf('.. member %s:\n', memberNames{j});
% $$$                     replab.infra.writeDocTests(doctestPath, package.member(memberNames{j}));
% $$$                 end
% $$$             end
% $$$         end
% $$$         
% $$$         function writeEnrichedSource(self, docSrcPath)
% $$$         % Writes the enriched source with the TOC elements
% $$$         %
% $$$         % Args:
% $$$         %   docSrcPath (charstring): Path of the enriched source, folder must exist, without trailing separator
% $$$             names = fieldnames(self.packages);
% $$$             for i = 1:length(names)
% $$$                 package = self.packages.(names{i});
% $$$                 fprintf('Writing enriched source for package %s:\n', package.fullName);
% $$$                 memberNames = fieldnames(package.members);
% $$$                 for j = 1:length(memberNames)
% $$$                     fprintf('.. member %s:\n', memberNames{j});
% $$$                     replab.infra.writeEnrichedSource(self, docSrcPath, package.member(memberNames{j}));
% $$$                 end
% $$$             end
% $$$         end

        
    end
   
    methods (Static)
                       
        function c = crawl(rootFolder)
        % Crawls the RepLAB source repository
        %
        % Args:
        %   rootFolder (charstring): Absolute path of the source directory (usually '$REPLAB_ROOT/src')
            packageData = {};
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
                path = fullfile(rootFolder, subpath{:});
                children = dir(path);
                ownFunctions = {};
                ownClasses = {};
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
                        filename = fullfile(rootFolder, subpath{:}, name);
                        disp(filename);
                        data = replab.infra.CodeTokens.fromFile(filename).parse;
                        switch class(data)
                          case 'replab.infra.FunctionLikeData'
                            ownFunctions{1,end+1} = data;
                          case 'replab.infra.ClassData'
                            ownClasses{1,end+1} = data;
                          otherwise
                            error(sprintf('Unknown type %s', class(data)));
                        end
                    end
                end
                packageData{1,end+1} = replab.infra.PackageData(packageNameParts, ownFunctions, ownClasses);
            end
            c = replab.infra.CodeBase(rootFolder, packageData);
        end
        
    end
    
end
