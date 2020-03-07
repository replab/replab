classdef CodeBase < replab.Str
% Describes a code base
%
% Terminology:
%
% - ``element``: Generic term for a subobject (class in a package, method in a class, etc)
% - ``member``: Methods/properties of a class, seen across the inheritance tree
% - ``parent``, ``children``: Package had subpackage children, with functions/classes as children, classes have methods/properties children
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
                id = replab.infra.shm.encode(pkg.path);
                packages.(id) = pkg;
            end
            self.packages = packages;
        end

        function p = package(self, varargin)
            id = replab.infra.shm.encode(varargin);
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
        %   cell(\*,1) of `.Package`: Subpackages of the given package
            if isempty(self.subpackages_)
                s = struct;
                pkgNames = fieldnames(self.packages);
                for i = 1:length(pkgNames)
                    pkg = self.packages.(pkgNames{i});
                    pp = pkg.path;
                    if ~isempty(pp)
                        parentPath = pp(1:end-1);
                        s = replab.infra.shm.update(s, parentPath, @(c) horzcat(c, {pkg}), {pkg});
                    end
                end
                self.subpackages_ = s;
            end
            c = replab.infra.shm.look_up(self.subpackages_, package.path, {});
        end

        function p = allPackages(self)
        % Returns all packages present in this code base
        %
        % Returns:
        %   cell(\*,1) of `.Package`: Packages
            pkgNames = fieldnames(self.packages);
            p = cellfun(@(name) self.packages.(name), pkgNames, 'uniform', 0);
        end

        function e = allSourceElements(self)
        % Returns all source elements (functions, classes) in this code base
        %
        % Returns:
        %   cell(\*,1) of `+replab.+infra.SourceElement`: Source elements
            p = self.allPackages;
            e = {};
            for i = 1:length(p)
                e = horzcat(e, p{i}.ownClasses, p{i}.ownFunctions);
            end
        end

        function f = allFunctions(self)
        % Returns all functions present in this code base
        %
        % Returns:
        %   cell(\*,1) of `.Function`: All functions
            p = self.allPackages;
            f = {};
            for i = 1:length(p)
                f = horzcat(f, p{i}.ownFunctions);
            end
        end

        function c = allClasses(self)
        % Returns all classes present in this code base
        %
        % Returns:
        %   cell(\*,1) of `.Class`: All classes
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
        %   cell(1,\*) of `.Class`: Subclasses of the given class
            if isempty(self.subclasses_)
                s = struct;
                classes = self.allClasses;
                for i = 1:length(classes)
                    cl = classes{i};
                    sup = cl.ownSuperclasses;
                    for j = 1:length(sup)
                        s = replab.infra.shm.update(s, sup{j}.path, @(c) horzcat(c, {cl}), {cl});
                    end
                end
                self.subclasses_ = s;
            end
            c = replab.infra.shm.look_up(self.subclasses_, cls.path, {});
        end

        function e = getIdentifier(self, id)
            path = strsplit(id, '.');
            e = self.get(path{:});
        end

        function e = get(self, varargin)
            e = self.root.get(varargin{:});
        end

        function p = root(self)
            p = self.packages.(replab.infra.shm.encode({}));
        end

    end

end
