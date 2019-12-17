classdef Package < replab.infra.Element
    
    properties
        packagePath % row cell vector of charstring: Package path
        elements % struct-based hash map: Package elements that are not subpackages
        ownFunctions 
        ownClasses
    end

    methods
       
        function self = Package(codeBase, packageData)
        % Constructs a package instance
        %
        % Args:
        %   codeBase (`.CodeBase`): Code base this object is part of
        %   packageData (`.PackageData`): Data corresponding to this package
            if isempty(packageData.path)
                name = [];
            else
                name = packageData.path{end};
            end
            self = self@replab.infra.Element(codeBase, name);
            self.packagePath = packageData.path;
            elements = struct;
            ownFunctions = {};
            ownClasses = {};
            for i = 1:length(packageData.ownFunctions)
                fd = packageData.ownFunctions{i};
                fun = replab.infra.Function(codeBase, self, fd);
                ownFunctions{1,end+1} = fun;
                elements.(fd.name) = fun;
            end
            for i = 1:length(packageData.ownClasses)
                cd = packageData.ownClasses{i};
                cls = replab.infra.Class(codeBase, self, cd);
                ownClasses{1,end+1} = cls;
                elements.(cls.name) = cls;
            end
            self.elements = elements;
            self.ownFunctions = ownFunctions;
            self.ownClasses = ownClasses;
        end
        
        function [packagePath elementPath] = splitPath(self)
            packagePath = self.packagePath;
            elementPath = cell(1, 0);
        end
        
        function e = lookup(self, id)
            if isfield(self.elements, id)
                pkgel = self.elements.(id);
            else
                pkgel = [];
            end
            p = horzcat(self.path, {id});
            subpkg = self.codeBase.package(p{:});
            if ~isempty(pkgel) && ~isempty(subpkg)
                error('Package %s has both an element and a subpackage named %s', self.fullIdentifier, id);
            end
            e = [];
            if ~isempty(pkgel)
                e = pkgel;
            end
            if ~isempty(subpkg)
                e = subpkg;
            end
        end

        function c = childrenNames(self)
        % Returns the names of all direct children of this package
        %
        % This includes its subpackages, the classes and functions it contains.
            spn = cellfun(@(x) x.name, self.codeBase.subpackages(self), 'uniform', 0);
            fn = fieldnames(self.elements);
            fn = fn(:).';
            c = horzcat(spn, fn);
        end
        
        function c = children(self)
        % Returns all the direct children of this package
        %
        % Children includes its subpackages, the classes and functions it contains.
            c = horzcat(self.codeBase.subpackages(self), self.ownFunctions, self.ownClasses);
        end
        
    end
    
end
