classdef Class < replab.infra.SourceElement

    properties
        superclassIdentifiers
        ownElements
    end

    properties (Access = protected)
        inheritedElements_
    end
    
    methods
        
        function self = Class(codeBase, package, classData)
            startLineNumber = 1;
            self = self@replab.infra.SourceElement(codeBase, package, classData.name, startLineNumber, classData.name, ...
                                                   classData.docLines, classData.docLineNumbers);
            sci = {};
            for i = 1:length(classData.superclassIdentifiers)
                id = classData.superclassIdentifiers{i};
                if ~isequal(id, 'handle')
                    sci{1,end+1} = id;
                end
            end
            self.superclassIdentifiers = sci;
            oe = struct;
            for i = 1:length(classData.ownMethods)
                md = classData.ownMethods{i};
                kind = 'method';
                m = replab.infra.ConcreteClassElement(codeBase, package, self, md.name, md.declarationLineNumber, ...
                                                      kind, md.declaration, md.attributes, ...
                                                      md.docLines, md.docLineNumbers);
                oe.(m.name) = m;
            end
            for i = 1:length(classData.ownProperties)
                pd = classData.ownProperties{i};
                kind = 'property';
                p = replab.infra.ConcreteClassElement(codeBase, package, self, pd.name, pd.declarationLineNumber, ...
                                                      kind, [], pd.attributes, ...
                                                      pd.docLines, pd.docLineNumbers);
                oe.(p.name) = p;
            end
            self.ownElements = oe;
        end
        
        % replab.infra.Element
        
        function c = childrenNames(self)
            c = vertcat(fieldnames(self.ownElements), fieldnames(self.inheritedElements));
            c = c(:).';
        end
        
        function e = lookup(self, id)
            if isfield(self.ownElements, id)
                e = self.ownElements.(id);
                return
            end
            ie = self.inheritedElements;
            if isfield(ie, id)
                e = ie.(id);
                return
            end
        end

        % replab.infra.SourceElement
        
        function str = name(self)
            str = self.sourceIdentifier;
        end
        
        function p = elementPath(self)
            p = {self.name};
        end
        
        % Own methods
        
        function ae = allElements(self)
        % Returns a struct whose fields contain 
            ae = replab.infra.shmMerge(self.ownElements, self.inheritedElements);
        end
        
        function am = allMethods(self)
            am = horzcat(self.ownMethods, self.inheritedMethods);
        end
        
        function ap = allProperties(self)
            ap = horzcat(self.ownProperties, self.inheritedProperties);
        end

        function om = ownMethods(self)
            oe = self.ownElements;
            om = oe(cellfun(@(x) isequal(x.kind, 'method'), oe));
        end
        
        function op = ownProperties(self)
            oe = self.ownElements;
            op = oe(cellfun(@(x) isequal(x.kind, 'property'), oe));
        end
        
        function im = inheritedMethods(self)
            ie = self.inheritedElements;
            im = oe(cellfun(@(x) isequal(x.kind, 'method'), ie));
        end
        
        function op = inheritedProperties(self)
            ie = self.inheritedElements;
            ip = ie(cellfun(@(x) isequal(x.kind, 'property'), ie));
        end
        
        function ie = inheritedElements(self)
            if isempty(self.inheritedElements_)
                ie = struct;
                already = struct;
                asc = self.allSuperclasses;
                for i = 1:length(asc)
                    sup = asc{i};
                    names = fieldnames(sup.ownElements);
                    mask = cellfun(@(n) ~isfield(already, n), names);
                    names = names(mask);
                    for j = 1:length(names)
                        name = names{j};
                        already.(name) = true;
                        se = sup.ownElements.(name);
                        iel = replab.infra.InheritedClassElement(self.codeBase, self, name, se.kind, se.attributes);
                        ie.(name) = iel;
                    end
                end
                self.inheritedElements_ = ie;
            end
            ie = self.inheritedElements_;
        end

        function asc = allSuperclasses(self)
        % Returns all superclasses of this class, not including itself
            asc = {};
            visitedClassIds = struct;
            queue = self.ownSuperclasses;
            while ~isempty(queue)
                c = queue{1};
                queue = queue(2:end);
                id = c.shmIdentifier;
                if ~isfield(visitedClassIds, id)
                    visitedClassIds.(id) = true;
                    asc{1,end+1} = c;
                    queue = horzcat(queue, c.ownSuperclasses);
                end
            end
        end
        
        function asc = allSubclasses(self)
        % Returns all subclasses of this class, not including itself
            asc = {};
            visitedClassIds = struct;
            queue = self.ownSubclasses;
            while ~isempty(queue)
                c = queue{1};
                queue = queue(2:end);
                id = c.shmIdentifier;
                if ~isfield(visitedClassIds, id)
                    visitedClassIds.(id) = true;
                    asc{1,end+1} = c;
                    queue = horzcat(queue, c.ownSubclasses);
                end
            end
        end
        
        function c = ownSubclasses(self)
            c = self.codeBase.subclasses(self);
        end
        
        function c = ownSuperclasses(self)
            c = cellfun(@(id) self.codeBase.getIdentifier(id), self.superclassIdentifiers, 'uniform', 0);
        end
        
    end
    
end
