classdef Class < replab.infra.SourceElement
% Describes a class in the code base

    properties
        superclassIdentifiers % cell(1,*) of charstring: Identifiers of superclasses
        ownElementsStruct % struct: Elements of this class of type `.ClassElement`
        propertyLines % integer(1,*): Line numbers of properties
    end

    properties (Access = protected)
        ownDocumentationElements_
        inheritedDocumentationElements_
        inheritedElementsStruct_
        allSuperclasses_
        allSubclasses_
    end

    methods

        function self = Class(codeBase, package, classData)
            startLineNumber = 1;
            self = self@replab.infra.SourceElement(codeBase, package, classData.name, startLineNumber, classData.name, ...
                                                   classData.docLines, classData.docLineNumbers, false);
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
                if ~isequal(md.name, classData.name)
                    kind = 'method';
                    m = replab.infra.ConcreteClassElement(codeBase, package, self, md.name, md.declarationLineNumber, ...
                                                          kind, md.declaration, md.attributes, ...
                                                          md.docLines, md.docLineNumbers);
                    oe.(m.name) = m;
                end
            end
            for i = 1:length(classData.ownProperties)
                pd = classData.ownProperties{i};
                kind = 'property';
                p = replab.infra.ConcreteClassElement(codeBase, package, self, pd.name, pd.declarationLineNumber, ...
                                                      kind, [], pd.attributes, ...
                                                      pd.docLines, pd.docLineNumbers);
                oe.(p.name) = p;
            end
            self.ownElementsStruct = oe;
            self.propertyLines = classData.propertyLines;
        end

        %% replab.infra.Element

        function c = childrenNames(self)
            c = vertcat(fieldnames(self.ownElementsStruct), fieldnames(self.inheritedElementsStruct));
            c = c(:).';
        end

        function e = lookup(self, id)
            if isfield(self.ownElementsStruct, id)
                e = self.ownElementsStruct.(id);
                return
            end
            ie = self.inheritedElementsStruct;
            if isfield(ie, id)
                e = ie.(id);
                return
            end
            e = [];
        end

        %% replab.infra.SourceElement

        function p = elementPath(self)
            p = {self.name};
        end

        %% Own methods

        function ae = allElements(self)
        % Returns all elements of this class (including inherited) as a (sorted by name) row cell vector of `.ClassElement`
            ae = struct2cell(orderfields(self.allElementsStruct));
            ae = ae(:).';
        end

        function ae = allElementsStruct(self)
        % Returns a struct whose fields contain all elements (including inherited)
            ae = replab.infra.shm.merge2(self.ownElementsStruct, self.inheritedElementsStruct);
        end

        function am = allMethods(self)
        % Returns a row cell vector containing all methods, inc. inherited, sorted by name
            ae = self.allElements;
            am = ae(cellfun(@(x) x.isMethod, ae));
        end

        function ap = allProperties(self)
        % Returns a row cell vector containing all properties, inc. inherited, sorted by name
            ae = self.allElements;
            ap = ae(cellfun(@(x) x.isProperty, ae));
        end

        function oe = ownDocumentationElements(self)
        % Returns all elements that have documentation in this class
            if isempty(self.ownDocumentationElements_)
                self.ownDocumentationElements_ = self.computeOwnDocumentationElements;
            end
            oe = self.ownDocumentationElements_;
        end

        function op = ownDocumentationProperties(self)
            oe = self.ownDocumentationElements;
            op = oe(cellfun(@(x) x.isProperty, oe));
        end

        function om = ownDocumentationMethods(self)
            oe = self.ownDocumentationElements;
            om = oe(cellfun(@(x) x.isProperty, oe));
        end

        function ie = inheritedDocumentationElements(self)
        % Returns all elements that have documentation in a parent class
            if isempty(self.inheritedDocumentationElements_)
                self.inheritedDocumentationElements_ = self.computeInheritedDocumentation;
            end
            ie = self.inheritedDocumentationElements_;
        end

        function ip = inheritedDocumentationProperties(self)
            ie = self.inheritedDocumentationElements;
            ip = ip(cellfun(@(x) x.isProperty, ie));
        end

        function im = inheritedDocumentationMethods(self)
            ie = self.inheritedDocumentationElements;
            im = im(cellfun(@(x) x.isMethod, ie));
        end

        function oe = ownElements(self)
        % Returns all elements declared in this class as a (sorted by name) row cell vector of `.ClassElement`
            oe = struct2cell(orderfields(self.ownElementsStruct));
            oe = oe(:).';
        end

        function om = ownMethods(self)
            oe = self.ownElements;
            om = oe(cellfun(@(x) x.isMethod, oe));
        end

        function op = ownProperties(self)
            oe = self.ownElements;
            op = oe(cellfun(@(x) x.isProperty, oe));
        end

        function ie = inheritedElementsStruct(self)
            if isempty(self.inheritedElementsStruct_)
                ie = struct;
                already = struct;
                names = fieldnames(self.ownElementsStruct);
                for i = 1:length(names)
                    name = names{i};
                    already.(name) = true;
                end
                asc = self.allSuperclasses;
                for i = 1:length(asc)
                    sup = asc{i};
                    names = fieldnames(sup.ownElementsStruct);
                    mask = cellfun(@(n) ~isfield(already, n), names);
                    names = names(mask);
                    for j = 1:length(names)
                        name = names{j};
                        already.(name) = true;
                        se = sup.ownElementsStruct.(name);
                        iel = replab.infra.InheritedClassElement(self.codeBase, self, name, se.kind, se.attributes);
                        ie.(name) = iel;
                    end
                end
                self.inheritedElementsStruct_ = ie;
            end
            ie = self.inheritedElementsStruct_;
        end

        function ie = inheritedElements(self)
            ie = struct2cell(orderfields(self.inheritedElementsStruct));
            ie = ie(:).';
        end

        function im = inheritedMethods(self)
            ie = self.inheritedElements;
            im = ie(cellfun(@(x) x.isMethod, ie));
        end

        function ip = inheritedProperties(self)
            ie = self.inheritedElements;
            ip = ie(cellfun(@(x) x.isProperty, ie));
        end

        function c = ownSubclasses(self)
        % Returns all direct subclasses of this class
        %
        % Returns:
        %   row cell vector of `.Class`: Subclasses
            c = self.codeBase.subclasses(self);
        end

        function c = ownSuperclasses(self)
        % Returns all direct superclasses of this class
        %
        % Returns:
        %   row cell vector of `.Class`: Superclasses
            c = cellfun(@(id) self.codeBase.getIdentifier(id), self.superclassIdentifiers, 'uniform', 0);
        end

        function asc = allSuperclasses(self)
        % Returns all superclasses of this class, not including itself
            if isempty(self.allSuperclasses_)
                self.allSuperclasses_ = self.computeAllSuperclasses;
            end
            asc = self.allSuperclasses_;
        end

        function asc = allSubclasses(self)
        % Returns all subclasses of this class, not including itself
            if isempty(self.allSubclasses_)
                self.allSubclasses_ = self.computeAllSubclasses;
            end
            asc = self.allSubclasses_;
        end

    end

    methods (Access = protected)

        function asc = computeAllSuperclasses(self)
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

        function asc = computeAllSubclasses(self)
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

        function oe = computeOwnDocumentationElements(self)
            ae = self.allElements;
            mask = cellfun(@(el) isa(el, 'replab.infra.ConcreteClassElement') && ~el.doc.isempty, ae);
            oe = ae(mask);
        end

        function ie = computeInheritedDocumentation(self)
            ae = self.allElements;
            mask = cellfun(@(el) isa(el, 'replab.infra.ConcreteClassElement') && ~el.doc.isempty, ae);
            ie = ae(~mask);
        end

    end

end
