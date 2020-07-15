classdef CharacterTable < replab.Obj
% Describes the character table of a group

    properties (SetAccess = protected)
        group % replab.Group: group represented by character table
        irreps % cell of replab.Rep: irreducible representation
        classes % cell of replab.ConjugacyClass: conjugacy classes
        chars % matrix(nirreps,nclasses): matrix of characters
        table % replab.str.Table: table object for display purposes
    end
    
    methods

        function self = CharacterTable(group)
            self.group = group;
            % if we have Sn, use Andrew's function
%             if isa(group, 'replab.SymmetricGroup')
%                 [classes, irreps, chars] = ...
%             end
            [classes, irreps, chars] = self.generateTable;
            self.classes = classes;
            self.irreps = irreps;
            self.chars = round(chars); % add precise character method later
            colnames = cellfun(@(v) v.representative, self.classes, 'UniformOutput', false);
            nirreps = length(self.irreps);
            rownames = cellfun(@(n)['X.', num2str(n)], num2cell(1:nirreps), 'UniformOutput', false);
            self.table = replab.str.Table(self.chars);
            self.table.addRowNames(rownames);
            self.table.addColumnNames(colnames);
        end

        function s = headerStr(self)
            group_str = replab.headerStr(self.group);
            if length(group_str) > 1
                group_str = [lower(group_str(1)), group_str(2:end)];
            end
            s = sprintf(['Character table for ', group_str]);
        end
        
        function [classes, irreps, chars] = generateTable(self)
            ord = double(self.group.order);
            decomp = self.group.naturalRep.decomposition.nice;
            irreps = decomp.components;
            classes = self.group.conjugacyClasses;
            k = length(classes);
            ccreps = cell(1, length(classes));
            cclens = cell(1, length(classes));
            for i = 1:k
                ccreps{i} = classes{i}.representative;
                cclens{i} = classes{i}.size;
            end
            nirreps = length(irreps);
            chars = cell(nirreps, k);
            for i = 1:nirreps
                irrep = irreps{i};
                chars(i,:) = cellfun(@(x) trace(irrep.image(x)), ccreps, 'UniformOutput', false);
            end
            chars = cell2mat(chars);
            start_irreps = 1;
            while nirreps ~= k
                for i = start_irreps:nirreps
                    for j = start_irreps:nirreps
                        new_rep = irreps{i}.kron(irreps{j});
                        ss = 0;
                        for n = 1:k
                            char = trace(new_rep.image(ccreps{n}));
                            ss = ss + cclens{n}*(char)^2;
                        end
                        if ss ~= ord
                            new_irreps = new_rep.decomposition.nice;
                            start_irreps = nirreps + 1;
                            for m = 1:new_irreps.nComponents
                                new_irrep = new_irreps.component(m);
                                new_char = cellfun(@(x) trace(new_irrep.image(x)), ...
                                                    ccreps, 'UniformOutput', false);
                                new_char = cell2mat(new_char);
                                if ~any(all(abs(chars - new_char) < 1/2, 2))
                                    irreps{nirreps + 1} = new_irrep;
                                    chars(nirreps + 1, :) = new_char;
                                    nirreps = nirreps + 1;
                                end
                            end
                        end
                    end
                end
            end
        end

        
        function sym = mulliken(rep)
            % returns the mulliken symbol of a representation
        end
        
        function setIrrepLabels(self, labels)
        % set the labels for the irreducible representation in the table
        %
        % Args:
        %   labels (cell(1,nirreps)): cell array of character string labels
            self.table.addRowNames(labels)
        end
        
        function setClassLabels(self, labels)
        % set the labels for the conjugacy classes in the table
        %
        % Args:
        %   labels (cell(1,nclasses)): cell array of character string labels
            self.table.addColumnNames(labels)
        end
        
        
    end
end