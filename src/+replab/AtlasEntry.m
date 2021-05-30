classdef AtlasEntry < replab.Obj
% Describes an Atlas entry available in the replab/atlas directory
%
% Stores
% Note: the JSON file must contain only ASCII characters; check if our parser `+replab.+util.parseJSON` for support
% of escape sequences.

    properties (SetAccess = protected)
        filename % (charstring): Filename, excluding path, including ``.json`` extension, omitted if not read from a file
        md5 % (charstring): MD5 hash in hexadecimal notation, without spaces and letters uppercase
        json % (charstring): File contents, omitted if not read from a file
        order % (vpi): Group order
        derivedSeriesAbelianInvariants % (cell(1,\*) of integer(1,\*)): Abelian invariants of the group derived series
    end

    methods (Static)

        function A = fromCache(s)
        %
        % Struct must have the md5
            error('TODO');
        end

        function A = fromFileContents(filename, json)
            data = replab
        end

        function A = fromAbstractGroup(name, group)
            error('TODO');
        end

    end

    methods (Access = protected)

        function group = computeGroup(self)
            data = replab.util.parseJSON(json);
            group = replab.atlas.parseGroup(data.group);
            if isfield(data, 'realCharacterTable')
                realCharacterTable = replab.atlas.parseCharacterTable(group, 'R', data.realCharacterTable);
                group.cache('realCharacterTable', realCharacterTable, 'error');
            end
            if isfield(data, 'complexCharacterTable')
                complexCharacterTable = replab.atlas.parseCharacterTable(group, 'C', data.complexCharacterTable);
                group.cache('complexCharacterTable', complexCharacterTable, 'error');
            end
            [name, derivedSeriesAbelianInvariants] = replab.AtlasEntry.additionalInfo(data, group);
            replab.AtlasEntry(filename, json, name, group.order, derivedSeriesAbelianInvariants, 'group', group);
        end

    end

    methods

        function self = AtlasEntry(filename, json, md5, order, derivedSeriesAbelianInvariants, varargin)
            args = struct('group', [], 'data', []);
            args = replab.util.populateStruct(args, varargin);
            self.filename = filename;
            self.json = json;
            self.md5 = md5;
            self.order = order;
            self.derivedSeriesAbelianInvariants = derivedSeriesAbelianInvariants;
            if ~isempty(args.group)
                self.cache('group', args.group, 'error');
            end
            if ~isempty(args.data)
                self.cache('data', args.data, 'error');
            end
        end

        function res = canMatch(self, group)
        % Returns whether this atlas entry can match the given group
            res = false;
            gAE = cellfun(@(d) d.abelianInvariants, group.derivedSeries, 'uniform', 0);
            myAE = cellfun(@(d) d.abelianInvariants, self.group.permutationGroup.derivedSeries, 'uniform', 0);
            if ~isequal(gAE, myAE)
                return
            end
            res = true;
            % TODO: tests based on conjugacy classes
        end

        function m = isomorphism(self, group)
            m = self.group.findIsomorphism(group);
        end

        function r = match(self, group)
            m = self.isomorphism(group);
            if isempty(m)
                r = [];
                return
            end
            r = replab.AtlasResult(group, self, m);
        end

        function laws = laws(self)
            laws = replab.laws.AtlasEntryLaws(self);
        end

    end

end
