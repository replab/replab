classdef AtlasEntry < replab.Obj
% Describes an Atlas entry available in the replab/atlas directory
%
% Stores
% Note: the JSON file must contain only ASCII characters; check if our parser `+replab.+util.parseJSON` for support
% of escape sequences.

    properties (SetAccess = protected)
        filename % (charstring): Filename, excluding path, including ``.json`` extension
        md5 % (charstring): MD5 hash in hexadecimal notation, without spaces and letters uppercase
        order % (vpi): Group order
        derivedSeriesAbelianInvariants % (cell(1,\*) of integer(1,\*)): Abelian invariants of the group derived series
    end

    methods (Static)

        function A = fromFile(filename)
        % Constructs an `.AtlasEntry` from a JSON file
        %
        % This method parses the JSON enough to read the group order and the abelian invariants. The group itself
        % is constructed lazily when requested.
        %
        % Args:
        %   filename (charstring): Filename, excluding path but including the ``.json`` extension
        %
        % Returns:
        %   `.AtlasEntry`: The parsed atlas entry
            json = fileread(fullfile(replab.globals.replabPath, 'atlas', filename));
            data = replab.util.parseJSON(json);
            derivedSeriesAbelianInvariants = cellfun(@cell2mat, data.group.derivedSeriesAbelianInvariants, 'uniform', 0);
            md5 = replab.util.md5(json);
            order = vpi(data.group.order);
            A = replab.AtlasEntry(filename, md5, order, derivedSeriesAbelianInvariants, 'json', json, 'data', data);
        end

    end

    methods (Access = protected)

        function G = computeGroup(self)
        % Computes the group stored in this entry from the known JSON data
            data = self.data;
            G = replab.atlas.parseGroup(data.group);
            if isfield(data, 'realCharacterTable')
                realCharacterTable = replab.atlas.parseCharacterTable(G, 'R', data.realCharacterTable);
                G.cache('realCharacterTable', realCharacterTable, 'error');
            end
            if isfield(data, 'complexCharacterTable')
                complexCharacterTable = replab.atlas.parseCharacterTable(G, 'C', data.complexCharacterTable);
                G.cache('complexCharacterTable', complexCharacterTable, 'error');
            end
        end

    end

    methods

        function self = AtlasEntry(filename, md5, order, derivedSeriesAbelianInvariants, varargin)
            args = struct('group', [], 'data', [], 'json', []);
            args = replab.util.populateStruct(args, varargin);
            self.filename = filename;
            self.md5 = md5;
            self.order = order;
            self.derivedSeriesAbelianInvariants = derivedSeriesAbelianInvariants;
            if ~isempty(args.group)
                self.cache('group', args.group, 'error');
            end
            if ~isempty(args.json)
                self.cache('json', args.json, 'error');
            end
            if ~isempty(args.data)
                self.cache('data', args.data, 'error');
            end
        end

        function J = json(self)
        % Returns the JSON data as text
        %
        % Returns:
        %   charstring: Data
            path = fullfile(replab.globals.replabPath, 'atlas', self.filename);
            data = self.cached('json', @() fileread(path));
        end

        function D = data(self)
        % Returns the parsed data
        %
        % Returns:
        %   struct or cell array: Parsed JSON data
            D = self.cached('data', @() replab.util.parseJSON(self.json));
        end

        function G = group(self)
        % Returns the group corresponding to this atlas entry
        %
        % Returns:
        %   `.AbstractGroup`: Group
            G = self.cached('group', @() self.computeGroup);
        end

        function res = canMatch(self, group)
        % Returns whether this atlas entry can match the given group
            res = false;
            x = self.derivedSeriesAbelianInvariants;
            y = cellfun(@(d) d.abelianInvariants, group.derivedSeries, 'uniform', 0);
            if length(x) ~= length(y)
                return
            end
            res = all(arrayfun(@(i) length(x{i}) == length(y{i}) && all(x{i} == y{i}), 1:length(x)));
            % TODO: tests based on conjugacy classes
        end

        function m = isomorphism(self, group)
        % Finds the isomorphism from the stored atlas group to the given group
        %
        % Returns:
        %   `.FiniteIsomorphism` or ``[]``: The isomorphism if it exists or ``[]``

            m = self.group.findIsomorphism(group);
        end

        function r = match(self, group)
        % Attemps to match the given group to this group
            m = self.isomorphism(group);
            if isempty(m)
                r = [];
                return
            end
            r = replab.AtlasResult(m);
        end

    end

    methods % Implementations

        function laws = laws(self)
            laws = replab.laws.AtlasEntryLaws(self);
        end

    end

end
