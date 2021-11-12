classdef IdentifyPointGroup < replab.Str
% Gives information about a molecule
%
% References:
%   Beruski, O., Vidal, L. N. (2013). Algorithms for computer detection of symmetry elements
%   in molecular systems. Journal of Computational Chemistry, Volume 35, Issue 4.
%   https://doi-org/10.1002/jcc.23493
%
%   Vallance, C. (2020). Symmetry Classification of Molecules: Point Groups.
%   Chemistry LibreTexts. https://chem.libretexts.org


    properties
        principalAxis % (double(3, 1)): principal axis of rotation of the molecule
        ord % (integer): highest order rotation of the molecule that preserves symmetry
        pointGroup % (charstring): characters representing the point group
        table % (`+replab.ComplexCharacterTable`): the character table for the molecule's point group
        ndec % (integer): number of decimal places to check for equality
    end

    methods (Static)

        function ct = findCharacterTable(pointGroup)
        % Determines the character table of a molecule with a certain point group
        %
        % Args:
        %   pointGroup (charstring): name of point group (i.e. C3v or D4d)
        %
        % Returns:
        %   ct (`+replab.CharacterTable`): character table of point group
            pg = split(pointGroup, '');
            pg = pg(2:end-1);
            if pg{1} == 'T'
                if length(pg) > 1
                    if pg{2} == 'h'
                        % direct product C2 and alternating group 4
                        Ci = replab.ct.CyclicCharacterTable(2);
                        A4 = replab.ct.A4CharacterTable;
                        ct = Ci.directProduct(A4);
                    elseif pg{2} == 'd'
                        % symmetric group 4
                        ct = replab.ct.S4CharacterTable;
                    end
                else
                    % alternating group 4
                    ct = replab.ct.A4CharacterTable;
                end
            elseif pg{1} == 'O'
                % symmetric group 4
                if length(pg) > 1 && pg{2} == 'h'
                    % direct product C2 and S4
                    Ci = replab.ct.CyclicCharacterTable(2);
                    S4 = replab.ct.S4CharacterTable;
                    ct = Ci.directProduct(S4);
                else
                    ct = replab.ct.S4CharacterTable;
                end
            elseif pg{1} == 'I'
                % alternating group 5
                if length(pg) > 1 && pg{2} == 'h'
                    % direct product C2 and alternating group 5
                    Ci = replab.ct.CyclicCharacterTable(2);
                    A5 = replab.ct.A5CharacterTable;
                    ct = Ci.directProduct(A5);
                else
                    % alternating group 5
                    ct = replab.ct.A5CharacterTable;
                end
            elseif pg{1} == 'C'
                if length(pg) > 2
                    n = str2double(pg{2});
                    if pg{3} == 'v'
                        ct = replab.ct.DihedralCharacterTable(n);
                    elseif pg{3} == 'h'
                        Ci = replab.ct.CyclicCharacterTable(2);
                        Cn = replab.ct.CyclicCharacterTable(n);
                        ct = Ci.directProduct(Cn);
                    end
                else
                    if pg{2} == 's' || pg{2} == 'i'
                        ct = replab.ct.CyclicCharacterTable(2);
                    else
                        n = str2double(pg{2});
                        ct = replab.ct.CyclicCharacterTable(n);
                    end
                end
            elseif pg{1} == 'D'
                n = str2double(pg{2});
                if length(pg) > 2
                    if pg{3} == 'd'
                        ct = replab.ct.DihedralCharacterTable(2*n);
                    elseif pg{3} == 'h'
                        Ci = replab.ct.CyclicCharacterTable(2);
                        Dn = replab.ct.DihedralCharacterTable(n);
                        ct = Ci.directProduct(Dn);
                    end
                else
                    ct = replab.ct.DihedralCharacterTable(n);
                end
            elseif pg{1} == 'S'
                n = str2double(pg{2});
                ct = replab.ct.CyclicCharacterTable(n);
            end
        end

        function SEAs = groupAtoms(D, ndec)
        % Finds atoms of same type with same distances
        %
        % Args:
        %   D (double(natoms, natoms)): distance matrix
        %   ndec (integer): number of decimal places to test equality
        %
        % Returns
        %   SEAs (cell(1, \*) of integer(1, \*)): array of lists of same atoms with same distances
            s = round(sum(D, 2), ndec);
            u = unique(s);
            DSort = sort(D, ndec);
            SEAs = cell(1, length(u));
            i = 1;
            while i <= length(u)
                f = find(s == u(i));
                if ~all(all(iszero(round(DSort(f,:) - DSort(f(1), :), ndec))))
                    subu = unique(round(DSort(f, :), ndec), 'rows');
                    dimu = size(subu);
                    SEAs(end+1:end+1+dimu(1)-ndec) = {{}};
                    for j = 1:dimu(1)
                        subf = find(all(round(DSort(f, :), ndec) == subu(j,:), ndec));
                        SEAs{i} = subf;
                        i = i + 1;
                    end
                else
                    SEAs{i} = f;
                    i = i + 1;
                end
            end
        end

        function [I, rcm] = InertiaTensor(coords, masses)
        % calculates the inertia tensor for a group of atoms
        %
        % Args:
        %   coords (double(3, \*)): coordinates of atoms
        %   masses (double(1, \*)): masses of atoms
        %
        % Returns:
        %   I (double(3, 3)): inertia tensor of given coordinates
        %   rcm (double(1, 3)): center of mass vector
            rcm = sum(diag(masses)*coords, 1)/sum(masses);
            rcm_coords = coords - rcm;
            wcs = sum(diag(masses)*rcm_coords.^2, 1);
            I = diag([wcs(2)+wcs(3), wcs(1)+wcs(3), wcs(1)+wcs(2)]);
            for i = 1:3
                for j = i+1:3
                    Ixy = sum(-diag(masses)*(rcm_coords(:,i).*rcm_coords(:,j)), 1);
                    I(i,j) = Ixy;
                    I(j,i) = Ixy;
                end
            end
        end

        function [ok, inds] = Rotation(k, axis, coords, symbols, ndec)
        % Apply rotation of 2*pi/k to coords about an axis
        %
        % Args:
        %   k (integer): number of rotations which should keep coordinates the same
        %   axis (double(1, 3)): axis about which to rotate coordinates
        %   coords (double(3, natoms)): atomic coordinates
        %   symbols (cell(1, natoms) of charstring): atomic symbols of coords
        %   ndec (integer): number of decimals to round
        %
        % Returns:
        %   ok (logical): whether reflection keeps coordinates the same
        %   inds (integer(1, natoms)): action of rotation on coordinates
            theta = 2*pi/k;
            axis(round(axis, ndec) == 0) = 0;
            a = axis / norm(axis);
            R = [cos(theta)+a(1)^2*(1-cos(theta)), a(1)*a(2)*(1-cos(theta))-a(3)*sin(theta), a(1)*a(3)*(1-cos(theta))+a(2)*sin(theta);
                 a(1)*a(2)*(1-cos(theta))+a(3)*sin(theta), cos(theta)+a(2)^2*(1-cos(theta)), a(2)*a(3)*(1-cos(theta))-a(1)*sin(theta);
                 a(1)*a(3)*(1-cos(theta))-a(2)*sin(theta), a(2)*a(3)*(1-cos(theta))+a(1)*sin(theta), cos(theta)+a(3)^2*(1-cos(theta))];
            new_coords = R * coords;
            [ok, inds] = replab.IdentifyPointGroup.coordsEqual(coords, new_coords, symbols, ndec);
        end

        function [ok, inds] = Reflection(normal, coords, symbols, ndec)
        % Reflect coords through a plane with a given normal
        %
        % Args:
        %   normal (double(1, 3)): normal of plane through which to reflect
        %   coords (double(3, natoms)): atomic coordinates
        %   symbols (cell(1, natoms) of charstring): atomic symbols of coords
        %   ndec (integer): number of decimals to round
        %
        % Returns:
        %   ok (logical): whether reflection keeps coordinates the same
        %   inds (integer(1, natoms)): action of reflection on coordinates
            normal(round(normal, ndec) == 0) = 0;
            normal = normal / norm(normal); % normalize normal vector
            sigma = [1-2*normal(1)^2, -2*normal(1)*normal(2), -2*normal(1)*normal(3);
                     -2*normal(1)*normal(2), 1-2*normal(2)^2, -2*normal(2)*normal(3);
                     -2*normal(1)*normal(3), -2*normal(2)*normal(3), 1-2*normal(3)^2];
            new_coords = sigma * coords;
            [ok, inds] = replab.IdentifyPointGroup.coordsEqual(coords, new_coords, symbols, ndec);
        end

        function [ok, inds] = ImproperRotation(k, axis, coords, symbols, ndec)
        % Apply rotation of 2*pi/k and then reflection through rotation plane
        %
        % Args:
        %   k (integer): number of rotations then reflection which should keep coordinates the same
        %   axis (double(1, 3)): axis about which to rotate coordinates and normal for reflection plane
        %   coords (double(3, natoms)): atomic coordinates
        %   symbols (cell(1, natoms) of charstring): atomic symbols of coords
        %   ndec (integer): number of decimals to round
        %
        % Returns:
        %   ok (logical): whether reflection keeps coordinates the same
        %   inds (integer(1, natoms)): action of improper rotation on coordinates
            theta = 2*pi/k;
            axis(round(axis, ndec) == 0) = 0;
            a = axis / norm(axis);
            R = [cos(theta)+a(1)^2*(1-cos(theta)), a(1)*a(2)*(1-cos(theta))-a(3)*sin(theta), a(1)*a(3)*(1-cos(theta))+a(2)*sin(theta);
                 a(1)*a(2)*(1-cos(theta))+a(3)*sin(theta), cos(theta)+a(2)^2*(1-cos(theta)), a(2)*a(3)*(1-cos(theta))-a(1)*sin(theta);
                 a(1)*a(3)*(1-cos(theta))-a(2)*sin(theta), a(2)*a(3)*(1-cos(theta))+a(1)*sin(theta), cos(theta)+a(3)^2*(1-cos(theta))];
            rot_coords = R * coords;
            sigma = [1-2*axis(1)^2, -2*axis(1)*axis(2), -2*axis(1)*axis(3);
                     -2*axis(1)*axis(2), 1-2*axis(2)^2, -2*axis(2)*axis(3);
                     -2*axis(1)*axis(3), -2*axis(2)*axis(3), 1-2*axis(3)^2];
            new_coords = sigma * rot_coords;
            [ok, inds] = replab.IdentifyPointGroup.coordsEqual(coords, new_coords, symbols, ndec);
        end

        function [coords, i1] = sortCoords(coords, ndec)
        % Sorts coordinates in ascending order
        %
        % Convention: atoms are sorted by increasing x-coordinate, then
        %             y-coordinate, and then z-coordinate
        %
        % Args:
        %   coords (double(3, natoms)): set of coordinates
        %   ndec (integer): number of decimals to check equality
        %
        % Returns:
        %   coords (double(3, natoms)): sorted coordinates
        %   i1 (integer(1, natoms)): indices that revert coords to original order
            natoms = size(coords, 2);
            [~, i1] = sort(coords(1, :));
            coords = coords(:, i1);
            [~, v1] = unique(round(coords(1, :), ndec), 'first');
            degen1 = 1:natoms;
            degen1(v1) = [];
            while ~isempty(degen1)
                f1 = find(round(coords(1, :), ndec) == round(coords(1, degen1(1)), ndec));
                degen1(1:length(f1)-1) = [];
                [~, i2] = sort(coords(2, f1));
                i1(f1) = i1(f1(i2));
                coords(:, f1) = coords(:, f1(i2));

                [~, v2] = unique(round(coords(2, f1), ndec), 'first');
                degen2 = 1:length(f1);
                degen2(:, v2) = [];
                while ~isempty(degen2)
                    f2 = find(round(coords(2, f1), ndec) == round(coords(2, f1(degen2(1))), ndec));
                    degen2(1:length(f2)-1) = [];
                    [~, i3] = sort(coords(3, f1(f2)));
                    i1(f1(f2)) = i1(f1(f2(i3)));
                    coords(:, f1(f2)) = coords(:, f1(f2(i3)));
                end
            end
        end

        function [ok, inds] = coordsEqual(old_coords, new_coords, symbols, ndec)
        % Checks whether a new set of coordinates is equivalent to original ones
        %
        % Convention: assumes self.coords is already sorted
        %
        % Args:
        %   old_coords (double(3, natoms)): reference coordinates
        %   new_coords (double(3, natoms)): new set of coordinates to check for equality
        %   symbols (cell(1, natoms) of charstring): molecular symbols of coordinates
        %   ndec (integer): number of decimal places to round for equality testing
        %
        % Returns:
        %   ok (logical): whether new coordinates are equal to self.coords
        %   inds (integer(1, natoms)): vector that will sort old_coords to new_coords
            [new_coords, inds] = replab.IdentifyPointGroup.sortCoords(new_coords, ndec);
            ok = isequal(round(old_coords, ndec), round(new_coords, ndec));
            if ok
                ok = isequal(symbols, symbols(inds));
            end
        end

        function D = distanceMatrix(coords)
        % Determine the distance matrix of the molecules
        %
        % Convention: entries of row i are the distances of atom i from
        %             atoms 1, 2, ..., natoms so D(i, i) = 0
        %
        % Args:
        %   coords (double(3, natoms)): molecular coordinates
        %
        % Returns:
        %   D (double(natoms, natoms)): distance matrix
            natoms = size(coords, 2);
            D = zeros(natoms);
            for i = 1:natoms
                for j = 1:natoms
                    D(i,j) = sqrt(sum((coords(:,i)-coords(:,j)).^2));
                end
            end
        end

    end

    methods

        function self = IdentifyPointGroup(coords, symbols, masses, ndec)
        % Creates a molecule object
        %
        % Args:
        %   coords (double(3, natoms)): a set of coordinates
        %   symbols (cell(1, natoms) of charstring): atomic symbols of coordinates
        %   masses (cell(1, natoms) of charstring): atomic masses of coordinates
        %   ndec (integer): number of decimal places to round when checking for equality
            self.ndec = ndec;
            SEAs = self.groupAtoms(self.distanceMatrix(coords), ndec);
            self.findPointGroup(coords, symbols, masses, SEAs);
        end


        function [order, axis] = findRotations(self, coords, symbols, masses, inds)
        % Determine the rotation order and axis for a group of atoms
        %
        % From Beruski
        %
        % Args:
        %   coords (double(3, natoms)): coordinates of atoms
        %   symbols (cell(1, natoms) of charstring): atomic symbols of coordinates
        %   masses (cell(1, natoms) of charstring): atomic masses of coordinates
        %   inds (integer(1, \*)): indices of group of atoms whose rotations we find
        %
        % Returns:
        %   order (integer): highest order rotation possible
        %   axis (double(1, 3)): axis of rotation
            k = length(inds);
            if k < 3
                return
            end
            [I_atom, ~] = self.InertiaTensor(coords(:, inds)', masses(inds));
            [evecs, evals] = eig(I_atom);
            moments = sort(round(diag(evals), self.ndec));
            u = unique(moments);
            if abs(sum(moments(1:2)) - moments(3)) < self.ndec
                if length(u) == 2
                    % Ck rotations
                    n = k;
                else
                    % Ck/2 rotations or next biggest factor of Ck
                    n = floor(k/2);
                    while rem(k, n) ~= 0
                        n = n - 1;
                    end
                end
                % Must check rotations about each eigenvector (but usually 3rd evec)
                evec = 3;
                same = self.Rotation(n, evecs(:, evec), coords, symbols, self.ndec);
                if ~same && length(u) ~= 2
                    evec = 2;
                    same = self.Rotation(n, evecs(:, evec), coords, symbols, self.ndec);
                    if ~same
                        evec = 1;
                        same = self.Rotation(n, evecs(:, evec), coords, symbols, self.ndec);
                    end
                    if ~same
                        evec = 3;
                    end
                end
                while ~same
                    % Don't have rotation of order k so check factors of k
                    n = min(n - 1, floor(k/2));
                    while rem(k, n) ~= 0
                        n = n - 1;
                    end
                    same = self.Rotation(n, evecs(:, evec), coords, symbols, self.ndec);
                end
                order = n;
                axis = evecs(:, evec);
            else
                if length(u) == 3
                    % Have a rotation of order 2 along some eigenvector
                    order = 2;
                    i = 1;
                    same = self.Rotation(2, evecs(:, i), coords, symbols, self.ndec);
                    while ~same
                        i = 1 + 1;
                        same = self.Rotation(2, evecs(:, i), coords, symbols, self.ndec);
                    end
                    axis = evecs(:, i);
                else
                    % Have a rotation of at most k/2 along different eigenvector
                    n = k / 2;
                    if length(u) == 2
                        if moments(2) ~= u(1) && moments(3) ~= u(1)
                            naxis = 1;
                        elseif moments(1) ~= u(2) && moments(3) ~= u(2)
                            naxis = 2;
                        else
                            naxis = 3;
                        end
                    else
                        naxis = 3;
                    end
                    same = self.Rotation(n, evecs(:, naxis), coords, symbols, self.ndec);
                    while ~same
                        % Check rotations that are factors of k/2
                        n = min(n - 1, floor(k/2));
                        while rem(k, n) ~= 0
                            n = n - 1;
                        end
                        same = self.Rotation(n, evecs(:, naxis), coords, symbols, self.ndec);
                    end
                    order = n;
                    axis = evecs(:, naxis);
                end
            end
        end

        function [ok, axis] = findC2Axes(self, coords, symbols, SEAs, rcm)
        % Determines whether molecule has a C2 axis and returns first one found
        %
        % Args:
        %   coords (double(3, natoms)): atomic coordinates
        %   symbols (cell(1, natoms) of charstring): atomic symbols
        %   SEAs (cell(1, \*) of double(\*, 1)): cell array of vectors of the indices of similar atoms
        %   rcm (double(1, 3)): centre of mass vector
        %
        % Returns:
        %   ok (logical): whether molecule contains C2 axis
        %   axis (double(1, 3)): axis of C2 rotation
            axis = [0,0,0];
            ok = false;
            nSEA = cellfun(@length, SEAs);
            diatomic = find(nSEA == 2);
            if length(diatomic) > 1
                rij = coords(:, SEAs{diatomic(1)}(1)) - coords(:, SEAs{diatomic(1)}(2));
                rkl = coords(:, SEAs{diatomic(2)}(1)) - coords(:, SEAs{diatomic(2)}(2));
                rC2 = cross(rij, rkl);
                same = self.Rotation(2, rC2, coords, symbols, self.ndec);
                if same
                    ok = true;
                    axis = rC2;
                end
            else
                for i = 1:length(SEAs)
                    SEA = SEAs{i};
                    ok = false;
                    if length(SEA) > 2
                        continue
                    end
                    if length(SEA) == 2
                        mid = (coords(:, SEA(1)) + coords(:, SEA(2))) / 2;
                        if ~isequal(round(mid, self.ndec), round(rcm', self.ndec))
                            rC2 = mid - rcm';
                            same = self.Rotation(2, rC2, coords, symbols, self.ndec);
                            if same
                                ok = true;
                            end
                        end
                        if ~ok
                            rC2 = coords(:, SEA(2)) - rcm';
                            same = self.Rotation(2, rC2, coords, symbols, self.ndec);
                            if same
                                ok = true;
                            end
                        end
                    end
                    if ~ok
                        rC2 = coords(:, SEA(1)) - rcm';
                        same = self.Rotation(2, rC2, coords, symbols, self.ndec);
                        if same
                            ok = true;
                        end
                    end

                    if ok
                        axis = rC2;
                        break
                    end
                end
            end
        end

        function ok = findCnAxes(self, n, max_rotations, coords, symbols, SEAs)
        % Find whether there are multiple higher order rotation axes
        %
        % Returns:
        %   n (integer): order of rotation
        %   max_rotations (integer): check whether there are at least this many
        %                            rotation axes of order n
        %   coords (double(3, natoms)): atomic coordinates
        %   symbols (cell(1, natoms) of charstring): atomic symbols
        %   SEAs (cell(1, \*) of double(\*, 1)): cell array of vectors of the indices of similar atoms
        %
        % Returns
        %   ok (logical): whether there are at least max_rotations rotation axes of order n
            nrot = 0;
            ok = false;
            nSEA = cellfun(@length, SEAs);
            f = find(nSEA > n);
            for i = 1:length(f)
                locs = SEAs{f(i)};
                k = length(locs);
                for p = 2:k
                    for q = 1:p-1
                        ri = coords(:, locs(p));
                        rj = coords(:, locs(q));
                        rij = ri - rj;
                        for m = 1:q-1
                            rk = coords(:, locs(m));
                            rik = ri - rk;
                            if round(norm(rij), self.ndec) ~= round(norm(rik), self.ndec)
                                continue
                            end
                            same = self.Rotation(n, cross(rij, rik), coords, symbols, self.ndec);
                            if same
                                nrot = nrot + 1;
                                if nrot == max_rotations
                                    ok = true;
                                    return
                                end
                            end
                        end
                    end
                end
            end
        end

        function [ind, perpAxes] = findPerpendicularAxes(self, axes, ords, rcm, SEAs, coords, symbols)
        % Finds all C2 axes perpendicular to the principal axis
        %
        % Args:
        %   axes (double(3, \*)): any already known axes
        %   ords (integer(1, \*)): orders of known axes of rotation
        %   rcm (double(1, 3)): centre of mass vector
        %   SEAs (cell(1, \*) of double(\*, 1)): cell array of vectors of the indices of similar atoms
        %   coords (double(3, natoms)): atomic coordinates
        %   symbols (cell(1, natoms) of charstring): atomic symbols
        %
        % Returns:
        %   ind (integer): number of perpendicular axes
        %   perpAxes (double(3, \*)): axes of rotation perpendicular to principal axis
            n = self.ord;
            perpAxes = zeros(3, n);
            ind = 1;
            for m = 1:length(ords)
                if ords(m) == 2 && isequal(round(dot(self.principalAxis, axes(:,m)), self.ndec), 0)
                    perpAxes(:, ind) = axes(:, m);
                    ind = ind + 1;
                end
            end
            for i = 1:length(SEAs)
                SEA = SEAs{i};
                k = length(SEA);
                for p = 2:k
                    for q = 1:p-1
                        midpoint = (coords(:, SEA(p)) + coords(:, SEA(q))) / 2;
                        if isequal(round(midpoint, self.ndec), round(rcm', self.ndec))
                            break
                        end
                        rperp = (midpoint - rcm') / norm(midpoint - rcm');
                        if round(abs(dot(rperp, self.principalAxis)), self.ndec) == 1
                            break
                        end
                        f = all(abs(perpAxes) - abs(rperp) < 10^(-1*self.ndec));
                        if any(f)
                            if any(all(sign(perpAxes(:, f)) == sign(rperp))) || any(all(sign(perpAxes(:, f)) == -1*sign(rperp)))
                                break
                            end
                        end
                        if ~isequal(round(dot(self.principalAxis, rperp), self.ndec), 0)
                            break
                        end
                        ok = self.Rotation(2, rperp, coords, symbols, self.ndec);
                        if ok
                            perpAxes(:, ind) = rperp;
                            ind = ind + 1;
                        end
                    end
                end
            end
        end

        function nsigmad = findDihedralMirrorPlanes(self, coords, symbols, perpAxes)
        % Determines number of dihedral mirror planes in the molecule
        %
        % Args:
        %   coords (double(3, natoms)): atomic coordinates
        %   symbols (cell(1, natoms) of charstring): atomic symbols
        %   perpAxes (double(3, \*)): axes of rotation perpendicular to principal axis
        %
        % Returns:
        %   nsigmad (integer): number of dihedral mirror planes
            nsigmad = 0;
            n = self.ord;
            for p = 2:n
                for q = 1:p-1
                    mid = (perpAxes(:, p) + perpAxes(:, q)) / 2;
                    sigmad = cross(self.principalAxis, mid);
                    ok = self.Reflection(sigmad, coords, symbols, self.ndec);
                    if ok
                        nsigmad = nsigmad + 1;
                    end
                    % We get a second dihedral mirror plane from the same pair
                    % of axes (checking this may be redundant)
                    mid = (perpAxes(:, p) + -1*perpAxes(:, q)) / 2;
                    sigmad = cross(self.principalAxis, mid);
                    ok = self.Reflection(sigmad, coords, symbols, self.ndec);
                    if ok
                        nsigmad = nsigmad + 1;
                    end
                    if nsigmad == self.ord
                        break
                    end
                end
            end
        end

        function [planar, normal] = checkPlanar(self, coords)
        % Determines if molecule is planar
        %
        % Args:
        %   coords (double(3, natoms)): atomic coordinates
        %
        % Returns:
        %   planar (logical): whether molecule is planar
        %   normal (double(3, 1)): normal of plane formed by first three atoms
            natoms = size(coords, 2);
            normal = cross(coords(:, 1) - coords(:, 2), coords(:, 1) - coords(:, 3));
            planar = true;
            for i = min(natoms, 4):natoms
                if round(dot(coords(:, i), normal), self.ndec) ~= 0
                    planar = false;
                    break
                end
            end
        end

        function [horizontal, nvertical] = findReflections(self, coords, symbols, SEAs)
        % Finds reflection planes in the molecule
        %
        % Args:
        %   coords (double(3, natoms)): atomic coordinates
        %   symbols (cell(1, natoms) of charstring): atomic symbols
        %   SEAs (cell(1, \*) of double(\*, 1)): cell array of vectors of the indices of similar atoms
        %
        % Returns:
        %   horizontal (logical): whether there is a horizontal reflection plane
        %   nv (integer): number of vertical reflection planes
            horizontal = false;
            n = self.ord;
            vertical = zeros(3, n);
            nvertical = 0;
            for i = 1:length(SEAs)
                SEA = SEAs{i};
                k = length(SEA);
                if k > 1
                    for p = 2:k
                        for q = 1:p-1
                            normal = coords(:, SEA(p)) - coords(:, SEA(q));
                            ok = self.Reflection(normal, coords, symbols, self.ndec);
                            if ok
                                dotprod = round(dot(normal, self.principalAxis)/norm(normal)/norm(self.principalAxis), self.ndec);
                                if dotprod == 1
                                    horizontal = true;
                                elseif dotprod == 0
                                    if any(all(abs(vertical - normal) < 10^(-1*self.ndec)))
                                        continue
                                    end
                                    nvertical = nvertical + 1;
                                    vertical(:, nvertical) = normal;
                                end
                            end
                        end
                        if horizontal || nvertical == n
                            break
                        end
                    end
                end
            end
        end

        function determineSphericalGroup(self, coords, symbols, SEAs)
        % Determine spherical group of the molecule
        %
        % Determines whether molecule has tetrahedral, octahedral, or icosahedral symmetry
        % and assigns point group
        %
        % Args:
        %   coords (double(3, natoms)): atomic coordinates
        %   symbols (cell(1, natoms) of charstring): atomic symbols
        %   SEAs (cell(1, \*) of double(\*, 1)): cell array of vectors of the indices of similar atoms
            % Check number of C5 axes greater than 2
            c5 = self.findCnAxes(5, 2, coords, symbols, SEAs);
            if c5
                if self.coordsEqual(coords, coords*-1, symbols, self.ndec)
                    self.pointGroup = 'Ih';
                    return
                else
                    self.pointGroup = 'I';
                    return
                end
            end
            % Check number of C4 axes greater than 2
            c4 = self.findCnAxes(4, 2, coords, symbols, SEAs);
            if c4
                if self.coordsEqual(coords, coords*-1, symbols, self.ndec)
                    self.pointGroup = 'Oh';
                    return
                else
                    self.pointGroup = 'O';
                    return
                end
            end
            % Otherwise tetrahedral
            for i = 1:length(SEAs)
                SEA = SEAs{i};
                k = length(SEA);
                if k > 1
                    for p = 2:k
                        for q = 1:p-1
                            normal = coords(:, SEA(p)) - coords(:, SEA(q));
                            ok = self.Reflection(normal, coords, symbols, self.ndec);
                            if ok
                                if self.coordsEqual(coords, -1*coords, symbols, self.ndec)
                                    self.pointGroup = 'Th';
                                    return
                                else
                                    self.pointGroup = 'Td';
                                    return
                                end
                            end
                        end
                    end
                end
            end
            self.pointGroup = 'T';
        end

        function determineLowSymmetry(self, coords, symbols, SEAs)
        % Determines point group on molecules without Cn axis
        %
        % Args:
        %   coords (double(3, natoms)): atomic coordinates
        %   symbols (cell(1, natoms) of charstring): atomic symbols
        %   SEAs (cell(1, \*) of double(\*, 1)): cell array of vectors of the indices of similar atoms
            natoms = size(coords, 2);
            nplane = cross(coords(:, 1) - coords(:, 2), coords(:, 1) - coords(:, 3));
            mirror = true;
            for i = min(natoms, 4):natoms
                if round(dot(coords(:, i), nplane), self.ndec) ~= 0
                    mirror = false;
                    break
                end
            end
            if mirror
                self.pointGroup = 'Cs';
                return
            else
                for i = 1:length(SEAs)
                    SEA = SEAs{i};
                    k = length(SEA);
                    if k > 1
                        for p = 2:k
                            for q = 1:p-1
                                normal = coords(:, SEA(p)) - coords(:, SEA(q));
                                ok = self.Reflection(normal, coords, symbols, self.ndec);
                                if ok
                                    self.pointGroup = 'Cs';
                                    return
                                end
                            end
                        end
                    end
                end
            end
            if ~mirror
                if self.coordsEqual(coords, -1*coords, symbols, self.ndec)
                    self.pointGroup = 'Ci';
                else
                    self.pointGroup = 'C1';
                end
            end
        end

        function findPointGroup(self, coords, symbols, masses, SEAs)
        % Determine the point group of molecule by symmetry operation flowchart
        %
        % Algorithm from Beruski and flowchart from Vallance
        %
        % Args:
        %   coords (double(3, natoms)): atomic coordinates
        %   symbols (cell(1, natoms) of charstring): atomic symbols
        %   masses (double(1, natoms)): atomic masses
        %   SEAs (cell(1, \*) of double(\*, 1)): cell array of vectors of the indices of similar atoms
            [I, rcm] = self.InertiaTensor(coords', masses);
            evals = round(eig(I), self.ndec);
            neval = length(unique(evals));
            if any(iszero(evals)) && length(unique(evals(~iszero(evals)))) == 1
                % if I1 == 0 and I2 == I3 then the molecule is linear
                if self.coordsEqual(coords, -1*coords, symbols, self.ndec)
                    self.pointGroup = 'Dinfv';
                else
                    self.pointGroup = 'Cinfv';
                end
            elseif neval == 1
                % if I1 == I2 == I3 then molecule is spherical
                % (tetrahedral, octahedral, or icosahedral)
                % http://faculty.otterbein.edu/djohnston/sym/common/images/flowchart.pdf
                self.determineSphericalGroup(coords, symbols, SEAs)
            else
                % otherwise molecule is symmetric or assymmetric rotor
                % Determine axis of highest rotation
                Cnaxes = zeros(3, length(SEAs));
                ords = zeros(1, length(SEAs));
                for i = 1:length(SEAs)
                    if length(SEAs{i}) > 2
                        [ords(i), Cnaxes(:, i)] = self.findRotations(coords, symbols, masses, SEAs{i});
                    else
                        ords(i) = 1;
                    end
                end
                [n, i] = max(ords);
                if n == 1
                    self.principalAxis = [0,0,0];
                else
                    self.principalAxis = Cnaxes(:, i);
                end
                if max(ords) == 1
                    % we need to check for axes of order 2
                    [ok, axis] = self.findC2Axes(coords, symbols, SEAs, rcm);
                    if ok
                        n = 2;
                        self.principalAxis = axis;
                    end
                end
                self.ord = n;
                if n == 1
                % molecule doesn't have a Cn axis so check for inversion and mirror planes
                    self.determineLowSymmetry(coords, symbols, SEAs)
                else
                % Check for C2 axes perpendicular to principal axis
                    [ind, perpAxes] = self.findPerpendicularAxes(Cnaxes, ords, rcm, SEAs, coords, symbols);
                    if ind > n
                    % dihedral group
                        % check for horizontal mirror plane (normal is principal axis)
                        ok = self.Reflection(self.principalAxis, coords, symbols, self.ndec);
                        if ok
                            self.pointGroup = sprintf('D%dh', n);
                        else
                        % check for dihedral mirror planes
                            nsigmad = self.findDihedralMirrorPlanes(coords, symbols, perpAxes);
                            if nsigmad == n
                                self.pointGroup = sprintf('D%dd', n);
                            else
                                self.pointGroup = sprintf('D%d', n);
                            end
                        end
                    else
                    % Cn, Cnh, Cnv, or Sn
                        % Check for planar molecules
                        [planar, nplane] = self.checkPlanar(coords);
                        nv = 0;
                        if planar
                            dotprod = round(abs(dot(nplane, self.principalAxis)/norm(nplane)/norm(self.principalAxis)), self.ndec);
                            if dotprod == 1
                                % horizontal mirror plane
                                self.pointGroup = sprintf('C%dh', n);
                                return
                            elseif dotprod == 0
                                % another vertical mirror plane
                                nv = 1;
                            end
                        end
                        % Find planes of reflection
                        [horizontal, nvertical] = self.findReflections(coords, symbols, SEAs);
                        nv = nv + nvertical;
                        if horizontal
                            self.pointGroup = sprintf('C%dh', n);
                        elseif nv >= n
                            self.pointGroup = sprintf('C%dv', n);
                        else
                        % Check for S2n axis collinear to Cn axis
                            ok = self.ImproperRotation(2*n, self.principalAxis, coords, symbols, self.ndec);
                            if ok
                                self.pointGroup = sprintf('S%d', 2*n);
                            else
                                self.pointGroup = sprintf('C%d', n);
                            end
                        end

                    end
                end
            end
        end

    end
end
