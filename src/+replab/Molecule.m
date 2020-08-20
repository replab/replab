classdef Molecule < replab.Str
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
        coords % (double(3, natoms)): the coordinates of atoms in the molecule
               %                                   coords are sorted by increasing values  
        symbols % (cell(1, natoms) of charstring): symbols of atoms corresponding to coordinates
        masses % (double(1, natoms)): masses of atoms corresponding to coordinates
        natoms % (integer): total number of atoms
        SEAs % (cell(1, \*) of double(\*, 1)): cell array of vectors of the indices of similar atoms
        principalAxis % (double(3, 1)): principal axis of rotation of the molecule
        ord % (integer): highest order rotation of the molecule that preserves symmetry
        pointGroup % (charstring): characters representing the point group
        table % (`+replab.CharacterTable`): the character table for the molecule's point group
        ndec % (integer): number of decimal places to check for equality
    end
    
    methods (Static)
        
        function mol = createMolecule(filename)
        % Creates a molecule object from a file
        %
        % Args:
        %   filename (charstring): name of file with molecular information
        %
        % Returns:
        %   mol (`+replab.Molecule`): molecule object
            [coords, sym] = replab.Molecule.readXYZ(filename);
            [coords, sym] = replab.Molecule.sortCoords(coords, sym, 2);
            mol = replab.Molecule(coords, sym, 2);
            mol.masses = replab.Molecule.getMasses(mol.symbols);
            mol.SEAs = replab.Molecule.groupAtoms(mol.distanceMatrix, mol.ndec);
            mol.findPointGroup
            mol.table = mol.findCharacterTable(mol.pointGroup);
        end
        
        function [coords, symbols] = readXYZ(filename)
         % Reads the coordinates and atomic symbols from an xyz file
         %
         % Args:
         %   filename (charstring): name and location of xyz file
         %
         % Returns
         %   coords (double(3, natoms)): coordinates from xyz file
         %   symbols (cell(1, natoms) of charstring): symbols of atoms from xyz file
            fileID = fopen(filename, 'r');
            [A, ~] = fscanf(fileID,'%c');
            atoms = split(A, newline);
            natoms = str2double(atoms{1});
            atoms = atoms(3:end);
            atoms(cellfun(@isempty, atoms)) = [];
            coords = zeros(natoms, 3);
            symbols = cell(1, natoms);
            for i = 1:natoms
                str = split(atoms{i});
                symbols{i} = str{1};
                coords(i,:) = cellfun(@str2double, str(2:end));
            end
            coords = coords';
        end
        
        function ct = findCharacterTable(pointGroup)
        % Determines the character table of a molecule with a certain point group
            pg = split(pointGroup, '');
            pg = pg(2:end-1);
            if pg{1} == 'T'
                if length(pg) > 1
                    if pg{2} == 'h'
                        % direct product C2 and alternating group 4
                    elseif pg{2} == 'd'
                        % symmetric group 4
                        ct = replab.ct.PermutationCharacterTable(replab.S(4));
                    end
                else
                    % alternating group 4
                end
            elseif pg{1} == 'O'
                % symmetric group 4
                if length(pg) > 1 && pg{2} == 'h'
                    % direct product C2 and S4
                else
                    ct = replab.ct.PermutationCharacterTable(replab.S(4));
                end
            elseif pg{1} == 'I'
                % alternating group 5
                if length(pg) > 1 && pg{2} == 'h'
                    % direct product C2 and alternating group 5 
                else
                    % alternating group 5
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
        
        function masses = getMasses(symbols, atomicMasses)
        % Get the atomic masses for the atom represented by each symbol
        %
        % Args:
        %   symbols (cell(1, natoms) of charstring): symbols of atoms
        %   atomicMasses (optional struct): masses of symbols given
        %
        % Returns
        %   masses (cell(1, natoms) of double): masses of atoms
            if nargin < 2
                atomicMasses = struct('Xe', 131.293, 'F', 18.998403, 'N', 14.003074, ...
                                      'H', 1.007825, 'C', 12.000000, 'O', 15.999, 'S', 32.065, ...
                                      'Cl', 35.453, 'B', 10.811, 'Br', 79.904);
            end
            masses = cellfun(@(s) atomicMasses.(s), symbols);
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
        
        function new_coords = Rotation(k, axis, coords, ndec)
        % Apply rotation of 2*pi/k to coords about an axis
        %
        % Args: 
        %   k (integer): number of rotations which should keep coordinates the same
        %   axis (double(1, 3)): axis about which to rotate coordinates
        %   coords (double(3, natoms)): original coordinates
        %   ndec (integer): number of decimals to round
        %
        % Returns:
        %   new_coords (double(3, natoms)): rotated coordinates
            theta = 2*pi/k;
            axis(round(axis, ndec) == 0) = 0;
            a = axis / norm(axis);
            R = [cos(theta)+a(1)^2*(1-cos(theta)), a(1)*a(2)*(1-cos(theta))-a(3)*sin(theta), a(1)*a(3)*(1-cos(theta))+a(2)*sin(theta);
                 a(1)*a(2)*(1-cos(theta))+a(3)*sin(theta), cos(theta)+a(2)^2*(1-cos(theta)), a(2)*a(3)*(1-cos(theta))-a(1)*sin(theta);
                 a(1)*a(3)*(1-cos(theta))-a(2)*sin(theta), a(2)*a(3)*(1-cos(theta))+a(1)*sin(theta), cos(theta)+a(3)^2*(1-cos(theta))];
            new_coords = R * coords;
        end

        function new_coords = Reflection(normal, coords, ndec)
        % Reflect coords through a plane with a given normal
        %
        % Args:
        %   normal (double(1, 3)): normal of plane through which to reflect
        %   coords (double(3, natoms)): original coordinates
        %   ndec (integer): number of decimals to round
        %
        % Returns:
        %   new_coords (double(4, natoms)): reflected coordinates
            normal(round(normal, ndec) == 0) = 0;
            normal = normal / norm(normal); % normalize normal vector
            sigma = [1-2*normal(1)^2, -2*normal(1)*normal(2), -2*normal(1)*normal(3);
                     -2*normal(1)*normal(2), 1-2*normal(2)^2, -2*normal(2)*normal(3);
                     -2*normal(1)*normal(3), -2*normal(2)*normal(3), 1-2*normal(3)^2];
            new_coords = sigma * coords;
        end
        
        function [coords, symbols] = sortCoords(coords, symbols, ndec)
        % Sorts coordinates in ascending order
        %
        % Convention: atoms are sorted by increasing x-coordinate, then
        %             y-coordinate, and then z-coordinate
        %
        % Args:
        %   coords (double(3, natoms)): set of coordinates
        %   symbols (cell(1, natoms) of charstring): atomic symbols
        %   ndec (integer): number of decimals to check equality
        %
        % Returns:
        %   coords (double(3, natoms)): sorted coordinates
        %   symbols (cell(1, natoms) of charstring): symbols sorted in same order as coords
            natoms = size(coords, 2);
            [~, i1] = sort(coords(1, :));
            coords = coords(:, i1);
            symbols = symbols(i1);
            [~, v1] = unique(round(coords(1, :), ndec));
            degen1 = 1:natoms;
            degen1(v1) = [];
            while ~isempty(degen1)
                f1 = find(round(coords(1, :), ndec) == round(coords(1, degen1(1)), ndec));
                degen1(1:length(f1)-1) = [];
                [~, i2] = sort(coords(2, f1));
                coords(:, f1) = coords(:, f1(i2));
                symbols(:, f1) = symbols(f1(i2));

                [~, v2] = unique(round(coords(2, f1), ndec));
                degen2 = 1:length(f1);
                degen2(:, v2) = [];
                while ~isempty(degen2)
                    f2 = find(round(coords(2, f1), ndec) == round(coords(2, f1(degen2(1))), ndec));
                    degen2(1:length(f2)-1) = [];
                    [~, i3] = sort(coords(3, f1(f2)));
                    coords(:, f1(f2)) = coords(:, f1(f2(i3)));
                    symbols(:, f1(f2)) = symbols(f1(f2(i3)));
                end
            end
        end
        
    end
    
    methods
        
        function self = Molecule(coords, symbols, ndec)
        % Creates a molecule object
        %
        % Args:
        %   coords (double(3, natoms)): a set of coordinates
        %   symbols (cell(1, natoms) of charstring): atomic symbols of coordinates
        %   ndec (integer): number of decimal places to round when checking for equality
            self.coords = coords;
            self.symbols = symbols;
            self.natoms = size(coords, 2);
            self.ndec = ndec;
        end
        
        function D = distanceMatrix(self)
        % Determine the distance matrix of the molecules
        %
        % Convention: entries of row i are the distances of atom i from
        %             atoms 1, 2, ..., natoms so D(i, i) = 0
        %
        % Returns:
        %   D (double(natoms, natoms)): distance matrix
            D = zeros(self.natoms);
            for i = 1:self.natoms
                for j = 1:self.natoms
                    D(i,j) = sqrt(sum((self.coords(:,i)-self.coords(:,j)).^2));
                end
            end
        end
        
        function new_coords = SnAxis(self, k, axis)
        % Rotate and reflect plane with the same axis
        %
        % Args:
        %   k (integer): rotation order (2*pi/k)
        %   axis (double(1, 3)): rotate about this axis and reflect through plane 
        %                        with axis as the normal
        %
        % Returns:
        %   new_coords (double(3, natoms)): rotated and reflected coordinates
            rot_coords = self.Rotation(k, axis, self.coords, self.ndec, self.ndec);
            new_coords = self.Reflection(axis, rot_coords, self.ndec);
        end
        
        function ok = coordsEqual(self, new_coords)
        % Checks whether a new set of coordinates is equivalent to original ones
        %
        % Convention: assumes self.coords is already sorted
        %
        % Args:
        %   new_coords (double(3, natoms)): set of coordinates to check for equality
        %
        % Returns:
        %   ok (logical): whether new coordinates are equal to self.coords
            [new_coords, new_symbols] = self.sortCoords(new_coords, self.symbols, self.ndec);
            ok = isequal(round(self.coords, self.ndec), round(new_coords, self.ndec));
            if ok
                ok = isequal(self.symbols, new_symbols);
            end
        end
        
        function [order, axis] = findRotations(self, inds)
        % Determine the rotation order and axis for a group of atoms
        %
        % Args:
        %   inds (integer(1, \*)): indices of atoms
        %
        % Returns:
        %   order (integer): highest order rotation possible
        %   axis (double(1, 3)): axis of rotation
            k = length(inds);
            if k < 3
                return
            end
            [I_atom, ~] = self.InertiaTensor(self.coords(:, inds)', self.masses(inds));
            [evecs, evals] = eig(I_atom);
            moments = sort(round(diag(evals), self.ndec));
            u = unique(moments);
            if abs(sum(moments(1:2)) - moments(3)) < self.ndec
                if length(u) == 2
                    % Ck rotations
                    n = k;
                else 
                    n = floor(k/2);
                    while rem(k, n) ~= 0
                        n = n - 1;
                    end
                end
                evec = 3;
                new_coords = self.Rotation(n, evecs(:, evec), self.coords, self.ndec);
                same = self.coordsEqual(new_coords);
                if ~same && length(u) ~= 2
                    evec = 2;
                    new_coords = self.Rotation(n, evecs(:, evec), self.coords, self.ndec);
                    same = self.coordsEqual(new_coords);
                    if ~same
                        evec = 1;
                        new_coords = self.Rotation(n, evecs(:, evec), self.coords, self.ndec);
                        same = self.coordsEqual(new_coords);
                    end
                    if ~same
                        evec = 3;
                    end
                end
                while ~same
                    n = min(n - 1, floor(k/2));
                    while rem(k, n) ~= 0
                        n = n - 1;
                    end
                    new_coords = self.Rotation(n, evecs(:, evec), self.coords, self.ndec);
                    same = self.coordsEqual(new_coords);
                end
                order = n;
                axis = evecs(:, evec);
            else
                if length(u) == 3
                    order = 2;
                    i = 1;
                    new_coords = self.Rotation(2, evecs(:, i), self.coords, self.ndec);
                    same = self.coordsEqual(new_coords);
                    while ~same
                        i = 1 + 1;
                        new_coords = self.Rotation(2, evecs(:, i), self.coords, self.ndec);
                        same = self.coordsEqual(new_coords);
                    end
                    axis = evecs(:, i);
                else
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
                    new_coords = self.Rotation(n, evecs(:, naxis), self.coords, self.ndec);
                    same = self.coordsEqual(new_coords);
                    while ~same
                        n = min(n - 1, floor(k/2));
                        while rem(k, n) ~= 0
                            n = n - 1;
                        end
                        new_coords = self.Rotation(n, evecs(:, naxis), self.coords, self.ndec);
                        same = self.coordsEqual(new_coords);
                    end
                    order = n;
                    axis = evecs(:, naxis);
                end
            end
        end
        
        function ok = findSphericalArrangement(self, n, max_rotations)
        % Find whether there are multiple higher order rotation axes
        %
        % Returns:
        %   n (integer): order of rotation
        %   max_rotations (integer): check whether there are at least this many 
        %                            rotation axes of order n
        %
        % Returns
        %   ok (logical): whether there are at least max_rotations rotation axes of order n
            nrot = 0;
            ok = false;
            nSEA = cellfun(@length, self.SEAs);
            f = find(nSEA > n);
            for i = 1:length(f)
                locs = self.SEAs{f(i)};
                k = length(locs);
                for p = 2:k
                    for q = 1:p-1
                        ri = self.coords(:, locs(p));
                        rj = self.coords(:, locs(q));
                        rij = ri - rj;
                        for m = 1:q-1
                            rk = self.coords(:, locs(m));
                            rik = ri - rk;
                            if round(norm(rij), self.ndec) ~= round(norm(rik), self.ndec)
                                continue
                            end
                            new_coords = self.Rotation(n, cross(rij, rik), self.coords, self.ndec);
                            if self.coordsEqual(new_coords)
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
        
        function findPointGroup(self)
        % Determine the point group of molecule by symmetry operation flowchart
        %
        % Algorithm from Beruski and flowchart from Vallance
            [I, rcm] = self.InertiaTensor(self.coords', self.masses);
            evals = round(eig(I), self.ndec);
            neval = length(unique(evals));
            if any(iszero(evals)) && length(unique(evals(~iszero(evals)))) == 1
                % linear
                if self.coordsEqual(-1*self.coords)
                    self.pointGroup = 'Dinfv';
                else
                    self.pointGroup = 'Cinfv';
                end
            elseif neval == 1
                % spherical
                % http://faculty.otterbein.edu/djohnston/sym/common/images/flowchart.pdf
                % Check number of C5 axes
                c5 = self.findSphericalArrangement(5, 2);
                if c5
                    if self.coordsEqual(self.coords*-1)
                        self.pointGroup = 'Ih';
                        return
                    else
                        self.pointGroup = 'I';
                        return
                    end
                end
                % Check number of C4 axes
                c4 = self.findSphericalArrangement(4, 2);
                if c4
                    if self.coordsEqual(self.coords*-1)
                        self.pointGroup = 'Oh';
                        return
                    else
                        self.pointGroup = 'O';
                        return
                    end
                end
                mirror = false;
                for i = 1:length(self.SEAs)
                    SEA = self.SEAs{i};
                    k = length(SEA);
                    if k > 1
                        for p = 2:k
                            for q = 1:p-1
                                normal = self.coords(:, SEA(p)) - self.coords(:, SEA(q));
                                new_coords = self.Reflection(normal, self.coords, self.ndec);
                                if self.coordsEqual(new_coords)
                                    mirror = true;
                                    break
                                end
                            end
                            if mirror
                                break
                            end
                        end
                    end
                    if mirror
                        break
                    end
                end
                if mirror
                    if self.coordsEqual(-1*self.coords)
                        self.pointGroup = 'Th';
                    else
                        self.pointGroup = 'Td';
                    end
                else
                    self.pointGroup = 'T';
                end
            else
                % symmetric or asymmetric
                Cnaxes = zeros(3, length(self.SEAs));
                ords = zeros(1, length(self.SEAs));
                for i = 1:length(self.SEAs)
                    if length(self.SEAs{i}) > 2
                        [ords(i), Cnaxes(:, i)] = self.findRotations(self.SEAs{i});
                    else
                        ords(i) = 1;
                    end
                end
                [n, i] = max(ords);
                self.ord = n;
                if n == 1
                    self.principalAxis = [0,0,0]; 
                else
                    self.principalAxis = Cnaxes(:, i);
                end
                if max(ords) == 1
                    % check for axes of order 2
                    nSEA = cellfun(@length, self.SEAs);
                    diatomic = find(nSEA == 2);
                    if length(diatomic) > 1
                        rij = self.coords(:, self.SEAs{diatomic(1)}(1)) - self.coords(:, self.SEAs{diatomic(1)}(2));
                        rkl = self.coords(:, self.SEAs{diatomic(2)}(1)) - self.coords(:, self.SEAs{diatomic(2)}(2));
                        rC2 = cross(rij, rkl);
                        new_coords = self.Rotation(2, rC2, self.coords, self.ndec);
                        if self.coordsEqual(new_coords)
                            n = 2;
                            self.principalAxis = rC2;
                        end 
                    else
                        for i = 1:length(self.SEAs)
                            SEA = self.SEAs{i};
                            ok = false;
                            if length(SEA) > 2
                                continue
                            end
                            if length(SEA) == 2
                                mid = (self.coords(:, SEA(1)) + self.coords(:, SEA(2))) / 2;
                                if ~isequal(round(mid, self.ndec), round(rcm', self.ndec))
                                    rC2 = mid - rcm';
                                    new_coords = self.Rotation(2, rC2, self.coords, self.ndec);
                                    if self.coordsEqual(new_coords)
                                        ok = true;
                                    end
                                end
                                if ~ok
                                    rC2 = self.coords(:, SEA(2)) - rcm';
                                    new_coords = self.Rotation(2, rC2, self.coords, self.ndec);
                                    if self.coordsEqual(new_coords)
                                        ok = true;
                                    end
                                end
                            end
                            if ~ok  
                                rC2 = self.coords(:, SEA(1)) - rcm';
                                new_coords = self.Rotation(2, rC2, self.coords, self.ndec);
                                if self.coordsEqual(new_coords)
                                    ok = true;
                                end
                            end

                            if ok
                                n = 2;
                                self.principalAxis = rC2;
                                break
                            end
                        end
                    end
                end
                if n == 1 
                % molecule doesn't have a Cn axis so check for inversion and mirror planes
                    nplane = cross(self.coords(:, 1) - self.coords(:, 2), self.coords(:, 1) - self.coords(:, 3));
                    mirror = true;
                    for i = min(self.natoms, 4):self.natoms
                        if round(dot(self.coords(:, i), nplane), self.ndec) ~= 0
                            mirror = false;
                            break
                        end
                    end
                    if mirror
                        self.pointGroup = 'Cs';
                        return
                    else
                        for i = 1:length(self.SEAs)
                            SEA = self.SEAs{i};
                            k = length(SEA);
                            if k > 1
                                for p = 2:k
                                    for q = 1:p-1
                                        normal = self.coords(:, SEA(p)) - self.coords(:, SEA(q));
                                        new_coords = self.Reflection(normal, self.coords, self.ndec);
                                        if self.coordsEqual(new_coords)
                                            self.pointGroup = 'Cs';
                                            return
                                        end
                                    end
                                end
                            end
                        end
                    end
                    if ~mirror
                        if self.coordsEqual(-1*self.coords)
                            self.pointGroup = 'Ci';
                        else
                            self.pointGroup = 'C1';
                        end
                    end
                else
                % Check for C2 axes perpendicular to principal axis
                    perpAxes = zeros(3, n);
                    ind = 1;
                    for m = 1:length(ords)
                        if ords(m) == 2 && isequal(round(dot(self.principalAxis, Cnaxes(:,m)), self.ndec), 0)
                            perpAxes(:, ind) = Cnaxes(:, m);
                            ind = ind + 1;
                        end
                    end
                    for i = 1:length(self.SEAs)
                        SEA = self.SEAs{i};
                        k = length(SEA);
                        for p = 2:k
                            for q = 1:p-1
                                midpoint = (self.coords(:, SEA(p)) + self.coords(:, SEA(q))) / 2;
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
                                rot_coords = self.Rotation(2, rperp, self.coords, self.ndec);
                                if self.coordsEqual(rot_coords)
                                    perpAxes(:, ind) = rperp;
                                    ind = ind + 1;
                                end
                            end
                        end
                    end
                    if ind > n
                    % dihedral group
                        % check for horizontal mirror plane (normal is principal axis)
                        reflect_coords = self.Reflection(self.principalAxis, self.coords, self.ndec);
                        if self.coordsEqual(reflect_coords)
                            self.pointGroup = sprintf('D%dh', n);
                        else
                        % check for dihedral mirror planes
                            nsigmad = 0;
                            for p = 2:n
                                for q = 1:p-1
                                    mid = (perpAxes(:, p) + perpAxes(:, q)) / 2;
                                    sigmad = cross(self.principalAxis, mid);
                                    new_coords = self.Reflection(sigmad, self.coords, self.ndec);
                                    if self.coordsEqual(new_coords)
                                        nsigmad = nsigmad + 1;
                                    end
                                    % We get a second dihedral mirror plane from the same pair 
                                    % of axes (checking this may be redundant)
                                    mid = (perpAxes(:, p) + -1*perpAxes(:, q)) / 2;
                                    sigmad = cross(self.principalAxis, mid);
                                    new_coords = self.Reflection(sigmad, self.coords, self.ndec);
                                    if self.coordsEqual(new_coords)
                                        nsigmad = nsigmad + 1;
                                    end
                                    if nsigmad == n
                                        break
                                    end
                                end
                            end
                            if nsigmad == n
                                self.pointGroup = sprintf('D%dd', n);
                            else
                                self.pointGroup = sprintf('D%d', n);
                            end
                        end
                    else
                    % Cn, Cnh, Cnv, or Sn
                        % Another reflection in planar molecules
                        nv = 0;
                        nplane = cross(self.coords(:, 1) - self.coords(:, 2), self.coords(:, 1) - self.coords(:, 3));
                        planar = true;
                        for i = min(self.natoms, 4):self.natoms
                            if round(dot(self.coords(:, i), nplane), self.ndec) ~= 0
                                planar = false;
                                break
                            end
                        end
                        if planar
                            dotprod = round(abs(dot(nplane, self.principalAxis)/norm(nplane)/norm(self.principalAxis)), self.ndec);
                            if dotprod == 1
                                self.pointGroup = sprintf('C%dh', n);
                                return
                            elseif dotprod == 0
                                nv = 1;
                            end
                        end
                        % Find planes of reflection
                        horizontal = false;
                        vertical = zeros(3, n);
                        for i = 1:length(self.SEAs)
                            SEA = self.SEAs{i};
                            k = length(SEA);
                            if k > 1
                                for p = 2:k
                                    for q = 1:p-1
                                        normal = self.coords(:, SEA(p)) - self.coords(:, SEA(q));
                                        new_coords = self.Reflection(normal, self.coords, self.ndec);
                                        if self.coordsEqual(new_coords)
                                            dotprod = round(dot(normal, self.principalAxis)/norm(normal)/norm(self.principalAxis), self.ndec);
                                            if dotprod == 1
                                                horizontal = true;
                                            elseif dotprod == 0
                                                if any(all(abs(vertical - normal) < 10^(-1*self.ndec)))
                                                    continue
                                                end
                                                nv = nv + 1;
                                                vertical(:, nv) = normal;
                                            end
                                        end
                                    end
                                    if horizontal || nv == n
                                        break
                                    end
                                end
                            end
                        end
                        if horizontal
                            self.pointGroup = sprintf('C%dh', n);
                        elseif nv >= n
                            self.pointGroup = sprintf('C%dv', n);
                        else
                        % Check for S2n axis collinear to Cn axis
                            rot_coords = self.Rotation(2*n, self.principalAxis, self.coords, self.ndec);
                            new_coords = self.Reflection(self.principalAxis, rot_coords, self.ndec);
                            if self.coordsEqual(new_coords)
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

