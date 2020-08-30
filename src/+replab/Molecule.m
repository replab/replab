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
        name % (charstring): name of `+replab.Molecule` object
        coords % (double(3, natoms)): the coordinates of atoms in the molecule
               %                                   coords are sorted by increasing values  
        symbols % (cell(1, natoms) of charstring): symbols of atoms corresponding to coordinates
        masses % (double(1, natoms)): masses of atoms corresponding to coordinates
        natoms % (integer): total number of atoms
        id % (`+replab.IdentifyPointGroup`): identification of molecular point group                                                                        through an axis on coordinates
    end
    
    methods (Static)
        
        function mol = fromFile(filename, atomicSymbols, atomicMasses)
        % Creates a molecule object from a file
        %
        % Args:
        %   filename (charstring): name of file with molecular information
        %   atomicSymbols (optional cell(1, \*) of charstring): atomic symbols of masses to 
        %                                                       replace or add to stored masses
        %   atomicMasses (optional cell(1, \*) of double): atomic masses corresponding to atomicSymbols
        %
        % Returns:
        %   mol (`+replab.Molecule`): molecule object
            [coords, sym] = replab.Molecule.readXYZ(filename);
            ndec = 2;
            mol = replab.Molecule(filename, coords, sym);
            if nargin > 1
                masses = replab.Molecule.getMasses(mol.symbols, atomicSymbols, atomicMasses);
            else
                masses = replab.Molecule.getMasses(mol.symbols);
            end
            mol.masses = masses;
            mol.id = replab.IdentifyPointGroup(mol.coords, mol.symbols, mol.masses, ndec);
%             mol.id.table = mol.id.findCharacterTable(mol.id.pointGroup);
        end
        
        function mol = fromAtoms(coords, symbols, atomicSymbols, atomicMasses)
        % Creates a molecule object atomic information
        %
        % Args:
        %   coords (double(3, natoms)): atomic coordinates (can be unsorted)
        %   symbols (cell(1, natoms) of charstring): atomic symbols
        %   atomicSymbols (optional cell(1, \*) of charstring): atomic symbols of masses to 
        %                                                       replace or add to stored masses
        %   atomicMasses (optional cell(1, \*) of double): atomic masses corresponding to atomicSymbols
        %
        % Returns:
        %   mol (`+replab.Molecule`): molecule object
            ndec = 2;
            mol = replab.Molecule(filename, coords, symbols);
            if nargin > 1
                masses = replab.Molecule.getMasses(mol.symbols, atomicSymbols, atomicMasses);
            else 
                masses = replab.Molecule.getMasses(mol.symbols);
            end
            mol.masses = masses;
            mol.id = replab.IdentifyPointGroup(mol.coords, mol.symbols, mol.masses, ndec);
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
        
        function masses = getMasses(symbols, atomicSymbols, atomicMasses)
        % Get the atomic masses for the atom represented by each symbol
        %
        % Args:
        %   symbols (cell(1, natoms) of charstring): symbols of atoms
        %   atomicSymbols (optional cell(1, \*) of charstring): atomic symbols of masses to 
        %                                                       replace or add to stored masses
        %   atomicMasses (optional cell(1, \*) of double): atomic masses corresponding to atomicSymbols
        %
        % Returns
        %   masses (cell(1, natoms) of double): masses of atoms
            storedMasses = struct('Xe', 131.293, 'F', 18.998403, 'N', 14.003074, ...
                                  'H', 1.007825, 'C', 12.000000, 'O', 15.999, 'S', 32.065, ...
                                  'Cl', 35.453, 'B', 10.811, 'Br', 79.904);
            if nargin > 1
                for i = 1:length(atomicSymbols)
                    storedMasses.(atomicSymbols{i}) = atomicMasses{i};
                end
            end
            masses = cellfun(@(s) storedMasses.(s), symbols);
        end
        
    end
    
    methods
        
        function self = Molecule(name, coords, symbols)
        % Creates a molecule object
        % 
        % Convention: automatically sort coordinates in ascending order
        % rounding to 2 decimal places for equality
        %
        % Args:
        %   coords (double(3, natoms)): a set of coordinates
        %   symbols (cell(1, natoms) of charstring): atomic symbols of coordinates
        %   ndec (integer): number of decimal places to round when checking for equality
            self.name = name;
            [coords, inds] = replab.IdentifyPointGroup.sortCoords(coords, 2);
            symbols = symbols(inds);
            self.coords = coords;
            self.symbols = symbols;
            self.natoms = size(coords, 2);
        end
        
        function action = RotationAction(self, k, axis, ndec)
        % Finds action of a rotation by 2pi/k on the molecule
        %
        % Args:
        %   k (integer): rotation by angle 2*pi/k
        %   axis (double(1, 3)): axis about which we rotate molecule
        %   ndec (integer) (optional): the number of decimals to round for equality
            if nargin < 4
                if ~isempty(self.id)
                    ndec = self.id.ndec;
                else
                    ndec = 2;
                end
            end
            [~, action] = replab.IdentifyPointGroup.Rotation(k, axis, self.coords, self.symbols, ndec);
        end
        
        function action = ReflectionAction(self, normal, ndec)
        % Finds action of a reflection by 2pi/k on the molecule
        %
        % Args:
        %   normal (double(1, 3)): normal of reflection plane
        %   ndec (integer) (optional): the number of decimals to round for equality
            if nargin < 4
                if ~isempty(self.id)
                    ndec = self.id.ndec;
                else
                    ndec = 2;
                end
            end
            [~, action] = replab.IdentifyPointGroup.Reflection(normal, self.coords, self.symbols, ndec);
        end
        
        function action = ImproperRotationAction(self, k, axis, ndec)
        % Finds action of a rotation by 2pi/k then a reflection on the molecule
        %
        % Args:
        %   k (integer): rotation by angle 2*pi/k
        %   axis (double(1, 3)): rotation axis and normal of reflection plane
        %   ndec (integer) (optional): the number of decimals to round for equality
            if nargin < 4
                if ~isempty(self.id)
                    ndec = self.id.ndec;
                else
                    ndec = 2;
                end
            end
            [~, action] = replab.IdentifyPointGroup.ImproperRotation(k, axis, self.coords, self.symbols, ndec);
        end
        
    end
end

