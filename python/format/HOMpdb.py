import __init__
import molecule.HOMmolecule


class PDBformat:
	pass

	# PDB class
class PDB:

    def __init__ (self):
        self.AllAtoms={}

    # Loads a PDB from a file
    # Requires: FileName, a string containing the filename
    def LoadPDB(self, FileName):

        autoindex = 1

        self.__init__()

        # Now load the file into a list
        file = open(FileName,"r")
        lines = file.readlines()
        file.close()

        for t in range(0,len(lines)):
            line=lines[t]
            if len(line) >= 7:
                if line[0:4]=="ATOM" or line[0:6]=="HETATM": # Load atom data (coordinates, etc.)
                    TempAtom = atom()
                    TempAtom.ReadPDBLine(line)

                    self.AllAtoms[autoindex] = TempAtom # because points files have no indecies
                    autoindex = autoindex + 1

    # Saves the PDB object to a PDB file
    # Requires: filename to be saved
    def SavePDB(self, filename):

        if len(self.AllAtoms) > 0: # so the pdb is not empty (if it is empty, don't save)

            file = open(filename,"w")

            # write coordinates
            for atomindex in self.AllAtoms:
                file.write(self.AllAtoms[atomindex].CreatePDBLine() + "\n")

            file.close()

    # identifies the greatest distance between any two atoms of the PDB
    # Used for quickly comparing two PDBs so RMSD alignment not necessary if they're different
    # Returns float = distance
    def farthest_away_atoms_dist(self):
        dist_max = 0.0
        keys = self.AllAtoms.keys()
        for key_index1 in range(0,len(keys)-1):
            key1 = keys[key_index1]
            atom1 = self.AllAtoms[key1]
            for key_index2 in range(key_index1+1,len(keys)):
                key2 = keys[key_index2]
                atom2 = self.AllAtoms[key2]
                dist = atom1.coordinates.dist_to(atom2.coordinates)
                if dist > dist_max: dist_max = dist
        self.max_inter_atom_distance = dist_max
        return dist_max

    # identifies the smallest distance between any two atoms of the PDB
    # Used for quickly comparing two PDBs so RMSD alignment not necessary if they're different
    # Returns float = distance
    def closest_atoms_dist(self):
        dist_min = 999999999.99
        keys = self.AllAtoms.keys()
        for key_index1 in range(0,len(keys)-1):
            key1 = keys[key_index1]
            atom1 = self.AllAtoms[key1]
            for key_index2 in range(key_index1+1,len(keys)):
                key2 = keys[key_index2]
                atom2 = self.AllAtoms[key2]
                dist = atom1.coordinates.dist_to(atom2.coordinates)
                if dist < dist_min: dist_min = dist
        self.min_inter_atom_distance = dist_min
        return dist_min

    # Creates a "distance fingerprint" so PDB's can be quickly compared w/o RMSD alignment
    # Returns: a list containing sorted distances
    def distance_fingerprint(self):
        fingerprint = []
        keys = self.AllAtoms.keys()
        for key_index1 in range(0,len(keys)-1):
            key1 = keys[key_index1]
            atom1 = self.AllAtoms[key1]
            for key_index2 in range(key_index1+1,len(keys)):
                key2 = keys[key_index2]
                atom2 = self.AllAtoms[key2]
                dist = atom1.coordinates.dist_to(atom2.coordinates)
                fingerprint.append(dist)
        fingerprint.sort()
        self.distance_fingerprint = fingerprint
        return fingerprint

    # To detect if two distance fingerprints are sufficiently similar
    # Requires: Two fingerprints (lists) and a tolerance (float)
    # Returns: True or False
    def fingerprint_same_as(self, fingerprint, tol):
        if len(self.distance_fingerprint) <> len(fingerprint): return False

        for index in range(len(self.distance_fingerprint)):
                item1 = self.distance_fingerprint[index]
                item2 = fingerprint[index]
                if math.fabs(item1-item2) > tol: return False
        return True

    # Print out info about the PDB
    def print_out_info(self):
        for index in self.AllAtoms:
            print self.AllAtoms[index].CreatePDBLine()

    # Translate the molecule
    def TranslateMolecule(self, x, y, z):
        for atomindex in self.AllAtoms:
            self.AllAtoms[atomindex].coordinates.x = self.AllAtoms[atomindex].coordinates.x + x
            self.AllAtoms[atomindex].coordinates.y = self.AllAtoms[atomindex].coordinates.y + y
            self.AllAtoms[atomindex].coordinates.z = self.AllAtoms[atomindex].coordinates.z + z

    # Rotate the PDB around it's center.
    # Requires: degrees to rotate, as float
    def RotatePDB(self, thetax, thetay, thetaz):

        # first, identify the geometric center
        x = 0.0
        y = 0.0
        z = 0.0
        count = 0
        for index in self.AllAtoms:
                if self.AllAtoms[index].element != "H":
                        count = count + 1
                        x = x + self.AllAtoms[index].coordinates.x
                        y = y + self.AllAtoms[index].coordinates.y
                        z = z + self.AllAtoms[index].coordinates.z
        x = x / count
        y = y / count
        z = z / count

        # now, move the pdb to the origin
        self.TranslateMolecule(-x, -y, -z)

        # now rotate
        sinx = math.sin(thetax)
        siny = math.sin(thetay)
        sinz = math.sin(thetaz)
        cosx = math.cos(thetax)
        cosy = math.cos(thetay)
        cosz = math.cos(thetaz)

        cosy_cosz = cosy * cosz
        sinx_siny_cosz_plus_cosx_sinz = sinx * siny * cosz + cosx * sinz
        sinx_sinz_minus_cosx_siny_cosz = sinx * sinz - cosx * siny * cosz
        cosy_sinz = cosy * sinz
        cosx_cosz_minus_sinx_siny_sinz = cosx * cosz - sinx * siny * sinz
        cosx_siny_sinz_plus_sinx_cosz = cosx * siny * sinz + sinx * cosz
        sinx_cosy = sinx * cosy
        cosx_cosy = cosx * cosy


        for atomindex in self.AllAtoms:
            vector = self.AllAtoms[atomindex].coordinates

            new_x = vector.x * cosy_cosz + vector.y * sinx_siny_cosz_plus_cosx_sinz + vector.z * sinx_sinz_minus_cosx_siny_cosz
            new_y = -vector.x * cosy_sinz + vector.y * cosx_cosz_minus_sinx_siny_sinz + vector.z * cosx_siny_sinz_plus_sinx_cosz
            new_z = vector.x * siny - vector.y * sinx_cosy + vector.z * cosx_cosy

            self.AllAtoms[atomindex].coordinates = point(new_x, new_y, new_z)

        # now move it back from the origin
        self.TranslateMolecule(x, y, z)

    # Same as with the atom class
    def SetUndoPoint(self): # you can restore all atom positions to some undo point. This sets that point.
        for atomindex in self.AllAtoms:
            self.AllAtoms[atomindex].SetUndoPoint()

    def Undo(self):
        for atomindex in self.AllAtoms:
            self.AllAtoms[atomindex].Undo()