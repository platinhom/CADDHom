#import numpy as np
import math

print "I'm in point"

class Point:
    x=99999.0
    y=99999.0
    z=99999.0    
    def __init__ (self, x, y ,z):
        self.x = x
        self.y = y
        self.z = z
        
    def coors(self):
        coor=(self.x,self.y,self.z)
        return coor
        
    def dist_to(self,apoint):
        return math.sqrt(math.pow(self.x - apoint.x,2) + math.pow(self.y - apoint.y,2) + math.pow(self.z - apoint.z,2))

    def CopyOf(self):
        return point(self.x, self.y, self.z)
    
    def average_with(self, other_point):
        return point((self.x + other_point.x) / 2.0, (self.y + other_point.y) / 2.0, (self.z + other_point.z) / 2.0)
    
    def dot_product_with(self, other_point):
        return self.x * other_point.x + self.y * other_point.y + self.z * other_point.z
    
    def length(self):
        return self.dist_to(point(0.0,0.0,0.0))
    
    def minus(self, other_point):
        return point(self.x - other_point.x, self.y - other_point.y, self.z - other_point.z)
    
    def CreatePDBLine(self):
        #if len(self.atomname) > 1: self.atomname = self.atomname[:1].upper() + self.atomname[1:].lower()
        output = "ATOM "
        #output = output + str(index).rjust(6) + self.atomname.rjust(5) + self.residue.rjust(4)
        output = output + "5".rjust(6) + "X".rjust(5) + "XXX".rjust(4)
        output = output + ("%.3f" % self.x).rjust(18)
        output = output + ("%.3f" % self.y).rjust(8)
        output = output + ("%.3f" % self.z).rjust(8)
        output = output + "X".rjust(24) # + "   " + str(uniqueID) #This last part must be removed
        return output  



#######basic function is put here#####
def coors_to_point(coors):
    point=Point(coors[0],coors[1],coors[2])
    return point