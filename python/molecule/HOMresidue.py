import __init__
import molecule.HOMatom

class Residue:
    def __init__(self):
        self.resname=''
        self.resid=0
        self.chain=''
        self.CA=Atom()
        self.attr=''
        self.Allatoms={}