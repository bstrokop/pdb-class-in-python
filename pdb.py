import numpy as np
class Pdb():
    

    def __init__(self, name):
        self.name = name
        dd = self.pdb_init()
        print(dd)
        self.read_pdb()
        #charar = np.chararray(10, itemsize=5, unicode=True)
        self.remark         = np.chararray(shape=(self.nrem,),    itemsize=80, unicode=True)
        #atoms
        self.label          = np.chararray(shape=(self.natom,),   itemsize=6,   unicode=True)
        self.atom_number    = np.ndarray  (shape=(self.natom,),   dtype=int)
        self.atom_name      = np.chararray(shape=(self.natom,),   itemsize=4,  unicode=True)
        self.altloc         = np.chararray(shape=(self.natom,),   itemsize=1,  unicode=True)
        self.residue_name   = np.chararray(shape=(self.natom,),   itemsize=3,  unicode=True)
        self.chain_id       = np.chararray(shape=(self.natom,),   itemsize=1,  unicode=True)
        self.residue_number = np.ndarray  (shape=(self.natom,),   dtype=int)
        self.insertion_code = np.chararray(shape=(self.natom,),   itemsize=1,  unicode=True)
        self.xyz            = np.ndarray  (shape=(3, self.natom),   dtype=np.float64)
        self.occ            = np.ndarray  (shape=(self.natom,),     dtype=np.float64)
        #print (np.size(self.occ), self.natom)
        #exit(1)
        self.biso           = np.ndarray(shape=(self.natom,),     dtype=np.float64)
        self.segid          = np.chararray(shape=(self.natom,), itemsize=4, unicode=True)
        self.elsymb         = np.chararray(shape=(self.natom,), itemsize=2, unicode=True)
        self.charge         = np.chararray(shape=(self.natom,), itemsize=2, unicode=True)       

        self.ulabel         = np.chararray(shape=(self.natom,), itemsize=6, unicode=True)

#       atom blocks     
        self.indxyz         = np.ndarray(shape=(self.natom,),    dtype=int)
        self.atom_type      = np.ndarray(shape=(self.natom,),    dtype=int)
        self.refinable      = np.ndarray(shape=(11, self.natom), dtype=bool)

       

#      group business FIXME need more data
#       self.ngroups       = ngroup
#       self.residue_start  = np.array((2,), dtype=int)
#       self.chain_start    = np.array(nchains+1, dtype=int)
        self.cell           = np.ndarray(shape=(6,), dtype=np.float64)
        self.ORT            = np.ndarray(shape=(3,3), dtype=np.float64)
        self.DEORT          = np.ndarray(shape=(3,3), dtype=np.float64)

#       continue a bit with aniso
        self.aniso_in_pdb   = False
        self.mixed  = False
        if self.naniso > 0:
            self.aniso_in_pdb = True
            self.u = np.ndarray(shape=(6, self.natom), dtype=np.float64)
            self.mixed = self.natom > self.naniso
       
        print(len(self.label), len(self.remark), 'mixed_ref=', self.mixed)

#       Test float64 and higher:
        a = 1 + np.finfo(np.longdouble).eps
        b = 1 + np.finfo(np.double).eps

        if  a != 1.0 and b != 1.0 :
            print('The world of numpy double and long double floats is good.')
        else:
            print('The world of numpy double and long double floats is not so good.')
            print(a,b)
            exit(2)

        iat = 0
        for line in self.lines:
            line = line.replace('\n','') 

            if len(line) > 5:
                if line[0:6] == 'ATOM  ' or line[0:6] == 'HETATM':
                    self.label[iat] = line[0:6]
                    #print(' iat=', iat, self.label[iat])
                    for key in dd:

                        d1, d2 = dd[key]

                        # check line length:
                        if  d2  > len(line):
                            print('Error in line:')
                            print(line)
                            print('Wanna access field no.', d2, 'but line length is just',len(line))
                            exit(2)

                        d1 = d1 - 1
                        if key == 'ALABEL':
                            self.label[iat] = line[d1:d2]
                        elif key == 'ATOM_NUMBER':
                            self.atom_number[iat] = line[d1:d2]
                        elif key == 'ATOM_NAME':
                            self.atom_name[iat] = line[d1:d2] 
                        elif key == 'ALTLOC':
                            self.altloc[iat] = line[d1:d2]
                        elif key == 'RESNAM':
                            self.residue_name[iat] = line[d1:d2]
                        elif key == 'CHAIN_ID':
                            self.chain_id[iat] = line[d1:d2]
                        elif key == 'RESNUM':
                            self.residue_number[iat] = line[d1:d2]
                        elif key == 'ICODE':
                            self.insertion_code[iat] = line[d1:d2]
                        elif key == 'X':
                            self.xyz[0,iat] = line[d1:d2]
                        elif key == 'Y':
                            self.xyz[1,iat] = line[d1:d2]
                        elif key == 'Z':
                            self.xyz[2,iat] = line[d1:d2]
                        elif key == 'OCC':
                            #print(line[d1:d2])
                            self.occ[iat] = float(line[d1:d2])
                            #exit(1)
                        elif key == 'BISO':
                            self.biso[iat] = line[d1:d2]
                        elif key == 'SEGID':
                            self.segid[iat] = line[d1:d2]
                        elif key == 'ELEMENT':
                            self.elsymb[iat] = line[d1:d2]
                        elif key == 'CHARGE':
                            self.charge[iat] = line[d1:d2]

                    #print('iat=',iat, "[",self.charge[iat],"]", type(self.occ[iat]), self.natom, len(self.occ), len(self.atom_number))
                    iat += 1

                elif line[0:6] == 'ANISOU':
                    
                    self.ulabel[iat] = line[0:6]                    

                    for key in dd:
                        d1, d2 = dd[key] 
                        d1 = d1 - 1
                        
                        if  key == 'U11':
                            self.u[1,iat] = line[d1:d2]
                        elif key == 'U22':
                            self.u[2,iat] = line[d1:d2]
                        elif key == 'U33':
                            self.u[3,iat] = line[d1:d2]
                        elif key == 'U12':
                            self.u[4,iat] = line[d1:d2]
                        elif key == 'U13':
                            self.u[5,iat] = line[d1:d2]
                        elif key == 'U23':
                            self.u[6,iat] = line[d1:d2]


                          

    def pdb_init(self):

        au_dict = {    "ALABEL"         : [1,6],     # ATOM / HETATM label                 String
                       "ATOM_NUMBER"    : [7,11],    # atom serial number                  Int
                       "ATOM_NAME"      : [13,16],   # atom name                           String
                       "ALTLOC"         : [17,17],   # alt location A/B business           String
                       "RESNAM"         : [18,20],   # residue name                        String
                       "CHAIN_ID"       : [22,22],   # chain id                            String
                       "RESNUM"         : [23,26],   # residue sequence number             Int
                       "ICODE"          : [27,27],   # Code for insertion of residues      String
                       "X"              : [31,38],   # x coordinate (orthogonal)           Float
                       "Y"              : [39,46],   # y coordinate                        Float
                       "Z"              : [47,54],   # z coordinate                        Float
                       "OCC"            : [55,60],   # occupancy record                    Float
                       "BISO"           : [61,66],   # temperature factor (isotropic)      Float
                       "SEGID"          : [73,76],   # segment id (for MR and rigid-body)  String
                       "ELEMENT"        : [77,78],   # element symbol (right justified)    String
                       "CHARGE"         : [79,80],   # atom charge (optional)              String
                       "U11"            : [29,35],   # U11 (orthogonal)                    Int
                       "U22"            : [36,42],   # U22                                 Int
                       "U33"            : [43,49],   # U33                                 Int
                       "U12"            : [50,56],   # U12                                 Int
                       "U13"            : [57,63],   # U13                                 Int
                       "U23"            : [64,70],   # U23                                 Int
                  }

        return au_dict

    def write_pdb(self, fileout):
        print('Gonna print PDB file:', fileout)
        #super long new style PDB format:
        atom_fmt = '{0:<6s}{1:>5d} {2:4s}{3:1s}{4:3s} {5:1s}{6:>4d}{7:1s}   {8:8.3f}{9:8.3f}{10:8.3f}{11:6.2f}{12:6.2f}      {13:4s}{14:>2s}{15:2s}' + '\n'
        # need to learn aniso fmt yet: FIXME

        with open(fileout, 'wt') as fout:
            iat = 0
            for line in self.lines:
                #print(line)
                if len(line) > 5:
                    if line[0:6] == 'ATOM  ' or line[0:6] == 'HETATM':
                         fout.write(atom_fmt.format(self.label[iat],          #0
                                                    self.atom_number[iat],    #1 
                                                    self.atom_name[iat],      #2
                                                    self.altloc[iat],         #3
                                                    self.residue_name[iat],   #4
                                                    self.chain_id[iat],       #5
                                                    self.residue_number[iat], #6
                                                    self.insertion_code[iat], #7
                                                    self.xyz[0,iat],          #8 
                                                    self.xyz[1,iat],          #9
                                                    self.xyz[2,iat],          #10
                                                    self.occ[iat],            #11
                                                    self.biso[iat],           #12
                                                    self.segid[iat],          #13
                                                    self.elsymb[iat],         #14
                                                    self.charge[iat]          #15
                             ) )  
                
                     
                         iat += 1
                    else:
                        # print all other lines except short ones
                        fout.write(line)
                else:
                    # very short line printing like 'TER' or 'END' to preserver structure of the pdb file
                    fout.write(line)

        fout.close()
        print('Closing pdbout file...')


    def read_pdb(self):
        fin = open(self.name, 'rt')
        self.lines = fin.readlines()
        fin.close()
        print(len(self.lines), 'lines read')
        natom = 0
        naniso = 0
        nsigatm = 0
        nsiguij = 0
        nrem = 0

        for line in self.lines:

            #cleaning line a bit
            line = line.replace('\n','') 
            line = line.replace('\t','') 

            #print(line, "[", len(line),"]")
            s6 = ''.join(line[:6])  # upper index in silces should be a bit larger (by 1)

            if  s6 == 'ATOM  ' or s6 == 'HETATM':
                print(line, "[", len(line),"]")
                natom += 1
            elif s6 == 'ANISOU':
                print(line, "[", len(line),"]")
                naniso += 1
            elif s6 == 'SIGATM':
                nsigatm += 1
            elif s6 == 'SIGUIJ':
                nsiguij += 1
            else:
                nrem += 1  # all the rest will goto remarks 
                #print('No match:', s6, "|", line[0:5], "|",len(s6), type('ATOM  '), len('ATOM  '))
                #if len(line) > 5:
                #    print(line[0], line[5])
        print()

        self.natom = natom
        self.naniso = naniso
        self.nsigatm = nsigatm
        self.nsiguij = nsiguij       
        self.nrem    = nrem
        print(self.natom,'atoms has been read.')
        print(self.nrem,'remarks has been saved.')


def main():
    import argparse
    from datetime import datetime as dt

    parser = argparse.ArgumentParser(description="Analyze and read I/O files")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-v", "--verbose", action="count", default=0)
    group.add_argument("-q", "--quiet", action="store_true")
    parser.add_argument("-pdb", "--pdbin", help="Input PDB file")
    parser.add_argument("-pdbout", "--pdbout", help="Output PDB file")
    args = parser.parse_args()
#    answer = args.x**args.y
    print("pdbin={} pdbout={}".format(args.pdbin, args.pdbout))

#    if args.quiet:
#        print(answer)
#    elif args.verbose:
#        print("{} to the power {} equals {}".format(args.x, args.y, answer))
#    else:
#        print("{}^{} == {}".format(args.x, args.y, answer))


#if __name__ == "__main__":
#main(sys.argv[1:])

    now = dt.now()
    print(f"Текущее время {now:%d.%m.%Y %H:%M}")
    pdb_1 = Pdb(args.pdbin)
    pdb_1.write_pdb(args.pdbout)
#    read_pdb(args.pdbin)

main()
