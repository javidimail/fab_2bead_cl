#! /usr/bin/python
import math
import errno, sys
import time
import socket
import getpass
from go_functions import readcoor, reademin, resconvert, writenbfix, writertf, writenbonded

reportfile=open('./report.txt','w+') # A report file to write some of the selected results

localtime = time.asctime(time.localtime(time.time()))
print >> reportfile, "Local current time :", localtime, '\nComputar host name:', socket.gethostname(), '\nUser:', getpass.getuser(), '\n'

ca=[]; cb=[]; eminb=[]  

ca=readcoor('cl_ca_cb.cor','CA',reportfile) # Extracting CA atoms from coodinate file (6 elements), e.g. [id, residu, 'CA', x, y, z, CH1]
cb=readcoor('cl_ca_cb.cor','CB',reportfile) # Same as CA,e.g. [id, residu, 'CB', x, y, z, CH1]
eminb=reademin('emin_b.dat',reportfile) # Reading stat potential file

resconvert(ca,reportfile) # Assigning a code for each residu, e.g, [id, residu, 'CB', x, y, z, ***CODE***]
resconvert(cb,reportfile) 

nbonded=open('./nbonded.prm', 'w+')
writenbonded(ca,cb,nbonded)

nbfix=open('./nbfix.prm','w+') # it opens nbfix.prm to write NBFIX section
writenbfix(ca,ca,eminb,8,nbfix,reportfile) # Writing NBFIX section for CA-CA's with cutoff=8 A 
writenbfix(cb,cb,eminb,8,nbfix,reportfile) # doing the same for CB-CB's
writenbfix(cb,ca,eminb,8,nbfix,reportfile) # doing the same for CB-CA's, **** Smaller c, should be first ****


rtf=open('./ch1_2bead.rtf','w+') # it opens a file to write the topology file 
writertf(ca,100,cb,200,rtf) # Writing RTF file, MASS CA=100, MASS CB=200

# Merge all parts of parameter file
prmfiles = ['bond.prm', 'angle.prm', 'dihedral.prm','imph.prm','nbonded.prm','nbfix.prm']
with open('cl_2bead.prm', 'w') as outfile:
    for fname in prmfiles:
        with open(fname, 'rb') as infile:
            fname2=fname.replace('.prm','')
            print >> outfile, '\n', fname2.upper(), '\n'
            for line in infile:
              outfile.write(line)

