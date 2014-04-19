#===========================================================
def readcoor(filename,type,reportfile):

 f1=open(filename,'rb') # Open coordinate file
 array =[row.strip().split('\t') for row in f1] # put coordinates in a list

 del array[0:4] #removing 4 first rows of the list (Comment elements of coordinate file)
 atom=[]; matrix=[]; count = 0;i=0

 while (count<len(array)):
    atom.append(array[count][0].split())
    if atom[count][3]==type:matrix.append(atom[count][1:8]);i=i+1 # Extracting of [id, residu, 'CB', x, y, z, CH1] from Coor file
    count = count +1

 print >> reportfile,'\n\nRead atoms of type',type, ':',i, 'Out of', len(array), '\n'

 for row in matrix: reportfile.write(' '.join(row) + '\n') # Print a copy of extracted matrix into the report file
 print "Coordinations of", type, "atoms opened successfully!"
 return matrix
#============================================================
def reademin(filename,reportfile):
# this function reads emin values
 f1=open(filename,'rb') # Open emin file
 array =[row.strip().split('\t') for row in f1]
 
 emin=[]; count=0

 while (count<len(array)):
     emin.append(array[count][0].split())
     count=count+1

 print >> reportfile, "\n", filename, "values:\n"
 for row in emin: reportfile.write(' '.join(row) + '\n') # Print a copy of extracted matrix into the report file

 print "emin file opened successfully!"
 return emin
#============================================================
def resconvert(array,reportfile):
# this fucntion converts residue names to an index according to the kgs paper, table 3
# e.g. input of this function: [id, residu, 'CB', x, y, z, CH1]; residu index=1, CH1 index=6 
# output of function: [id, residu, 'CB', x, y, z, ***CODE***]
 count=0
 i=0
 gly=0;ala=0;ser=0;cys=0;val=0;thr=0;ile=0;pro=0;met=0;asp=0;asn=0;leu=0;lys=0;glu=0;gln=0;arg=0;hsd=0;his=0;phe=0;typ=0;trp=0;cyx=0
 print >> reportfile, '\nConverting residue types to codes ...' 
 while (count<len(array)):
    if array[count][1]=='GLY': array[count][6]= 0; gly=gly+1; i=i+1
    elif array[count][1]=='ALA': array[count][6]= 1; ala=ala+1; i=i+1
    elif array[count][1]=='SER': array[count][6]= 2; ser=ser+1; i=i+1
    elif array[count][1]=='CYS': array[count][6]= 3; cys=cys+1; i=i+1
    elif array[count][1]=='VAL': array[count][6]= 4; val=val+1; i=i+1
    elif array[count][1]=='THR': array[count][6]= 5; thr=thr+1; i=i+1
    elif array[count][1]=='ILE': array[count][6]= 6; ile=ile+1; i=i+1
    elif array[count][1]=='PRO': array[count][6]= 7; pro=pro+1; i=i+1
    elif array[count][1]=='MET': array[count][6]= 8; met=met+1; i=i+1
    elif array[count][1]=='ASP': array[count][6]= 9; asp=asp+1; i=i+1
    elif array[count][1]=='ASN': array[count][6]= 10; asn=asn+1; i=i+1
    elif array[count][1]=='LEU': array[count][6]= 11; leu=leu+1; i=i+1
    elif array[count][1]=='LYS': array[count][6]= 12; lys=lys+1; i=i+1
    elif array[count][1]=='GLU': array[count][6]= 13; glu=glu+1; i=i+1
    elif array[count][1]=='GLN': array[count][6]= 14; gln=gln+1; i=i+1
    elif array[count][1]=='ARG': array[count][6]= 15; arg=arg+1; i=i+1
    elif array[count][1]=='HSD': array[count][6]= 16; hsd=hsd+1; i=i+1
    elif array[count][1]=='HIS': array[count][6]= 16; his=his+1; i=i+1
    elif array[count][1]=='PHE': array[count][6]= 17; phe=phe+1; i=i+1
    elif array[count][1]=='TYR': array[count][6]= 18; typ=typ+1; i=i+1
    elif array[count][1]=='TRP': array[count][6]= 19; trp=trp+1; i=i+1
    elif array[count][1]=='CYX': array[count][6]= 20; cyx=cyx+1; i=i+1
    else: print 'ERROR: Nondefined residues in the coor file!';sys.exit(errno.EACCES)
    count=count+1
 print >> reportfile, 'Total number of converted residues to codes:', i
 print >> reportfile, 'Nomber of converted residues to code:', '\n0-GLY:', gly , '\n'# '\n1-ALA:', ala, '\n2-SER:', ser, '\n3-CYS', cys
 print "Index values assigned to all residues successfully!"
#============================================================
def writenbfix(array1,array2,statpot,cutoff,nbfix,reportfile):
 # This fucntion writes NBFIX section. It also considers GLY as a single residue.
 import math
 i=0;j=0;
 print >> reportfile, '\nDetailed version of NBFIX section between atoms of', array1[0][2], 'and', array2[0][2], ':\n'

 for i in range(len(array1)):
   for j in range(i, len(array2)):
      
       if (i==0 and j==0):     #Writing some reports to report file just to make sure everything is okay!!!
          print >> reportfile, "i counter for",  array1[0][2],  "atoms and j counters for",  array2[0][2], "atoms started counting from 0..."
          print >> reportfile, "\nExpected values of counter i:\n", range(len(array1)), "\n"
          print >> reportfile, "\nExpected values of counter j:\n", range(len(array2)), "\n"
          
       m=int(array1[i][0]);n=int(array2[j][0]) 
       if abs(m-n)>3: 
            dist= math.sqrt(math.pow(float(array1[i][3])-float(array2[j][3]),2) + \
                            math.pow(float(array1[i][4])-float(array2[j][4]),2) + \
                            math.pow(float(array1[i][5])-float(array2[j][5]),2))
            if dist<=cutoff:
           
              s = float(statpot[array1[i][6]][array2[j][6]]);
              e_min = round((1.+0.6*abs(s))*(-2.),2);
              e_min_rep = round(e_min*(2./3.),2);
              dist = round(dist,3);

              print >> nbfix, array1[i][2].replace("C","")+'%i' %(m), array2[j][2].replace("C","")+'%i' %(n), e_min, dist
#              elif s<0: print >> nbfix, array1[i][2].replace("C","")+'%i' %(m), array2[j][2].replace("C","")+'%i' %(n), e_min_rep, dist,'-1'
#              else: print >> nbfix, array1[i][2].replace("C","")+'%i' %(m), array2[j][2].replace("C","")+'%i' %(n), e_min, dist,'0'
              
             #Writing a report of NBFIX section!

              print >> reportfile, array1[i][2].replace("C","")+'%i' %(m), array2[j][2].replace("C","")+'%i' %(n), e_min, dist, \
                               array1[i][1], array2[j][1], array1[i][6], array2[j][6], '--> emin:', s, '--> E_min:',e_min  
#              elif s<0: print >> reportfile, array1[i][2].replace("C","")+'%i' %(m), array2[j][2].replace("C","")+'%i' %(n), e_min_rep, dist,'-1', \
#                               array1[i][1], array2[j][1],'--> emin:', s, '--> E_min_rep:', e_min_rep
#              else: print >> reportfile, array1[i][2].replace("C","")+'%i' %(m), array2[j][2].replace("C","")+'%i' %(n), e_min, dist,'0',  \
#                              array1[i][1], array2[j][1], '--> emin:', s, '--> E_min:',e_min

 print "NBFIX interactios of", array1[0][2], 'and', array2[0][2], "wrote successfully!"             
#=============================================================
def writenbonded(array1,array2,prmfile):

 prmfile=open('./nbonded.prm', 'w+')
 i=0;j=0;m=0;n=0
 for i in range(len(array1)):print >> prmfile, 'A%i'%(i+1),'0.0','-1e-12','22.7'
 for j in range(len(array2)):print >> prmfile, 'B%i'%int(array2[j][0]),'0.0','-1e-12','22.7'

#==============================================================
def writertf(array1,mass1,array2,mass2,topfile):
  i=0;j=0;m=0;n=0;residue="";reslist=[]

  print >> topfile, '* Topology file for a single bead GO model system'
  print >> topfile, '28 6'
  for i in range(len(array1)):
       m=int(array1[i][0])
       print >> topfile, 'MASS', (i+1), array1[i][2].replace("C","")+'%i' %(m), mass1

  for j in range(len(array1), len(array1) + len(array2)):
       n=int(array2[j-len(array1)][0])
       print >> topfile, 'MASS', (j+1), array2[j-len(array1)][2].replace("C","")+'%i' %(n), mass2
 
  print >> topfile, "\nDECL -CA" 
  print >> topfile, "DECL +CA" 
  print >> topfile, "DECL #CA ! '#'=second next atom\n"

  for i in range(len(array1)):
      residue=str(array1[i][1])
      m=int(array1[i][0])
      print >> topfile, "\nRESI", "P%i" %m, "0.00"
      print >> topfile, "GROUP"
      print >> topfile, "ATOM", "CA", "A%i" %m, "0.00"
      if residue <> "GLY": print >> topfile, "ATOM", "CB", "B%i" %m, "0.00" #GLY does NOT contain CB
      if i<len(array1)-1 :print >> topfile, "BOND", "CA", "+CA"
      if residue <> "GLY": print >> topfile, "BOND", "CA", "CB"
      if (0 < i < len(array1)-2) and (residue <> "GLY"):
         print >> topfile, "DIHE -CA CA +CA #CA"
         print >> topfile, "IMPR CA -CA +CA CB"
      if i==len(array1)-2: print >> topfile, "IMPR CA -CA +CA CB"
 
  print "Topology file is ready to use!"


