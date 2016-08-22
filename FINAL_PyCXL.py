
# coding: utf-8

# In[3]:

import re
from operator import itemgetter
from tabulate import tabulate #### must install! pip install tabulate
import os
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MultipleLocator

### GLOBAL SETTINGS

H = 1.00782503207
H2O = 18.0105650641
NH3 = 17.026549101

plt.close("all")
np.set_printoptions(suppress=True)


###################


A55 = ['A55','MEVNKKQLADIFGASIRTIQNWQEQGMPVLRGGGKGNEVLYDSAAVIKWYAERDA']
terminase = ['TerL','MNISNSQVNRLRHFVRAGLRSLFRPEPQTAVEWADANYYLPKESAYQEGRWETLPFQRAIMNAMGSDYIREVNVVKSARVGYSKMLLGVYAYFIEHKQRNTLIWLPTDGDAENFMKTHVEPTIRDIPSLLALAPWYGKKHRDNTLTMKRFTNGRGFWCLGGKAAKNYREKSVDVAGYDELAAFDDDIEQEGSPTFLGDKRIEGSVWPKSIRGSTPKVRGTCQIERAASESPHFMRFHVACPHCGEEQYLKFGDKETPFGLKWTPDDPSSVFYLCEHNACVIRQQELDFTDARYICEKTGIWTRDGILWFSSSGEEIEPPDSVTFHIWTAYSPFTTWVQIVKDWMKTKGDTGKRKTFVNTTLGETWEAKIGERPDAEVMAERKEHYSAPVPDRVAYLTAGIDSQLDRYEMRVWGWGPGEESWLIDRQIIMGRHDDEQTLLRVDEAINKTYTRRNGAEMSISRICWDTGGIDPTIVYERSKKHGLFRVIPIKGASVYGKPVASMPRKRNKNGVYLTEIGTDTAKEQIYNRFTLTPEGDEPLPGAVHFPNNPDIFDLTEAQQLTAEEQVEKWVDGRKKILWDSKKRRNEALDCFVYALAALRISISRWQLDLSALLASLQEEDGAATNKKTLADYARALSGEDE','TerS','MEVNKKQLADIFGASIRTIQNWQEQGMPVLRGGGKGNEVLYDSAAVIKWYAERDAEIENEKLRREVEELRQASEADLQPGTIEYERHRLTRAQADAQELKNARDSAEVVETAFCTFVLSRIAGEIASILDGLPLSVQRRFPELENRHVDFLKRDIIKAMNKAAALDELIPGLLSEYIEQSG']
IHF = ['IHF_B','MTKSELIERLATQQSHIPAKTVEDAVKEMLEHMASTLAQGERIEIRGFGSFSLHYRAPRTGRNPKTGDKVELEGKYVPHFKPGKELRDRANIYG','IHF_A', 'MALTKAEMSEYLFDKLGLSKRDAKELVELFFEEIRRALENGEQVKLSGFGNFDLRDKNQRPGRNPKTGEDIPITARRVVTFRPGQKLKSRVENASPKDE']
Cyt = ['CytB5', 'MAEQSDEAVKYYTLEEIQKHNHSKSTWLILHHKVYDLTKFLEEHPGGEEVLREQAGGDATENFEDVGHSTDAREMSKTFIIGELHPDDRPKLNKPPETLITTIDSSSSWWTNWVIPAISAVAVALMYRLYMAED', 'Cyt2E1', 'MSALGVTVALLVWAAFLLLVSMWRQVHSSWNLPPGPFPLPIIGNLFQLELKNIPKSFTRLAQRFGPVFTLYVGSQRMVVMHGYKAVKEALLDYKDEFSGRGDLPAFHAHRDRGIIFNNGPTWKDIRRFSLTTLRNYGMGKQGNESRIQREAHFLLEALRKTQGQPFDPTFLIGCAPCNVIADILFRKHFDYNDEKFLRLMYLFNENFHLLSTPWLQLYNNFPSFLHYLPGSHRKVIKNVAEVKEYVSERVKEHHQSLDPNCPRDLTDCLLVEMEKEKHSAERLYTMDGITVTVADLFFAGTETTSTTLRYGLLILMKYPEIEEKLHEEIDRVIGPSRIPAIKDRQEMPYMDAVVHEIQRFITLVPSNLPHEATRDTIFRGYLIPKGTVVVPTLDSVLYDNQEFPDPEKFKPEHFLNENGKFKYSDYFKPFSTGKRVCAGEGLARMELFLLLCAILQHFNLKPLVDPKDIDLSPIHIGFGCIPPRYKLCVIPRS']



### USER INPUT

sequenceoverride = 'n'
wholesequence = 'GGGKGNEVLYDSAAVIKGGGKGNEVLYDSAAVIK'
mgffile = '/home/easyprot/Desktop/Windows_Shared/A55-IHF_final_cxl_biotin/A55_BS2/unmatched2/TerS_A55_lo_BS2_unmatched.mgf_wscans.mgf'
scantitle = 1364.4
protein = A55 ## not a string

isotope_tolerance = 0.005
ppmlimit = 2000
minpercentintense = 1

showspectra = 'y'
numspecs = 3
numpeaksperlist = 50


#### NOTE protein is NOT a string! END USER INPUT!!                          


##### Begin defininte pepA, pepB, cxl1, cxl2, linkmass

if re.findall(r'BS2', mgffile) != []:

    linkmass = 96.0211301283
    
if re.findall(r'EDC', mgffile) != []:
    
    linkmass = -18.0105650641
    
##splice out mods titles from EasyProt

sequence = ''

i = 0

for element in re.split(r'\W*', wholesequence):

    if element == 'BS2_K' or element == 'Oxidation_M' or element == 'Cys_CAM':
        
        i = i + 1
        
    else:
        
        sequence = sequence + element
        
        i = i + 1
        
#####################################################################################
### check which protein the first 3,4,5 aas in wholesequence match which protein ####
#####################################################################################

if sequenceoverride == 'n':

    begsequcheck = []

    i = 0

    for prot in protein:

        if re.findall(sequence[:3],prot) != []:

            begsequcheck.append(i)

        i = i + 1

    i = 0

    if len(begsequcheck) == 2:

        begsequcheck = []

        for prot in protein:

            if re.findall(sequence[:4],prot) != []:

                begsequcheck.append(i)

            i = i + 1

    i = 0

    if len(begsequcheck) == 2:

        begsequcheck = []

        for prot in protein:

            if re.findall(sequence[:5],prot) != []:

                begsequcheck.append(i)

            i = i + 1

    pepA = ''
    temppep = ''
    split = 0


    for aa in list(sequence):

        temppep = temppep + aa    

        if re.findall(temppep, protein[begsequcheck[0]]) == []:

            pepA = sequence[:split]

            break

        split = split + 1
    
##assign starting protein index and protein name to pepA

protA = protein[begsequcheck[0]-1]
startA = protein[begsequcheck[0]].find(pepA)

## double check pepA ends in K or R and if it isn't, it's C terminal

listsequA = list(pepA)

if re.findall('K', str(listsequA[len(pepA)-1])) == []:    

    if re.findall('R', str(listsequA[len(pepA)-1])) == []:
        
        pepAend = int(startA + len(pepA))
                
        if len(list(protein[begsequcheck[0]])) != pepAend:
                       
            pepA = pepA[:len(pepA)-1]

           
## assign pepB

pepB = ''.join(re.split(pepA, sequence))

### done with pepB, now assign protB, startB

i = 0

for prot in protein:
    
    if re.findall(pepB, prot) != []:
        
        protB = str(protein[i-1])
        startB = prot.find(pepB)
        
    i = i + 1

peptideA = pepA
peptideB = pepB

####################################################################    
############# now determine posA and posB options ##################
####################################################################

if re.findall(r'BS2', mgffile) != []:

    pepA_Ks = []

    i = 0

    for aa in list(pepA):

        if aa == 'K':

            pepA_Ks.append(i + 1)

        i = i + 1
        
    if startA + 1 == 1:
        
        pepA_Ks.append(1)


    pepB_Ks = []

    i = 0

    for aa in list(pepB):

        if aa == 'K':

            pepB_Ks.append(i + 1)

        i = i + 1
        
    if startB + 1 == 1:
        
        pepB_Ks.append(1)

    cxl1 = []
    cxl2 = []

    for K1 in pepA_Ks:

        for K2 in pepB_Ks:

            cxl1.append([K1, K2])
            
            if K1 != 1 and K2 != 1:
                
                cxl2.append(['K' + str(K1+int(startA)),'K' + str(K2+int(startB))])
                
            if K1 != 1 and K2 ==1:
                
                cxl2.append(['K' + str(K1+int(startA)),'nterm'])
                
            if K1 == 1 and K2 != 1:
                
                cxl2.append(['nterm','K' + str(K2+int(startB))])
                
            if K1 == 1 and K2 == 1:
                
                cxl2.append(['nterm','nterm'])

            
if re.findall(r'EDC', mgffile) != []:

    pepA_Ks = []
    pepA_Es = []
    pepA_Ds = []
    fullprotA = protein[protein.index(protA)+1]
    fullprotB = protein[protein.index(protB)+1]
    
    i = 0

    for aa in list(pepA):

        if aa == 'K':

            pepA_Ks.append(i + 1)
            
        if aa == 'E':
            
            pepA_Es.append(i + 1)
            
        if aa == 'D':
            
            pepA_Ds.append(i + 1)

        i = i + 1
        
    if (startA + 1) == 1:
        
        pepA_Ks.append(1)
        
    if (startA + len(pepA)) == len(fullprotA):
        
        pepA_Es.append(len(fullprotA))
        
    pepB_Ks = []
    pepB_Es = []
    pepB_Ds = []
    
    i = 0

    for aa in list(pepB):

        if aa == 'K':

            pepB_Ks.append(i + 1)
            
        if aa == 'E':
            
            pepB_Es.append(i + 1)
            
        if aa == 'D':
            
            pepB_Ds.append(i + 1)

        i = i + 1
        
    if startB + 1 == 1:
        
        pepB_Ks.append(-1)
        
    if startB + len(pepB) == len(fullprotB):
        
        pepB_Es.append(len(fullprotB))
    
    cxl1 = []
    cxl2 = []
    
    for aa1 in pepA_Ks:
        
        if aa1 == -1:
        
            for aa2 in pepB_Es:
                
                if aa2 == len(fullprotB):

                    cxl1.append([aa1, len(pepB)])
                    cxl2.append(['nterm', 'cterm'])
                    
                else:
                    
                    cxl1.append([aa1, aa2])
                    cxl2.append(['nterm', 'E' + str(aa2 + startB)])

            for aa3 in pepB_Ds:

                cxl1.append([aa1, aa3])
                cxl2.append(['nterm', 'D' + str(aa3 + startB)])
                
        else:
            
            i = 0
            
            for aa2 in pepB_Es:
                
                if aa2 == (len(fullprotB)):

                    cxl1.append([aa1, len(pepB)])
                    cxl2.append(['K'+ str(aa1 + startA), 'cterm'])
                
                else:
                    
                    cxl1.append([aa1, aa2])
                    cxl2.append(['K'+ str(aa1 + startA), 'E' + str(aa2 + startB)])
                    
            for aa3 in pepB_Ds:
                
                cxl1.append([aa1, aa3])
                cxl2.append(['K'+ str(aa1 + startA), 'D' + str(aa3 + startB)])
            
    for aa1 in pepB_Ks:
        
        if aa1 == -1:
        
            for aa2 in pepA_Es:
                
                if aa2 == (len(fullprotA)):

                    cxl1.append([len(pepA), aa1])
                    cxl2.append(['cterm', 'nterm'])
                    
                else:
                    
                    cxl1.append([aa2, aa1])
                    cxl2.append(['E' + str(aa2 + startA), 'nterm'])

            for aa3 in pepA_Ds:

                cxl1.append([aa3, aa1])
                cxl2.append(['D' + str(aa3 + startA), 'nterm'])
                
        else:
            
            for aa2 in pepA_Es:
                
                if aa2 == (len(fullprotA)):

                    cxl1.append([len(pepA), aa1])
                    cxl2.append(['cterm', 'K' + str(aa1 + startB)])
                    
                else:
                    
                    cxl1.append([aa2, aa1])
                    cxl2.append(['E' + str(aa2 + startA), 'K' + str(aa1 + startB)])

            for aa3 in pepA_Ds:

                cxl1.append([aa3, aa1])
                cxl2.append(['D' + str(aa3 + startA), 'K' + str(aa1 + startB)])       
            
    #for aa1 in pepB_Ks:
        
        #for aa2 in pepA_Es:
            
            #cxl1.append([aa2, aa1])
            #cxl2.append(['E'+ str(aa2 + startA), 'K' + str(aa1 + startB)])
            
        #for aa3 in pepA_Ds:
            
            #cxl1.append([aa3, aa1])
            #cxl2.append(['D'+ str(aa3 + startA), 'K' + str(aa1 + startB)])
            
#print pepA
#print pepB

#print cxl1
#print cxl2

#print len(fullprotA)
#print len(fullprotB)

#print startA
#print startB

#print pepA_Es
        
###########################################################################
##### COLLECT PEPTIDE A AND B MASSES FOR THEOR FRAG GENERATION LATER ######

AAmass = 'A 71.037114,R 156.101111,N 114.042927,D 115.026943,C 103.009185,E 129.042593,Q 128.058578,G 57.021464,H 137.058912,I 113.084064,L 113.084064,K 128.094963,M 131.040485,F 147.068414,P 97.052764,S 87.032028,T 101.047679,W 186.079313,Y 163.06332,V 99.068414,'


pepBmass = 0
    
for let in list(peptideB):

    pepBmass = pepBmass + float(''.join(re.findall(let + r'\s(\S*?)\,', AAmass)))
    
pepAmass = 0
    
for let in list(peptideA):

    pepAmass = pepAmass + float(''.join(re.findall(let + r'\s(\S*?)\,', AAmass)))
    

#### FIGURE OUT MODS!!!!!!



f = open(mgffile)
mgf = f.read()

precursorinfo = '\n'.join(re.findall(r'.*\n.*\n'+ r'SCANS\=' + str(scantitle), mgf))

pcharge = float(''.join(re.findall(r'CHARGE\=(\d)\+', precursorinfo)))
pmz = float(''.join(re.findall(r'PEPMASS\=(\S*)', precursorinfo)))

pnm = (pmz * pcharge) - (pcharge * H)

nomodmass = pepAmass + pepBmass + 2 * H2O + linkmass

modmassdiff = abs(nomodmass-pnm)

print '\n'
print 'experimental neutral mass = ' + str(pnm)
print 'unmodified neutral mass = ' + str(nomodmass)
print 'precursor mass difference = ' + str(modmassdiff)

finalmod = []

if modmassdiff < float(1):
    
    countm = 0
    countc = 0

else:

    ms = re.findall(r'M', peptideA) + re.findall(r'M', peptideB)
    cs = re.findall(r'C', peptideA) + re.findall(r'C', peptideB)
    
    numms = range(len(ms))
    numcs = range(len(cs))
   
    
    if numms != []:
        
        for m in numms:
            
            finalmod.append([m+1,0,abs(((m+1)*15.99491)-modmassdiff)])
            
    if numcs != []:
        
        for c in numcs:
            
            finalmod.append([0,c+1,abs(((c+1)*57.02146)-modmassdiff)])
            
    if numms != [] and numcs != []:
        
        tempmass = 0
        
        for m in numms:
            
            for c in numcs:
                
                tempmass = (m+1)*15.99491 + (c+1)*57.02146
                
                finalmod.append([m+1,c+1,abs(tempmass-modmassdiff)])

    finalmod = sorted(finalmod, key=itemgetter(2), reverse = False)
    
    countm = finalmod[0][0]
    countc = finalmod[0][1]
            
print 'oxidized Met # = ' + str(countm)
print 'carbamidomethyl Cys # = ' + str(countc)
print '\n'

###########################################################################
########## DEISOTOPE SCAN #################################################
###########################################################################



prescan = '\n'.join(re.findall(r'SCANS\=' + str(scantitle) + r'.*?IONS', mgf, re.DOTALL))

prescanlines = prescan.split('\n')
del prescanlines[0:2]
del prescanlines[len(prescanlines)-1]

entries = len(prescanlines)

i = 0

for line in prescanlines:
    
    prescanlines[i] = line.replace(' ', '\t')
    
    i= i + 1

scan = np.fromstring('\n'.join(prescanlines), dtype=float, sep='\t').reshape(len(prescanlines), 2)

#### show mass spectrum

scan_trans = np.transpose(scan)

### deisotoped scan
z = np.zeros((entries, 1))

dscan = np.append(scan, z, axis = 1)

i = 0

while (i < len(prescanlines)):
    
    if dscan[i][2] != 0:
    
        i = i + 1
        
    else:
        
        peakcharge = 0
        chargestate = []
        
        for num in range(15): ##define charges to check
            
            chargestate.append(1/(float(num+1)))
        
        ################################### gather appropriate peaks to compare
        
        tempmz = dscan[i][0]
        tempintense = dscan[i][1]
        temp = [dscan[i].tolist()]
        peakindex = [0]
        
        j = 1
        
        while i < (len(prescanlines) -10) and j < 11:
                        
            if (dscan[i + j][1]) > (0.333 * tempintense) and (dscan[i + j][0] - tempmz)  < (1 + isotope_tolerance):
                
                temp.append(dscan[i+j].tolist())
                tempmz = dscan[i+j][0]
                tempintense = dscan[i+j][1]
                peakindex.append(j)

                j = j + 1

            else:

                j = j + 1
                
        while i >= (len(prescanlines) -10) and j < (len(prescanlines) - i):

            if (dscan[i + j][1]) > (0.333 * tempintense) and (dscan[i + j][0] - tempmz)  < (1 + isotope_tolerance):
                
                temp.append(dscan[i+j].tolist())
                tempmz = dscan[i+j][0]
                tempintense = dscan[i+j][1]
                peakindex.append(j)

                j = j + 1

            else:

                j = j + 1
                
        ### now find peak differences
        
        d = []
        
        j = 1
        
        while j < (len(temp)):
            
            d.append(temp[j][0]-temp[j-1][0])
            
            j = j + 1
            
        ### now check differences and assign charge
        
        tempcharge = 0
        d2 = []
        grouping = []
        peakcharge = 0
        
        if d == []:
            
            dscan[i][2] = -2 ### can't assign charge state
 
            
        elif len(d) == 1:
            
            for charge in chargestate:
                
                if d[0] < (charge + isotope_tolerance) and d[0] > (charge - isotope_tolerance):
                    
                    tempcharge = 1/charge
                    dscan[i][2] = 1/charge
                    dscan[i+peakindex[1]][2] = -1
                    
                    break
                    
                if tempcharge == 0:
                    
                    dscan[i][2] = -2
            
        else:  
            
            for num in range(len(d)-1):

                d2.append(d[num+1]-d[num])
                
            for item in d2:
                
                if abs(item) < isotope_tolerance:
                    
                    grouping.append(1)
                    
                else:
                    
                    grouping.append(0)
            j = 0
            
            if grouping[0] == 0:
                
                dscan[i][2] = -2
                
            if grouping[0] == 1:
                
                for charge in chargestate:
                
                        if d[0] < (charge + isotope_tolerance) and d[0] > (charge - isotope_tolerance):
                        
                            peakcharge = 1/charge
                    
                            dscan[i][2] = peakcharge
                            dscan[i + peakindex[1]][2] = - 1
                            dscan[i + peakindex[2]][2] = - 1
                            
                if peakcharge == 0:
                    
                    dscan[i][2] = -2
                
            if grouping[0] == 1 and len(grouping) > 1 and peakcharge != 0:                            
                             
                for num in range(1, len(grouping)):
                
                    if grouping[num] == 1:
    
                        dscan[i + peakindex[num+2]][2] = -1
                        
                    else: 
                        
                        break
        #print dscan[i]
        #print temp
        #print d
        #print d2
        
        i = i + 1

#print dscan        
        
dscanfinal = []

for entry in dscan:
    
    peak = entry.tolist()
    
    if peak[2] != -1:
        
        dscanfinal.append(peak)
        
dscanfinal = np.array(dscanfinal)

####################################################################
########### GENERATE LIST OF THEORETICAL IONS ######################
####################################################################

AAmass = 'A 71.037114,R 156.101111,N 114.042927,D 115.026943,C 103.009185,E 129.042593,Q 128.058578,G 57.021464,H 137.058912,I 113.084064,L 113.084064,K 128.094963,M 131.040485,F 147.068414,P 97.052764,S 87.032028,T 101.047679,W 186.079313,Y 163.06332,V 99.068414,'

iontable = []


masterpeaks = []

for pos in cxl1:
        
    masterpeaks.append('')
 
##########################################################################################    
##### BEGIN GIANT LOOP ###################################################################
##########################################################################################

k = 0

for pos in cxl1:

    #### gather b fragments without crosslink
    iontable = []
    
    posA = pos[0]
    posB = pos[1]

    i = 0
    temp = ''
    mz = 0


    if int(posA) > 1:

        for aa in list(peptideA):

            if i < (int(posA) - 1):

                nm = 0
                temp = temp + aa


                for let in list(temp):

                    nm = nm + float(''.join(re.findall(let + r'\s(\S*?)\,', AAmass)))

                j = 1
                mz = nm

                while mz >= 150 and mz <= 2000:

                    mz = (nm + (j * H))/float(j)

                    if mz >= float(150):

                        iontable.append([str(mz),'b',str(j),str(temp)])

                    j = j + 1

                i = i + 1

            else:

                break

    i = 0
    temp = ''
    mz = 0



    if int(posB) > 1:

        for aa in list(peptideB):

            if i < (int(posB) - 1):

                nm = 0
                temp = temp + aa


                for let in list(temp):

                    nm = nm + float(''.join(re.findall(let + r'\s(\S*?)\,', AAmass)))

                j = 1
                mz = nm

                while mz >= 150 and mz <= 2000:

                    mz = (nm + (j * H))/float(j)

                    if mz >= float(150):

                        iontable.append([str(mz),'b',str(j),str(temp)])

                    j = j + 1

                i = i + 1

            else:

                break

    ##### gather y fragments without crosslink            

    temp = ''

    if int(posA) < len(peptideA):

        pep = list(peptideA)
        del pep[0:posA]

        pep = pep[::-1]

        for aa in pep:

            nm = H2O
            temp = temp + aa
            label = list(temp)[::-1]


            for let in list(temp):

                nm = nm + float(''.join(re.findall(let + r'\s(\S*?)\,', AAmass)))

            j = 1
            mz = nm

            while mz >= 150 and mz <= 2000:

                mz = (nm + j * H)/float(j)

                if mz >= float(150):

                    iontable.append([str(mz),'y',str(j),str(''.join(label))])

                j = j + 1

    temp = ''

    if int(posB) < len(peptideB):

        pep = list(peptideB)
        del pep[0:posB]

        pep = pep[::-1]

        for aa in pep:

            nm = H2O
            temp = temp + aa
            label = list(temp)[::-1]


            for let in list(temp):

                nm = nm + float(''.join(re.findall(let + r'\s(\S*?)\,', AAmass)))

            j = 1
            mz = nm

            while mz >= 150 and mz <= 2000:

                mz = (nm + j * H)/float(j)

                if mz >= float(150):

                    iontable.append([str(mz),'y',str(j),str(''.join(label))])

                j = j + 1

    ##### ok, now time to gather fragments with crosslinks

    pepBmass = 0

    for let in list(peptideB):

        pepBmass = pepBmass + float(''.join(re.findall(let + r'\s(\S*?)\,', AAmass)))

    pepAmass = 0

    for let in list(peptideA):

        pepAmass = pepAmass + float(''.join(re.findall(let + r'\s(\S*?)\,', AAmass)))

    ################################################################
    ##### Frags with crosslinks#####################################
    ################################################################

    #### crosslink with pepA yfrags

    if posA > 1:

        pepAshmass = 0

        pepAsh = list(peptideA)[posA:len(peptideA)]

        for let in list(pepAsh):

            pepAshmass = pepAshmass + float(''.join(re.findall(let + r'\s(\S*?)\,', AAmass)))


        pep = list(peptideA)[0:posA]
        pep = pep[::-1]

        i = 0
        temp = '' 

        for aa in pep:

            nm = pepBmass + H2O + pepAshmass + H2O + linkmass  

            temp = temp + aa

            label = ''.join(list(temp)[::-1])

            for let in list(temp):

                nm = nm + float(''.join(re.findall(let + r'\s(\S*?)\,', AAmass)))

            j = 1
            mz = nm

            while mz >= 150:

                mz = (nm + j * H)/float(j)

                if mz >= float(150) and mz <= 2000:

                    iontable.append([str(mz),'y',str(j), str(''.join(label)) + str(''.join(pepAsh)) + ', ' + peptideB + '(' + str(i + 1) + ',' + str(posB)+ ')'])

                j = j + 1

            i = i + 1

    ##### crosslink with pepB yfrags

    temp = ''
    mz = 0



    if posB > 1:

        pepBshmass = 0

        pepBsh = list(peptideB)[posB:len(peptideB)]

        for let in list(pepBsh):

            pepBshmass = pepBshmass + float(''.join(re.findall(let + r'\s(\S*?)\,', AAmass)))

        pep = list(peptideB)[0:posB]
        pep = pep[::-1]   

        temp = ''
        i = 0

        for aa in pep:

            nm = pepAmass + H2O + pepBshmass + H2O + linkmass        
            temp = temp + aa

            label = ''.join(list(temp)[::-1])

            for let in list(temp):

                nm = nm + float(''.join(re.findall(let + r'\s(\S*?)\,', AAmass)))

            j = 1
            mz = nm

            while mz >= 150:

                mz = (nm + j * H)/float(j)

                if mz >= float(150) and mz <= 2000:

                    iontable.append([str(mz),'y',str(j), peptideA + ', ' + str(''.join(label)) + str(''.join(pepBsh))  + '(' + str(posA) + ',' + str(i + 1) + ')'])

                j = j + 1

            i = i + 1

    ##### crosslink with pepA bfrags


    temp = ''
    mz = 0

    if posA < len(peptideB):

        pepAfh = ''.join(list(peptideA)[0:posA])

        pepAsh = ['']
        pepAsh = pepAsh + list(peptideA)[posA:len(peptideA)-1]

        pepAfhmass = 0

        for let in list(pepAfh):

            pepAfhmass = pepAfhmass + float(''.join(re.findall(let + r'\s(\S*?)\,', AAmass)))

        for aa in pepAsh:

            nm = pepAfhmass + pepBmass + H2O + linkmass        
            temp = temp + aa

            for let in list(temp):

                nm = nm + float(''.join(re.findall(let + r'\s(\S*?)\,', AAmass)))

            j = 1
            mz = nm

            while mz >= 150:

                mz = (nm + j * H)/float(j)

                if mz >= float(150) and mz <= 2000:

                    iontable.append([str(mz),'b',str(j), pepAfh + str(''.join(temp)) + ', ' + peptideB + '(' + str(posA) + ',' + str(posB) + ')'])

                j = j + 1 


    ##### crosslink with pepB bfrags, seems good

    temp = ''
    mz = 0

    if posB < len(peptideB):

        pepBfh = ''.join(list(peptideB)[0:posB])

        pepBsh = ['']
        pepBsh = pepBsh + list(peptideB)[posB:len(peptideB)-1]

        pepBfhmass = 0

        for let in list(pepBfh):

            pepBfhmass = pepBfhmass + float(''.join(re.findall(let + r'\s(\S*?)\,', AAmass)))

        for aa in pepBsh:

            nm = pepBfhmass + pepAmass + H2O + linkmass        
            temp = temp + aa

            for let in list(temp):

                nm = nm + float(''.join(re.findall(let + r'\s(\S*?)\,', AAmass)))

            j = 1
            mz = nm

            while mz >= 150:

                mz = (nm + j * H)/float(j)

                if mz >= float(150) and mz <= 2000:

                    iontable.append([str(mz),'b',str(j), peptideA + ', ' + pepBfh + str(''.join(temp)) + '(' + str(posA) + ',' + str(posB) + ')'])

                j = j + 1

    #### REST OF MOD STUFF
    
    miontable = []
    
    if (countm + countc) != 0:
               
        for frag in iontable:
            
            temp = frag

            if re.findall('\(', temp[3]) != []:
                
                fs2 = ''.join(re.split(r'\(\d*\,\d*\)', temp[3]))
                ab = re.split(', ', fs2)
                fsfinal = ''.join(ab)
                
                r = re.findall(r'(.*)' + ab[0] + r'(.*)' + ab[1] + r'(.*)', sequence)
                
                remaining = ''.join(re.findall(r'\w*', str(r)))

            else:

                fsfinal = temp[3]
                remaining = ''.join(re.findall(r'(.*)' + fsfinal, sequence) + re.findall(fsfinal + r'(.*)', sequence))

            fragms = len(re.findall('M', fsfinal))
            fragcs = len(re.findall('C', fsfinal))
            
            otherms = len(re.findall('M', remaining))                
            unmodm = (fragms + otherms) - countm

            if fragms < countm:

                mopt = range(countm-otherms, fragms + 1)

            else: 

                mopt = mopt = range(countm-otherms, countm + 1)


            
            othercs = len(re.findall('C', remaining))
            unmodc = (fragcs + othercs) - countc

            if fragcs < countc:

                copt = range(countc-othercs, fragcs + 1)

            else: 

                copt = range(countc-othercs, countc + 1)
          

            ###############

            i = 0

            for item in mopt:

                if item <= 0:

                    mopt[i] = 0

                i = i + 1

            mopt = list(set(mopt))

            i = 0

            for item in copt:

                if item <= 0:

                    copt[i] = 0

                i = i + 1

            copt = list(set(copt))

            #################

            for m in mopt:

                for c in copt:

                    if m != 0:

                        fragsequ1 = re.sub('M', 'm', temp[3], count = m)

                    else:

                        fragsequ1 = temp[3]

                    if c != 0:

                        finalfragsequ = re.sub('C', 'c', fragsequ1, count = c)

                    else:

                        finalfragsequ = fragsequ1
                        
                    if (m + c) >0:
                        
                        newfragmass = float(temp[0])
                        
                    else:

                        newfragmass = float(temp[0])

                        newfragmass = ((float(c) * 57.021464 + float(m) * 15.99491)/(float(temp[2]))) + newfragmass

                    miontable.append([str(newfragmass), temp[1], temp[2], ''.join(finalfragsequ)])

    else:
        
        miontable = iontable
 
    ## now take -H2O and -NH2 into account
    
    fiontable = []
    
    for item in miontable:
        
        fiontable.append(item)
        
        fragmass = float(item[0])
        bory = item[1]
        chargestate = float(item[2])
        label = item[3]
        
        H2Ofragmass = fragmass + ((-H2O)/(float(item[2])))
        NH3fragmass = fragmass + ((-NH3)/(float(item[2])))
        
        fiontable.append([H2Ofragmass, bory, chargestate, label + '-H2O'])
        fiontable.append([NH3fragmass, bory, chargestate, label + '-NH3'])

        
    
    ### order by increasing m/z

    i = 0

    for item in fiontable:
        
        fiontable[i] = [float(item[0]), item[1], float(item[2]), item[3]]
        
        i = i + 1

    fiontable = sorted(fiontable, key=itemgetter(0), reverse = False)

    i = 0

    for item in fiontable:

        fiontable[i] = [float(item[0]), item[1], float(item[2]), item[3]]

        i = i + 1

    ###################################################################
    ########## MATCH SPECTRUM TO IONS #################################
    ###################################################################

    mz = 0
    candidates = []
    matched = []
    ppm = []

    for peak in dscanfinal:

        candidates = []

        ppm = []
        mz = float(peak[0])
        mzfind = ''.join(re.split('\.', str(peak[0]))[0])
        charge = int(peak[2])

        if charge == -2:

            for ion in fiontable:

                tempppm = 0            
                tempppm = ((abs(mz-float(ion[0])))/(float(ion[0]))*1000000)

                if tempppm < float(ppmlimit):

                    candidates.append(ion + [tempppm])

        else:

            for ion in fiontable:            

                tempppm = 0            
                tempppm = ((abs(mz-float(ion[0])))/(float(ion[0]))*1000000)

                if tempppm < float(ppmlimit) and charge == int(ion[2]):

                    candidates.append(ion + [tempppm])           

        candidates = sorted(candidates, key=itemgetter(4), reverse = False)

        if charge != -2:

            c=''

        else:

            c='u'

        if candidates != []:

            matched.append([mz, peak[1], candidates[0][4], str(candidates[0][1]) +'(+' + str(int(candidates[0][2])) + ')', candidates[0][3], c])

    #### remove duplicates

    matched = sorted(matched, key=itemgetter(3,4,5,2), reverse = False)
    
    finalmatched = []
    
    i = 0
    
    for row in matched:
        
        if i ==0:
            
            finalmatched.append(row)

        else:
            
            if i < (len(matched)-1):

                charge = row[3]
                sequ = row[4]

                prevcharge = matched[i-1][3]
                prevsequ = matched[i-1][4]

                if charge != prevcharge or sequ != prevsequ:

                    finalmatched.append(row)

        i = i + 1
 
    
    finalmatched = sorted(finalmatched, key=itemgetter(1), reverse = True)
    
    
    
    masterpeaks[k] = finalmatched
    
    #print str(posA) + ', ' + str(posB)
    
    #print tabulate(iontable)
    
    k = k + 1


#### GATHER TOTAL INTENSITY, AVG WEIGHT PPM, # PEAKS MATCHED

intensitysum = []
wavgppm = []
totpeaks = []
    
for item in masterpeaks:
    
    totintense = 0
    totpeaks.append(len(item))
    
    for entry in item:
        
        totintense = totintense + entry[1]
        
    intensitysum.append(totintense)
    
i = 0
    
for item in masterpeaks:
    
    tempwavgppm = 0
    
    for entry in item:
        
        tempwavgppm = tempwavgppm + entry[2]*(entry[1]/intensitysum[i])
        
    wavgppm.append(tempwavgppm)
        
    i = i + 1
    
intensitysumpos = []

i = 0

for x in intensitysum:
    
    intensitysumpos.append([x,i])
    
    i = i + 1
    
intensitysumpos = sorted(intensitysumpos, key=itemgetter(0), reverse = True)
        
bestmatchorder = []

for item in intensitysumpos:
    
    bestmatchorder.append(item[1])

mastertable = []
mastertable.append(['POSITION','PROT INDEX','TOT MATCHED INTENS','W AVG PPM','# PEAKS MATCHED'])


for x in bestmatchorder:

    mastertable.append([cxl1[x],cxl2[x], intensitysum[x], wavgppm[x], totpeaks[x]])

print 'scan number = ' + str(scantitle)
print 'total peak intensity = ' + str((np.transpose(dscanfinal))[1].sum())
print '\n'
print 'peptide A = ' + pepA + ', length = ' + str(len(pepA)) + ', prot index = ' + str(startA) + ', protein = ' + protA
print 'peptide B = ' + pepB + ', length = ' + str(len(pepB)) + ', prot index = ' + str(startB) + ', protein = ' + protB

consecutive = ''

if (startA + len(pepA)) == startB:
    
    consecutive = 'y'
    
elif (startB + len(pepB)) == startA:
    
    consecutive = 'y'
    
else:
    consecutive = 'n'
    
print 'adjacent? = ' + consecutive

print '\n'    
print tabulate(mastertable)
print '\n'

## now print matchedions


for x in bestmatchorder:
    
    temp = []
    
    print str(cxl1[x])
    
    temp.append(['M/Z', 'INTENSITY', 'PPM', 'ION', 'SEQUENCE', 'KNOWN CS'])
    
    temp = temp + masterpeaks[x][:numpeaksperlist]
    
    print tabulate(temp)
    print '\n'
    
#########################################################################
################### Visualize Spectrum ##################################
#########################################################################

if showspectra == 'y':
    
    seepeaks = bestmatchorder[:numspecs]
    
    for a in seepeaks:      
    
        fullscan = np.transpose(dscanfinal)

        figure = pylab.figure()

        ax = figure.add_subplot(1,1,1)

        ax.bar(fullscan[0], fullscan[1], width=0.1, facecolor='black')
    
        #creates figure title
        plt.title(str(cxl1[a]) + '  ' + str(cxl2[a]))

        ## sizes the plot appropriately
        matplotlib.pyplot.tight_layout()
        
        maxintense = max(fullscan[1])
        minintense = (float(minpercentintense)/100) * maxintense
        
        ### add labels

        mz = []
        intensity = []

        for item in masterpeaks[a]:

            mz.append(item[0])
            intensity.append(item[1])

        i = 0

        for x in mz:

            x = x
            y = intensity[i]
            y2 = maxintense * .1 + float(y)
            label = masterpeaks[a][i][3]+masterpeaks[a][i][5]


            if re.findall(r'\)', masterpeaks[a][i][4]) != []:

                label = label[:1]+'_'+label[1:]


            intensecheck = maxintense - maxintense * 0.2

            if showspectra == 'y':

                if float(y) >= intensecheck and float(y) >= minintense:

                    ax.annotate(label, xy=(float(x), float(y)), xytext=(float(x) + 10, float(y) - maxintense*0.1), arrowprops=dict(arrowstyle="->"), fontsize=8)

                else:

                    if float(y) >= minintense:

                        ax.annotate(label, xy=(float(x), float(y)), xytext=(float(x)-20, float(y2)), arrowprops=dict(arrowstyle="->"), fontsize=8)

            i = i + 1

print '\n'
print 'peaks in full scan = ' + str(len(scan))
print 'peaks in deisotoped scan = ' + str(len(dscanfinal))
print '\n'

#comparison = []

#match1 = []
#match3 = []

#for a in masterpeaks[bestmatchorder[0]]:
    
    #match1.append([a[0],a[1],a[2]])
    
#for a in masterpeaks[bestmatchorder[2]]:
    
    #match3.append([a[0],a[1],a[2]])

#for a in match1:
    
    #if a not in match3:
        
        #comparison.append(a)
        
#for a in match3:
    
    #if a not in match1:
        
        #comparison.append(a)
        
#print comparison
        



# In[ ]:




# In[53]:

matplotlib.pyplot.close("all")


# In[85]:




# In[ ]:



