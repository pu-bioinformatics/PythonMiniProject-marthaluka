#!/usr/bin/python/


"""
This module extracts the information from a PDB file and executes basic analysis on it. It considers PDB files that 
refer to 3D structures of proteins.
Martha Luka
"""

import os.path
import sys

# Option 1: Open a PDB file

def printOptions(my_file="None"):                            #print out options for the user
    "This provides the user with options. Its the first function to be called"
    print ("""
********************************************************************************
* PDB FILE ANALYZER                                                            *
********************************************************************************
* Select an option from below:                                                 *
*                                                                              *
* 1) Open a PDB File (O)                                                       *
* 2) Information (I)                                                           *
* 3) Show histogram of amino acids (H)                                         *
* 4) Display Secondary Structure (S)                                           *
* 5) Export PDB File (X)                                                       *
* 6) Exit (Q)                                                                  *
*                                                                              *
*                                                         Current PDB: %s  *
********************************************************************************
""" % (printTitle (my_file)))
    print("Please pick an option: ")
    option = input ()
    optionsMenu (option, my_file)



def printTitle (my_file):                                #get name of pdb file as saved by the user
    "Extract filename from the input path"
    
    sep = "/"
    if my_file=="None":
        fileName= "None"
    else:
        for i in my_file:
            i=i.replace(sep, " ")
        path=my_file.split("/")
        fileName = path[-1]
        return fileName
    
    
my_pdb_file=None
    
def askFileName (my_file="None"):                    #input path to file. Default is none
    "Allow user input of file and check whether its a valid file path"
    
    global my_pdb_file
    global path
    
    if my_file == "None":
        my_file=input ("Please input file name: ")
        path=str(my_file)
        if os.path.exists(my_file)==False:
            print ("Oops! We can't seem to find that file. Please try again")
            askFileName (my_file="None")
        else:
            my_pdb_file = open(my_file, "r")
            print ("The file %s has been successfully loaded" % printTitle (my_file))
            printOptions (my_file)
            
    else:
        loadAnotherFile (my_file)
                
        
        
def loadAnotherFile (my_file):                        #if user asks to load a file while another is already open
    "Has checks when the user wants to load another file"
    reply = input ("""
    A file is already open. 
    Would you like to replace the current file? 
    Please enter y for yes or n for no:
    """)                            
    if reply.lower() == "y":
        askFileName ()           
    else:
        printOptions (my_file)
        
# Option 2: Information

def printFile (my_pdb_file):                         #extract name of pdb file
    "Extract name of pdb file from inside the file, not based on how it is saved"
    while True:
        try:
            my_pdb_file.seek(0)
            for line in my_pdb_file: 
                if line.startswith('HEADER'):
                    print("PDB file: %s" %line[59:].strip())
                printTitleOfPdb ()
        except AttributeError:
            print("There's no file loaded. Load a file first")
            printOptions()
            break
           
        
def chains ():                         #chains in protein
    "Create a list of the chains in the proteins"
    global my_pdb_file
    while True:
        try: 
            my_pdb_file.seek(0)
            chainsSet = set()                                 ##to ensure only unique entries, use a set instead of a list
            for line in my_pdb_file:
                if line.startswith ("SEQRES"):
                    chainsSet.add(line[11])
            chainsList=list(chainsSet)                        #convert set to list
            chainsList.sort()                                 #sets are not ordered. This is to order the list created from the set
            return chainsList
        except AttributeError:
            print("No file loaded. kindly load a file first")
            printOptions()
            break

    
def printTitleOfPdb ():                           #Title of the pdb file
    "Title of the pdb"
    global my_pdb_file
    my_pdb_file.seek(0)
    empty_str =""
    for line in my_pdb_file:
        if line.startswith('TITLE'): 
            empty_str += line[5:].strip()
    print ("TITLE: %s" %empty_str)
    chainsInProtein=chains ()
    printChains(chainsInProtein)



def helicesPerChain (chainsInProtein):
    "Create a dictionary to count helices per chain"
    global my_pdb_file
    my_pdb_file.seek(0)
    helixDict ={}                              #initiating an empty dict for helices
    for i in chainsInProtein:
        helixDict[i]=0                         #default value is zero
    for line in my_pdb_file:
        if line.startswith('HELIX') and line[19] in chainsInProtein:
            helixDict[line[19]] += 1           #number of helices in the different chains of the protein
    return helixDict
    
    
def sheetsPerChain (chainsInProtein):
    "Create a dictionary to count sheets per chain"
    global my_pdb_file
    my_pdb_file.seek(0)
    sheetsDict ={}                             #initiating an empty dict for helices
    for i in chainsInProtein:
        sheetsDict[i]=0                         #default value is zero
    for line in my_pdb_file:
        if line.startswith('SHEET') and line[21] in chainsInProtein:
            sheetsDict[line[21]] += 1         #number of sheets in the different chains of the protein
    return sheetsDict
            

def AANumberPerChain (chainsInProtein):
    "Create a dictionary to count amino acid number per chain"
    global my_pdb_file
    my_pdb_file.seek(0)
    NumberOfAAPerChainList = []                           
    for line in my_pdb_file:
        if line.startswith ("SEQRES"):
            NumberOfAAPerChainList.append(line[12:17].strip())

    NumberOfAAPerChainList2=[]
    for i in NumberOfAAPerChainList:
        if i not in NumberOfAAPerChainList2:                           #unique numbers
            NumberOfAAPerChainList2.append(i)

    if len(NumberOfAAPerChainList2)<len(chainsInProtein):
        NumberOfAAPerChainList2.append(NumberOfAAPerChainList2[-1])
    my_dict=dict(zip(chainsInProtein,NumberOfAAPerChainList2))
    return my_dict

    
def AAsequencePerChain():
    "Create a dictionary to generate amino acid sequence per chain"
    global my_pdb_file
    chainsInProtein=chains()
    my_pdb_file.seek(0)
    aminoAcidsInChains ={}
    aaDict = {'CYS': 'C','ASP': 'D','SER': 'S','GLN': 'Q','LYS': 'K','ILE': 'I','PRO': 'P','THR': 'T','PHE': 'F','ASN': 'N','GLY': 'G','HIS': 'H','LEU': 'L','ARG': 'R','TRP': 'W','ALA': 'A','VAL': 'V','GLU': 'E','TYR': 'Y','MET': 'M'}
    for i in chainsInProtein:
        aminoAcidsInChains[i]=""

    for line in my_pdb_file:
        if line.startswith('SEQRES') and line[11] in chainsInProtein:
            aminoAcidsInChains[line[11]]+= (line[19:].strip())+" "   #extract 3 letter words from file

    aminoAcidsInChains2 ={}
    for i in chainsInProtein:                 
        a=aminoAcidsInChains.get(i)
        AAInChainList=a.split()             #create a list so as to apply the dictionary
        
        aminoAcidsInChains2[i]=""
        for aa in AAInChainList:            #extract one letter symbols using dict
            n=aaDict.get(aa)
            aminoAcidsInChains2[i]+="".join(n)   #Back to string
            
    for i in chainsInProtein:
        count=0
        multiLineSeq=""
        for aa in aminoAcidsInChains2[i]:
            count+=1
            multiLineSeq+=aa
            if count%50==0:
                multiLineSeq+="\n"
            aminoAcidsInChains2[i]=multiLineSeq      

    return (aminoAcidsInChains2)           

def printChains (chainsInProtein):
    "Sort out the information ready for and print"
    k=' and '.join(chainsInProtein)
    print ("CHAINS: %s" % k)
    chainsInProtein=chains ()
    
    global path
    aaNumber=AANumberPerChain (chainsInProtein)
    helixDict= helicesPerChain (chainsInProtein) 
    sheetsDict=sheetsPerChain (chainsInProtein)
    aminoAcidsInChains2= AAsequencePerChain()
    
    for chain in chainsInProtein:
        numberOfaa=aaNumber.get(chain)
        helicesInChain=helixDict.get(chain)
        sheetsInChain=sheetsDict.get(chain)
        aaSeq=aminoAcidsInChains2.get(chain)
        #printInfo(i,a,b,c,d)
        print ("""
    - CHAIN %s
     Number of amino acids: %s 
     Number of helix: %s
     Number of sheet: %s 
     Sequence: 
%s
                """ % (chain, numberOfaa, helicesInChain, sheetsInChain, aaSeq))
    printOptions(path)
 



#Option 3: Show histogram of amino acids

def getAllAminoAcids ():
    global my_pdb_file
    while True:
        try:
            aaPresentList=[]
            for line in my_pdb_file:
                if line.startswith('SEQRES'):
                    aaPresentList.append(line[19:].split())
            AminoAcidsInProtein=[]
            for aa in aaPresentList:                        #merge the list of lists into one list
                AminoAcidsInProtein+=aa

            ListOfPossibleAminoAcids = ['CYS','ASP','SER','GLN','LYS','ILE','PRO','THR','PHE','ASN','GLY','HIS','LEU','ARG','TRP','ALA','VAL','GLU','TYR','MET']
            aaFreqDict ={}
            for i in ListOfPossibleAminoAcids:              #create an empty dictionary using ListOfPossibleAminoAcids
                aaFreqDict[i]= (0)

            for i in AminoAcidsInProtein:
                if i in ListOfPossibleAminoAcids:          #count amino acids using dict
                    aaFreqDict[i]+=1
            return aaFreqDict
        except TypeError:
            print("No file loaded. Please load a file first")
            printOptions()
            break

def printOptions2 (aaFreqDict):                                        #print out options for the user
    "Prints out histogram options"
    
    print ("""
    Choose an option to order by:
         -number of amino acids - ascending (an)
         -number of amino acids - descending (dn)
         -alphabetically - ascending (aa)
         -alphabetically - descending (da)
    Order by:
     """)
    orderByOption = input ()
    optionsMenu2 (orderByOption, aaFreqDict)

    
def optionsMenu2 (orderByOption, aaFreqDict):                             #menu to handle the user input
    "Handles the user input"
    if orderByOption.lower() in ('an', 'dn', 'aa', 'da'):
        if orderByOption.lower() == 'an':                              
            numberAscending (aaFreqDict)
        elif orderByOption.lower() == 'dn':
            numberDescending(aaFreqDict)
        elif orderByOption.lower() == 'aa':
            alphabetAscending (aaFreqDict)
        elif orderByOption.lower() == 'da':
            alphabetDescending (aaFreqDict)             
    else:
        print ("Wrong choice! \nPlease try again")
        printOptions2 (aaFreqDict)

        
star= '*'

def numberAscending (aaFreqDict):
    "Arrange histogram based on ascending number of amino acids"
    global path
    for (key, value) in sorted(aaFreqDict.items(), key=lambda item:item[1]):
        if value>0:
            print ("%s (%3d): %s" %(key, value, star*int(value)))
            
    printOptions (path)

def numberDescending(aaFreqDict): 
    "Arrange histogram based on descending number of amino acids"
    global path
    for (key, value) in sorted(aaFreqDict.items(), key=lambda item:item[1], reverse=True):
        if value>0:
            print ("%s (%3d): %s" %(key, value, star*int(value)))
            
    printOptions (path)
    
def alphabetAscending (aaFreqDict):
    "Arrange histogram based on ascending alphabetical order"
    global path
    for key in sorted(aaFreqDict):
        if aaFreqDict[key] >0:
            print("%s (%3d): %s" %(key, aaFreqDict[key], star*int(aaFreqDict[key]))) 
            
    printOptions (path)

def alphabetDescending (aaFreqDict):
    "Arrange histogram based on descending alphabetical order"
    global path
    for key in sorted(aaFreqDict,reverse=True):
        if aaFreqDict[key] >0:
            print("%s (%3d): %s" %(key, aaFreqDict[key], star*int(aaFreqDict[key]))) 
            
    printOptions (path)
            


## Option 4: Display secondary structure

def multiToSingleLine ( aminoAcidsInChains2, chainsInProtein):
    "This converts the multiline aa sequence earlier generated (by option 2-information) to a single line aa sequence"
    
    AAInChainsSingleLine ={}
    for i in chainsInProtein:
        multiLine=aminoAcidsInChains2.get(i)
        singleLine=multiLine.replace('\n','')
        AAInChainsSingleLine[i]=singleLine       
    return AAInChainsSingleLine




def helixStructure (aaNumber, chainsInProtein):
    """ This function creates a dictionary for the chains in the protein and assigns helix symbols at respective 
    postions of the sec structures list"""
    
    global my_pdb_file
    my_pdb_file.seek(0)
    secStrListDict={}
    for chain in chainsInProtein:    
        secStrListDict[chain] = ['-']*int(aaNumber[chain])
        helixSymList=[]
        for line in my_pdb_file:
            if line.startswith("HELIX") and line[19] in chainsInProtein:
                helixSymList=["/"]*((int(line[32:37].strip())-int(line[20:27].strip()))+1)
                secStrListDict[line[19]][(int(line[20:27].strip())-1):int(line[32:37].strip())]=helixSymList
    
    return secStrListDict




def sheetStructure (aaNumber, chainsInProtein):
    """ This function creates a dictionary for the chains in the protein and assigns sheet symbols at respective 
    postions of the sec structures list"""
    
    global my_pdb_file
    secStrListDictHelices=helixStructure (aaNumber, chainsInProtein)
    my_pdb_file.seek(0)
    for chain in chainsInProtein:     
        sheetSymList=[]
        for line in my_pdb_file:
            if line.startswith("SHEET") and line[21] in chainsInProtein:
                sheetSymList=["|"]*((int(line[34:38].strip())-int(line[22:28].strip()))+1)
                secStrListDictHelices[line[21]][int(line[22:28].strip())-1:int(line[34:38].strip())]=sheetSymList
    
        secStrListDictHelices[chain]=''.join(secStrListDictHelices[chain])
    return secStrListDictHelices




def helixTag (chainsInProtein):
    "Give a suitable tag at the beginning of every helical secondary structure"
    
    global my_pdb_file
    aaNumber=AANumberPerChain (chainsInProtein)
    
    my_pdb_file.seek(0)
    tagListDict={}
    for chain in chainsInProtein: 
        tagListDict[chain]=[' ']*int(aaNumber[chain])
        tagList=[]
        for line in my_pdb_file:
            if line.startswith("HELIX") and line[19] in chainsInProtein:
                tag=line[11:15].strip()
                tagListDict[line[19]][int(line[20:27].strip())-1]=tag
    return tagListDict




def sheetTag (chainsInProtein, aaNumber):
    "Give a suitable tag at the beginning of every sheet secondary structure"
    
    global my_pdb_file
    #aaNumber=AANumberPerChain (chainsInProtein)
    
    tagListDict=helixTag (chainsInProtein)
    my_pdb_file.seek(0)
    for chain in chainsInProtein: 
        tagList=[]
       
        for line in my_pdb_file: 
            if line.startswith("SHEET") and line[21] in chainsInProtein:
                tag=''.join(line[5:15].split())
                tagListDict[line[21]][int(line[22:28].strip())-1:int(line[22:28].strip())-1+len(tag)]=tag
             
                
        tagListDict[chain]= ''.join(tagListDict[chain])
    return tagListDict





def fileName ():                           #Title of the pdb file
    "Extract the pdb id from the open file"
    global my_pdb_file
    my_pdb_file.seek(0)
    for line in my_pdb_file: 
        if line.startswith('HEADER'):
            return line[59:].strip()


        

def countLastPosition (chainsInProtein, aaNumber):
    "Create a label of the last aa position in every chain"
    lastPosDict={}
    for chain in chainsInProtein: 
        lastPosDict[chain]=[' ']*int(aaNumber[chain])
        lastPos='('+ aaNumber[chain]+')'
        lastPosDict[chain][-len(lastPos):-1]=lastPos
        lastPosDict[chain]=''.join(lastPosDict[chain])        
    return lastPosDict


def countFirstPosition (chainsInProtein, aaNumber):
    "Create a label of the first aa position in every chain"
    firstPosDict={}
    for chain in chainsInProtein: 
        firstPosDict[chain]=[' ']*int(aaNumber[chain])
        firstPos='(1)'
        firstPosDict[chain][0]=firstPos
        firstPosDict[chain]=''.join(firstPosDict[chain])
    return firstPosDict



def generateMultiLine (chainsInProtein, singleLine):
    """For ease of printing line for line, this function generates a multiline sequence from the necessary string objects
    (not more than 80)"""
    for chain in chainsInProtein:
        count=0
        multiLineSeq=""
        for i in singleLine[chain]:
            count+=1
            multiLineSeq+=i
            if count% 80==0:
                multiLineSeq+="\n"
            singleLine[chain]=multiLineSeq
            multiLine=singleLine
    
    return multiLine

            
def printSecStr ():
    """Print out the secondary structure"""
    global path
    chainsInProtein=chains ()
    aaNumber=AANumberPerChain (chainsInProtein)
    
    aminoAcidsInChains2= AAsequencePerChain()
    aminoAcidsSeqDict=multiToSingleLine (aminoAcidsInChains2, chainsInProtein)
    secStrListDictAll= sheetStructure (aaNumber, chainsInProtein)
    tagListDictAll = sheetTag(chainsInProtein, aaNumber)
    firstPosDict=countFirstPosition (chainsInProtein, aaNumber)
    lastPosDict=countLastPosition (chainsInProtein, aaNumber)
    
    
    aminoAcidsSeqDictMultiLine=generateMultiLine (chainsInProtein, aminoAcidsSeqDict)  #generate multiple lines of all the suitable dictionaries 
    secStrListDictAllMultiLine=generateMultiLine (chainsInProtein, secStrListDictAll)    #so as to print out the secondary structure at 80 AA per line
    tagListDictAllMultiLine=generateMultiLine (chainsInProtein, tagListDictAll)
    firstPosDictMultiLine=generateMultiLine (chainsInProtein, firstPosDict)
    lastPosDictMultiLine=generateMultiLine (chainsInProtein, lastPosDict)
    
    print("Secondary structure of the PDB ID: %s" %fileName ())
    for chain in chainsInProtein:
        print("Chain %s:" %chain)
        line1=aminoAcidsSeqDictMultiLine[chain].split('\n')
        line2=secStrListDictAllMultiLine[chain].split('\n')
        line3=tagListDictAllMultiLine[chain].split('\n')
        line4=firstPosDictMultiLine[chain].split('\n')
        line5=lastPosDictMultiLine[chain].split('\n')
        
        count=0
        for i in line1, line2, line3, line4, line5:
            while count < len(line1):
                print(line4[count])
                print(line1[count])
                print(line2[count])
                print(line3[count])
                print(line5[count])
                count+=1
    printOptions (path)



            
#Option 5: Export

def exportPDBFile ():
    "Export the pdb file"
    global my_pdb_file
    while True:
        try:
            my_pdb_file.seek(0)
            newFilePath=input("Please enter file path: ")
            with open(newFilePath, "w") as myNewFile: 
                for line in my_pdb_file:
                    myNewFile.writelines(line)
            print("The file %s has been succesfully exported" %newFilePath)
            printOptions ()
        except AttributeError:
            print ("There's no file loaded. Please load a file")
            printOptions()
            break
    

# Option 6: Exit

def exitFunction ():
    "Exit the module"
    global my_pdb_file
    choice = input("Do you want to exit(E) or do you want go back to the menu (M): ")
    if choice.lower() == "e":
        my_pdb_file.close()
        sys.exit()
    else:
        printOptions ()
        
        
#menu to handle the user input        
def optionsMenu (option, my_file):
    "Handles the user input in printOptions"
    if option.lower() in ('o', 'i', 'h', 's', 'x', 'q'):
        if option.lower() == 'o':                              #open pdb
            askFileName(my_file)
        elif option.lower() == 'i':
            printFile (my_pdb_file)
        elif option.lower() == 'h':
            aaFreqDict=getAllAminoAcids ()
            printOptions2 (aaFreqDict)
        elif option.lower() == 's':
            printSecStr ()
        elif option.lower() == 'x':
            exportPDBFile ()   
        elif option.lower() == 'q':                        
            exitFunction ()             
    else:
        print ("Wrong choice! \nPlease try again")
        printOptions(my_file)

        
def main():
    printOptions()

    
main()    #call the main function