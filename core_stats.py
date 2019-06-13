import sys, os, re, pickle, json,math
import pyproteins.container.Core as core


class CmapRes(object):
    """ This class allows to store residues from ccmap and their counts ( in contacts and plain )"""
    def __init__(self, datum, role=None):
        """ Takes datum : {'chainID' : , 'resID' :  }"""
        self.chainID = datum['chainID'].strip()
        self.resID = datum['resID'].strip()
        self.role = role
        self.count=0
        self.cCount=0

    @property
    def index(self):
        """This function compounds an index or ID string to target the residue
        example : 147HRec , residue's id is 147, from chain H in receptor protein """
        return str(str(self.resID) +self.chainID + self.role)

    def __hash__(self):
        return hash(self.chainID + str(self.resID) + self.role)

    def __str__(self):
        return self.resID + self.chainID + self.role

    def __repr__(self):
        return self.__str__()

    def __eq__(self, other):
        try:
            if self.index == other.index :
                return True
            else :
                return False
        except AttributeError :
            return False

    def increase_count(self, count='plain'):
        # Use count='plain' for residue occurence in poses
        # Use count='pond' for residue occurence in contacts
        if count=='plain':
            self.count+=1
        elif count=='pond':
            self.cCount+=1
        else :
            raise Exception("Unknown count method '" + count + "' for CmapRes.increase_count")
    def reset_all(self):
        self.count=0
        self.cCount=0

class ResStats(object):
    """This object allows the storage and transformation or the statistics on the residues of the contact maps
    in a zDock or MegaDock experiment"""
    def __init__(self,expSize,name=None):
        # dictionary of zDock objects ex : {"154HLig" : CmapResobj}
        self.rDict={}
        try :
            self.expSize=int(expSize)
        except ValueError:
            raise Exception(" The experiment size value must be an integer ")
        self.expName=name

    def setSize(self,n):
        self.expSize=n

    def addRes(self, resObj):
        if resObj.index not in self.rDict:
            self.rDict[resObj.index]=resObj

    def __return__(self):
        return self.rDict

    def __getitem__(self,key):
        try:
            return self.rDict[key]
        except KeyError :
            raise Exception('You Must add resObj ' + str(key) + 'to resStats' )

    @property
    def plainResDict(self):
        """ Returns a dictionnary {"154HLig" : 3} of plain counts.
        Ocurrence number of the residue in the contact maps """
        count={}
        for res in self.rDict :
            count[res]=self.rDict[res].count
        return count

    @property
    def resFreq(self):
        """ Returns a dictionnary {"154HLig" : 3} of frequences.
        Ocurrence number of the residue in the contact maps divided by the size of the experiment (number of complexes) """
        normd={}
        for key in self.rDict :
            _v=self.rDict[key].count
            normd[key]=_v/self.expSize
        return normd

    @property
    def pondResDict(self):
        """ Returns a dictionnary {"154HLig" : 3} of counts.
        Ocurrence number of the residue in all the contacts ( can appear more than once in one only contact map ) """
        count={}
        for res in self.rDict :
            count[res]=self.rDict[res].cCount
        return count

    def write(self,filename, type='freq'):
        """Export the Residues frequencies out to a tab delimited file """
        if type=='freq' :
            data="Residue Frequencies"
            dataset = self.resFreq
        elif type=='pond':
            data="Residue Contact Counts"
            dataset = self.pondResDict
        if self.expName :
            data+=f"({self.expName})"

        writeScore(dataset, size=self.expSize, filename=filename, title=data )

class ContactStats(core.mdTree):
    """ This class allows to store contacts counts from ccmap """
    def __init__(self,expSize,append=False,name=None):
        super().__init__(append=append)
        try :
            self.expSize=int(expSize)
        except ValueError :
            raise Exception("The experiment size must be an integer ")
        self.expName=name
    def incrMdTree(self, x, y):
        _v = self.get(x, y) if self.get(x, y) else 0
        _v += 1
        self.set(x, y, _v)

    def setSize(self,n):
        self.expSize=n

    def get(self, k1, k2): # Mutable ref returned
        x, y = self._digest(k1, k2)
        if x not in self.data:
            return 1/self.expSize
        if y not in self.data[x]:
            return 1/self.expSize
        return self.data[x][y]

    @property
    def contactFreq(self):
        """Returns normalised contact counts by the experiment size (number of complexes) """
        newCounts=ContactStats(self.expSize, name=self.expName)
        for x in self.data:
            for y in self.data[x]:
                _v = self.get(x, y)
                newCounts.set(x, y, _v/self.expSize )
        return newCounts

    @property
    def logFreq(self):
        """Returns a ContactStats object carrying the log of the frequencies of each contact"""
        logCounts=ContactStats(self.expSize)
        for x in self.data:
            for y in self.data[x]:
                _v = self.get(x, y)
                logCounts.set(x, y, math.log(_v/self.expSize))
        return logCounts

    def render_table(self, n=None):
        """Returns a table of size n**2 with counts in the intersection of columns and lines"""
        L,R=[],[]
        for x in self.data:
            if x[-3:]=='Lig':
                L.append(x)
                continue
            if x[-3:]=='Rec':
                R.append(x)
        a='\t'+'\t'.join([i for i in R[:n]])
        print(a)
        for l in range(len(L[:n])):
            b=L[l]+'\t'+'\t'.join([str(self.get(L[l],R[i])) if self.get(L[l],R[i])!=None else str(0) for i in range(len(R[:n]))])
            print(b)

    @property
    def all(self):
        """ Returns all contact counts in the dataset allowing to perform statistical analyses and overall view of the counts distribution """
        all_val=[]
        for x in self.data:
            for y in self.data[x]:
                all_val.append(self.get(x, y))
        return all_val

    def write(self,filename):
        """Export the Contact frequencies out to a tab delimited file """
        scores={}
        freq= self.contactFreq
        title="Contact Frequencies"
        if self.expName :
            title+=f"({self.expName})"
        for x in self.data :
            for y in self.data[x] :
                scores[x+"_"+y]=freq.get(x,y)
        writeScore(scores, size=self.expSize, filename=filename, title=title )


def writeScore(scores , size=1, filename="scores.tsv", title='Exp1'):
    e=True
    while e==True :
        if os.path.isfile(filename):
            re= input(f"Warning : File {filename} already exists, do you wish to continue anyway and replace it ? (yes/no)")
            if re=='yes':
                os.system("rm " + filename)
                e=False
                break
            elif re=='no':
                filename=input("new file name (or 'x' to exit):  ")
                if filename=="x" :
                    return
            else :
                print("Please answer 'yes' or 'no' ")
        else :
            e=False
    with open(filename,'a+') as f:
        f.write("# Title : " + title + "\n")
        f.write("# Experiment_size : " + str(size) + "\n")
        for i in scores :
            f.write(str(i)+ "\t" + str(scores[i]) + "\n")

    return filename

def writeScores(scores , size=1, filename="scores.tsv", title='Exp1', header=None):
    e=True
    if header==None :
        header=["score" + str(i) for i in range(len(scores[0]))]
    else:
        assert len(header)==len(scores[0])
    assert len(list(set([len(i) for i in scores])))==1
    while e==True :
        if os.path.isfile(filename):
            re= input(f"Warning : File {filename} already exists, do you wish to continue anyway and replace it ? (yes/no)")
            if re=='yes':
                os.system("rm " + filename)
                e=False
                break
            elif re=='no':
                filename=input("new file name (or 'x' to exit):  ")
                if filename=="x" :
                    return
            else :
                print("Please answer 'yes' or 'no' ")
        else :
            e=False
    with open(filename,'a+') as f:
        f.write("# Title : " + title + "\n")
        f.write("# Experiment_size : " + str(size) + "\n")
        f.write("Pose" + "\t" + "\t".join(header)+ "\n")

        for pose in range(len(scores)):
            f.write(str(pose)+"\t"+ "\t".join([str(i) for i in scores[pose]]) + "\n")

    return filename
