#! /usr/bin/env python3

import sys, os, re, pickle, json,math, numpy as np
path=["/Users/jprieto/Docking/modules/pyproteinsExt/src", "/Users/jprieto/Docking/scripts/dockingPP"]
for way in path :
    if path not in sys.path :
        sys.path.append(way)
import pyproteinsExt.structure.coordinates as PDB
import pyproteinsExt.structure.operations as PDBop
from core_stats import ResStats, ContactStats , CmapRes, writeScores
from rotation_utils import trans_matrix, eulerFromMatrix
import ccmap

from multiprocessing import Pool

parserPDB = PDB.Parser()

class Pose(object):
    """Object containing a pose of a Docking prediction output and calculating its coordinates to
    fit to the reference receptor's docking position"""
    def __init__(self, belongsTo, id, euler, tr):
        self.id = id
        self.euler = euler
        self.translate = tuple([ (-1) * t * belongsTo.step for t in tr ])
        self._ccmap = None
        self.belongsTo = belongsTo
        self.ligOffset = tuple( [ (-1) * o for o in belongsTo.baryLIG] )
        self.recOffset = tuple( [ (-1) * o for o in belongsTo.baryREC] )
        self.dictorizedReceptor = None
        self.dictorizedLigand = None
        self.rmsd=None

    def __str__(self):
        return str(self.id) + ') ' + str(self.euler) + ' ' + str(self.translate)
    def __repr__(self):
        return str(self.id) + ') ' + str(self.euler) + ' ' + str(self.translate)

    def set_RMSD(self,RMSD):
        self.rmsd=RMSD

    def ccmap(self, dist=5):
        #print self._ccmap
        if not self._ccmap:
            pdbObjRec = self.belongsTo.pdbObjReceptor
            pdbObjLig = self.belongsTo.pdbObjLigand
            #tmp_ccmap = ccmap.zmap(dist, self.euler)
            self.dictorizedReceptor = pdbObjRec.atomDictorize
            self.dictorizedLigand = pdbObjLig.atomDictorize
            #print dist
            ccmapAsString = ccmap.zmap( (self.dictorizedReceptor, self.dictorizedLigand), dist,
                                   self.euler, self.translate, self.recOffset,  self.ligOffset)
            self._ccmap = json.loads(ccmapAsString)
        return self._ccmap

    def dump(self):
        pdbObjRec = self.belongsTo.pdbObjReceptor.clone()
        pdbObjLig = self.belongsTo.pdbObjLigand.clone()

        if self.dictorizedReceptor:
            pdbObjRec.setCoordinateFromDictorize(self.dictorizedReceptor)
        if self.dictorizedLigand:
            pdbObjLig.setCoordinateFromDictorize(self.dictorizedLigand)

        return str(pdbObjRec) + str(pdbObjLig)

    @property
    def resMapList(self):
        """Returns a list of all residues implied in the contact between this Pose and the receptor"""
        self.has_ccmap()
        residues=[]
        for contact in self._ccmap['data']:
            residues.append(CmapRes(contact['root'], role='Rec').index)
            for partner in contact['partners']:
                residues.append(CmapRes(partner,role='Lig').index)
        return list(set(residues))

    @property
    def contactMapList(self):
        contacts=[]
        for contact in self._ccmap['data']:
            for partner in contact['partners']:
                contacts.append((CmapRes(contact['root'], role='Rec').index,CmapRes(partner,role='Lig').index))
        return contacts

    @property
    def resSize(self):
        """Returns the number of residues in the contact map between this pose and the receptor """
        self.has_ccmap()
        return len(self.resMapList)

    @property
    def conSize(self):
        """Returns the number of contacts between this pose and the receptor """
        self.has_ccmap()
        consize=0
        for contact in self._ccmap['data']:
            for partner in contact['partners']:
                consize+=1
        return consize

    def SumScore(self, resStats, method='plain'):
        """ This function returns the sum of the residues statistics,
         method can take three values : 'plain' , 'freq' or 'log' """
        self.has_ccmap()
        if method == 'plain' :
            counts=resStats.plainResDict
            return sum([counts[res] if res in counts else 1 for res in self.resMapList])
        elif method == 'freq' :
            freqs=resStats.resFreq
            return sum([freqs[res] if res in freqs else 1/(len(self.belongsTo.pList)) for res in self.resMapList])
        elif method == 'log' :
            freqs=resStats.resFreq
            return sum([math.log(freqs[res]) if res in freqs else math.log(1/(len(self.belongsTo.pList))) for res in self.resMapList])
        else :
            raise Exception("Unknown '" + method + "' method value")

    def MeanScore(self, resStats, method='plain'):
        """ This function returns the mean of the residues statistics,
         method can take two values : 'plain' , 'freq'  """
        self.has_ccmap()
        if method == 'plain' :
            counts=resStats.plainResDict
            return sum([counts[res] if res in counts else 1 for res in self.resMapList])/self.resSize
        elif method == 'freq' :
            freqs=resStats.resFreq
            return sum([freqs[res] if res in freqs else 1/(len(self.belongsTo.pList)) for res in self.resMapList])/self.resSize
        else :
            raise Exception("Unknown '" + method + "' method value")

    def SquareSumScore(self, resStats, method='plain'):
        """ This function returns the square sum of the residues statistics,
        method can take three values : 'plain' , 'freq' or 'log' """

        self.has_ccmap()
        if method == 'plain' :
            counts=resStats.plainResDict
            return sum([counts[res] if res in counts else 1 for res in self.resMapList])
        elif method == 'freq' :
            freqs=resStats.resFreq
            return sum([freqs[res] if res in freqs else (1/(len(self.belongsTo.pList)))**2 for res in self.resMapList])
        elif method == 'log' :
            freqs=resStats.resFreq
            return sum([math.log(freqs[res]) if res in freqs else (math.log(1/(len(self.belongsTo.pList))))**2 for res in self.resMapList])
        else :
            raise Exception("Unknown '" + method + "' method value")

    def cmapSumScore(self,conStats, method= 'plain'):
        """ This function returns the sum of the contacts statistics,
         method can take three values : 'plain' , 'freq' or 'log' """
        self.has_ccmap()
        score=0
        if method == 'plain' :
            for contact in self._ccmap['data']:
                for partner in contact['partners']:
                    score+=conStats.get(CmapRes(contact['root'], role='Rec').index,CmapRes(partner,role='Lig').index)
            return score
        elif method == 'freq' :
            freqs=conStats.contactFreq
            for contact in self._ccmap['data']:
                for partner in contact['partners']:
                    score+=freqs.get(CmapRes(contact['root'], role='Rec').index,CmapRes(partner,role='Lig').index)
            return score
        elif method == 'log' :
            freqs=conStats.contactFreq
            for contact in self._ccmap['data']:
                for partner in contact['partners']:
                    score+=math.log(freqs.get(CmapRes(contact['root'], role='Rec').index,CmapRes(partner,role='Lig').index))
            return score
        else :
            raise Exception("Unknown '" + method + "' method value")

    def cmapMeanScore(self,conStats, method= 'plain'):
        """ This function returns the sum of the contacts statistics,
         method can take three values : 'plain' , 'freq'"""
        self.has_ccmap()
        score=0
        if method == 'plain' :
            for contact in self._ccmap['data']:
                for partner in contact['partners']:
                    score+=conStats.get(CmapRes(contact['root'], role='Rec').index,CmapRes(partner,role='Lig').index)
            return score/self.conSize
        elif method == 'freq' :
            for contact in self._ccmap['data']:
                for partner in contact['partners']:
                    score+=conStats.contactFreq.get(CmapRes(contact['root'], role='Rec').index,CmapRes(partner,role='Lig').index)
            return score/self.conSize
        else :
            raise Exception("Unknown '" + method + "' method value")

    def cmapSquareSumScore(self,conStats, method= 'plain') :
        """ This function returns the square sum of the contacts statistics,
        method can take three values : 'plain' , 'freq' or 'log' """
        self.has_ccmap()
        score=0
        if method == 'plain' :
            for contact in self._ccmap['data']:
                for partner in contact['partners']:
                    score+=conStats.get(CmapRes(contact['root'], role='Rec').index,CmapRes(partner,role='Lig').index)**2
            return score
        elif method == 'freq' :
            for contact in self._ccmap['data']:
                for partner in contact['partners']:
                    score+=conStats.contactFreq.get(CmapRes(contact['root'], role='Rec').index,CmapRes(partner,role='Lig').index)**2
            return score
        elif method == 'log' :
            for contact in self._ccmap['data']:
                for partner in contact['partners']:
                    score+=math.log(conStats.get(CmapRes(contact['root'], role='Rec').index,CmapRes(partner,role='Lig').index))**2
            return score
        else :
            raise Exception("Unknown '" + method + "' method value")

    def has_ccmap(self):
        """ Check ccmap has been calculated for this pose """
        if self._ccmap != None :
            return True
        else :
            raise Exception("You must calculate ccmap before calculating scores, use : DockData.ccmap(start=...,stop=...,dist=...)")

class DockData(object):

    def dump(self):
        return pickle.dumps(self)


    def __init__(self, **kwargs):
        self._name="complex"
        self.step = None
        self.nCells = None
        self.eulerREC = None
        self.fileREC = None
        self.fileLIG = None
        self.baryREC = None
        self.baryLIG = None
        self.pList = []
        self.pdbObjLigand = None
        self.pdbObjReceptor = None
        self.truePos=[]


    def setComplexName(self, name):
        self._name=name

    @property
    def complexName(self):
        return self._name

    #Post processing / statistics operations on set of poses
    def contactSummary(self):
        pass # access individual pose through self.pList

    def ccmap(self, **kwargs):
        start = kwargs['start'] if 'start' in kwargs else 0
        stop  = kwargs['stop']  if 'stop' in kwargs else len(self.pList)
        dist  = kwargs['dist']  if 'dist' in kwargs else 5

        #default distance value

        self.threadPacketSize =  kwargs['pSize'] if 'pSize' in kwargs else 200
        self.ncpu = kwargs['ncpu'] if 'ncpu' in kwargs else None

        elem = []
        indices= slice(start, stop, self.threadPacketSize).indices(len(self.pList))
        for i in range(*indices):
            j = i + self.threadPacketSize if i + self.threadPacketSize <= stop else stop
            elem.append( (self.pList[slice(i, j)], dist) )

        print("Created " + str(len(elem)) + " data packets (" +  str(self.threadPacketSize) + " zObjects each) for process pool")
        pool = Pool(self.ncpu) if self.ncpu else Pool() #note the default will use the optimal number of workers
        mpData = pool.map(mpCcmap, elem)
        pool.close()
        self._unpackProcessPool(mpData)

    def _unpackProcessPool(self, poolResults):
         #print self.pList[0]._ccmap
        #print hope[0][0]._ccmap
        print ('unpacking')
        i=0
        for pack in poolResults:
            for poseClone in pack:
                # Commented because creates error when start != 0
                # if poseClone.id != self.pList[i].id:
                #     raise ValueError('mismatching id ' + str(self.pList[i].id) + ' '  + str(poseClone.id) )
                self.pList[i]._ccmap = poseClone._ccmap
                self.pList[i].dictorizedReceptor = poseClone.dictorizedReceptor
                self.pList[i].dictorizedLigand = poseClone.dictorizedLigand
                i += 1

    def push(self, *args):
        pose=Pose(self, *args)
        self.pList.append( pose )
        return pose

    def setReceptor(self, pdbFile): # Optional chain arguments ?
        # Get receptor's structure
        self.pdbObjReceptor = parserPDB.load(file=pdbFile)

    def setLigand(self, pdbFile):
        # Get ligand's structure
        self.pdbObjLigand = parserPDB.load(file=pdbFile)

    @property
    def dictPos(self):
        """Dictorize position data for visualization"""
        panel_dict={'x':[],'y':[],'z':[],'a1':[],'a2':[],'a3':[]}
        for pose in self.pList :
            panel_dict['x'] .append(pose.translate[0])
            panel_dict['y'] .append(pose.translate[1])
            panel_dict['z'] .append(pose.translate[2])
            panel_dict['a1'] .append(pose.euler[0])
            panel_dict['a2'] .append(pose.euler[1])
            panel_dict['a3'] .append(pose.euler[2])
        return panel_dict

    @property
    def getStats(self):
        """Return both ResStats and ContactStats objects for the experiment """
        n=len(self.pList)
        con_stats=ContactStats(n,append=False, name=self.complexName)
        res_stats=ResStats(n, name=self.complexName)


        for i in range(n):
            residues=[]
            try:
                for cclist in self[i]._ccmap['data']:
                    rootRes = CmapRes(cclist['root'], role='Rec')
                    res_stats.addRes(rootRes)
                    residues.append(rootRes.index)

                    for partner in cclist['partners']:
                        res_stats[rootRes.index].increase_count(count='pond')
                        partnerRes = CmapRes(partner,role='Lig')
                        residues.append(partnerRes.index)
                        res_stats.addRes(partnerRes)

                        res_stats[partnerRes.index].increase_count(count='pond')
                        con_stats.incrMdTree(rootRes.index, partnerRes.index)
                for i in list(set(residues)):
                    res_stats[i].increase_count(count='plain')

            except TypeError:
                pass
            #     print(f'Warning : only {i} poses could be analysed' )
            #     con_stats.setSize(int(i))
            #     res_stats.setSize(int(i))
            #     break

        return (res_stats,con_stats)

    @property
    def contactStats(self):
        """Return ContactStats object for the experiment """
        n=len(self.pList)
        con_stats=ContactStats(n, append=False, name=self.complexName)

        for i in range(n):
            try:
                for cclist in self[i]._ccmap['data']:
                    rootRes = CmapRes(cclist['root'], role='Rec')
                    for partner in cclist['partners']:
                        partnerRes = CmapRes(partner,role='Lig')
                        con_stats.incrMdTree(rootRes.index, partnerRes.index)

            except TypeError:
                pass
            #     print(f'Warning : only {i} poses could be analysed ' )
            #     con_stats.setSize(int(i))
            #     break
        return con_stats

    @property
    def resStats(self):
        """Return ResStats object for the experiment """
        n=len(self.pList)
        res_stats=ResStats(n, name=self.complexName)

        for i in range(n):
            residues=[]
            try:
                for cclist in self[i]._ccmap['data']:
                    rootRes = CmapRes(cclist['root'], role='Rec')
                    res_stats.addRes(rootRes)
                    residues.append(rootRes.index)

                    for partner in cclist['partners']:
                        partnerRes = CmapRes(partner,role='Lig')
                        residues.append(partnerRes.index)
                        res_stats.addRes(partnerRes)
                        res_stats[rootRes.index].increase_count(count='pond')
                        res_stats[partnerRes.index].increase_count(count='pond')

                for i in list(set(residues)):
                    res_stats[i].increase_count(count='plain')

            except TypeError:
                pass
            #     print(f'Warning : only {i} poses could be analysed ' )
            #     res_stats.setSize(int(i))
            #     break

        return res_stats

    def bestPoses(self, stats=None, n=10, criteria = 'residue' , function='sum', method = 'freq') :
        """ criteria can be either 'residue'(default) or 'contact',
        function can take values : 'sum'(default),'square',
        method can take three values : 'plain' , 'freq'(default) or 'log' """
        if criteria =='residue':
            _functions={'sum' : Pose.SumScore , 'square' : Pose.SquareSumScore}
            if stats==None : stats=self.resStats
        elif criteria=='contact':
            _functions={'sum' : Pose.cmapSumScore, 'square' : Pose.cmapSquareSumScore }
            if stats==None : stats=self.contactStats
        else:
            raise Exception("Criteria must be 'residue' or 'contact'")
        size=stats.expSize
        return sorted(self.pList[:size], key=lambda o:_functions[function](o,stats , method=method))[:n]

    def poseScores(self, stats=None, criteria = 'residue' , function='sum', method = 'freq') :
        """ criteria can be either 'residue'(default) or 'contact',
        function can take values : 'sum'(default),'square','mean'
        method can take three values : 'plain' , 'freq'(default) or 'log' """
        if criteria =='residue':
            _functions={'sum' : Pose.SumScore , 'square' : Pose.SquareSumScore, 'mean': Pose.MeanScore}
            if stats==None : stats=self.resStats
        elif criteria=='contact':
            _functions={'sum' : Pose.cmapSumScore, 'square' : Pose.cmapSquareSumScore, 'mean':Pose.cmapMeanScore }
            if stats==None : stats=self.contactStats
        else:
            raise Exception("Criteria must be 'residue' or 'contact'")
        size=stats.expSize
        poses= { i:_functions[function](p,stats , method=method) for i,p in enumerate(self.pList[:size])}
        # sposes={i:poses[i] for i in sorted(poses.keys(), key=lambda o:poses[o], reverse=True)}
        # return sorted((self.pList[:size],_functions[function](self.pList[:size],stats , method=method)) , key=lambda o:_functions[function](o,stats , method=method))
        return poses

    def write_all_scores(self, size=1 , filename="scores", title='Exp1', header=None, F=False) :
        header = ["Surface size", "Residue freq sum", "Residue mean freq", "Residue log sum", "Residue square sum", "Number of contacts", "Contact freq sum", "Contact mean freq", "Contact log sum", "Contact square sum"]

        resS , conS, scores=self.all_scores()
        resS.write(filename+"_resstats.tab")
        conS.write(filename+"_constats.tab")
        assert len(list(set([len(i) for i in scores])))==1
        if F :
            e=False
        else :
            e=True
        score_file=filename + '.tsv'
        while e==True :
            if os.path.isfile(score_file):
                re= input(f"Warning : File {score_file} already exists, do you wish to continue anyway and replace it ? (yes/no)")
                if re=='yes':
                    os.system("rm " + score_file)
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
        with open(score_file,'a+') as f:
            f.write("# Title : " + title + "\n")
            f.write("# Experiment_size : " + str(size) + "\n")
            f.write("Pose" + "\t" + "\t".join(header)+ "\n")

            for pose in range(len(scores)):
                f.write(str(pose)+"\t"+ "\t".join([str(i) for i in scores[pose]]) + "\n")

        return score_file

    def all_scores(self) :
        scores= []
        resS , conS = self.getStats
        # resS.write(filename+"_resstats.tab")
        # conS.write(filename+"_constats.tab")
        size=resS.expSize
        rfreqs=resS.resFreq
        cfreqs=conS.contactFreq
        for i,p in enumerate(self.pList[:size]) :
            p.has_ccmap()
            complex=[0,0,0,0,0,0,0,0,0,0]
            for res in p.resMapList :
                complex[1] += rfreqs[res] if res in rfreqs else 1/size
                complex[3] += math.log(rfreqs[res]) if res in rfreqs else math.log(1/size)
                complex[4] += rfreqs[res]**2 if res in rfreqs else 1/size**2
            complex[0] = p.resSize
            try :
                complex[2] = complex[1]/p.resSize  # mean freq
            except ZeroDivisionError:
                if complex[1]==0:
                    complex[2]==0
                else :
                    raise Exception("weird pose : " + str(p.id) + " with size 0 resSize and "+ str(p.ccmap))
            for contact in p.contactMapList:
                # print(contact)
                complex[6] += cfreqs.get(contact[0],contact[1])
                complex[8] += math.log(cfreqs.get(contact[0],contact[1]))
                complex[9] += (cfreqs.get(contact[0],contact[1]))**2
            complex[5] = p.conSize
            complex[7] = complex[6]/p.conSize if p.conSize !=0 else 0# mean contact freq
            scores.append(complex)
        return resS,conS,scores

    def __str__(self):
        return str({ 'step' : self.step, 'nCells' : self.nCells , 'EulerREC' : self.eulerREC,
                       'fileREC' : self.fileREC, 'baryREC' : self.baryREC,
                       'fileLIG' : self.fileLIG, 'baryLIG' : self.baryLIG,
                        'pList'  : self.pList
                   })

    def __getitem__(self, key):
        return self.pList[key]

    def __iter__(self):
        for Pose in self.pList:
            yield Pose


def parse(fileName, maxPose = 0):
    reL1 = r'^([\d]+)[\s]+([\d\.]+)[\s]*$'
    reL2 = r'^[\s]*([\.\d-]+)[\s]+([\d\.-]+)[\s]+([\.\d-]+)[\s]*$'
    reL3 = r'^[\s]*([\S]+)[\s]+([\.\d-]+)[\s]+([\d\.-]+)[\s]+([\.\d-]+)[\s]*$'
    reZPOSE = r'^[\s]*[\d]+:[\s]+\([\s]*([\.\d]+),[\s]*([\.\d]+),[\s]*([\.\d]+)\)[\s]+\([\s]*([\.\d]+),[\s]*([\.\d]+),[\s]*([\.\d]+)\)$'

    reZPOSE = r'^[\s]*([\d-]+):[\s]+\([\s]*([\.\d-]+),[\s]*([\.\d-]+),[\s]*([\.\d-]+)\)[\s]+\([\s]*([\d-]+),[\s]*([\d-]+),[\s]*([\d-]+)\)'

    dockdataObj = DockData()
    dockBool = False
    R=False
    with open(fileName, 'r') as f:
        for line in f:
            if '*** docking results ***' in line:
                dockBool = True

            m = re.match(reL1, line)
            if m :
                dockdataObj.nCells = int(m.groups()[0])
                dockdataObj.step = float(m.groups()[1])
                continue
            m = re.match(reL2, line)
            if m :
                dockdataObj.eulerREC = ( float(m.groups()[0]), float(m.groups()[1]),float(m.groups()[2]) )
                continue
            m = re.match(reL3, line)
            if m :
                if not dockdataObj.fileREC:
                    dockdataObj.fileREC = m.groups()[0]
                    dockdataObj.baryREC = (float(m.groups()[1]), float(m.groups()[2]), float(m.groups()[3]))
                    continue
                dockdataObj.fileLIG = m.groups()[0]
                dockdataObj.baryLIG = (float(m.groups()[1]), float(m.groups()[2]), float(m.groups()[3]))
                continue
            m = re.match(reZPOSE, line)
            if m:
                if dockBool == False:
                    euler = (float(m.groups()[1]), float(m.groups()[2]), float(m.groups()[3]))
                    tr = [int(m.groups()[4]), int(m.groups()[5]), int(m.groups()[6])]
                    tr=tuple([ t - dockdataObj.nCells if t > dockdataObj.nCells / 2 else t for t in tr ])
                    dockdataObj.truePos.append(Pose(dockdataObj, int(m.groups()[0]), euler, tr))
                if dockBool == True:
                    try:
                        RMSD=float("".join(line.split(" ")[-8:-5]).strip())
                        R=True
                    except ValueError:
                        print("RMSD=" + "".join(line.split(" ")[-9:-4]).strip())
                        raise
                    # print(RMSD)
                    euler = (float(m.groups()[1]), float(m.groups()[2]), float(m.groups()[3]))


                    tr = [int(m.groups()[4]), int(m.groups()[5]), int(m.groups()[6])]
                    tr=tuple([ t - dockdataObj.nCells if t > dockdataObj.nCells / 2 else t for t in tr ])
                    pose=dockdataObj.push( int(m.groups()[0]), euler, tr)
                    if R :
                        pose.set_RMSD(RMSD)
                        R=False
                    if len(dockdataObj.pList) == maxPose:
                        return dockdataObj
        return dockdataObj

def zParse(fileName, maxPose = 0):
    """zDock pose-line format : 2.932153\t2.830209\t2.734943\t11\t5\t10\t1046.365"""
    reL1 = r'^([\d]+)[\s]+([\d\.]+)[\s]*$'
    reL2 = r'^[\s]*([\.\d-]+)[\s]+([\d\.-]+)[\s]+([\.\d-]+)[\s]*$'
    reL3 = r'^[\s]*([\S]+)[\s]+([\.\d-]+)[\s]+([\d\.-]+)[\s]+([\.\d-]+)[\s]*$'

    reZPOSE =  r'^([\d\.-]+)[\s]+([\d\.-]+)[\s]+([\d\.-]+)[\s]+([\d]+)[\s]+([\d]+)[\s]+([\d]+)[\s]+.*'

    dockdataObj = DockData()
    dockBool = False
    with open(fileName, 'r') as f:
        x=1
        for line in f:
            m = re.match(reL1, line)
            if m :
                dockdataObj.nCells = int(m.groups()[0])
                dockdataObj.step = float(m.groups()[1])
                continue
            m = re.match(reL2, line)
            if m :
                dockdataObj.eulerREC = ( float(m.groups()[0]), float(m.groups()[1]),float(m.groups()[2]) )
                continue
            m = re.match(reL3, line)
            if m :
                if not dockdataObj.fileREC:
                    dockdataObj.fileREC = m.groups()[0]
                    dockdataObj.baryREC = (float(m.groups()[1]), float(m.groups()[2]), float(m.groups()[3]))
                    continue
                dockdataObj.fileLIG = m.groups()[0]
                dockdataObj.baryLIG = (float(m.groups()[1]), float(m.groups()[2]), float(m.groups()[3]))
                continue
            m = re.match(reZPOSE, line)
            if m:
                euler = (float(m.groups()[0]), float(m.groups()[1]), float(m.groups()[2]))
                # Make rotation matrices
                rand_rot=trans_matrix(*dockdataObj.eulerREC)
                pose_rot=trans_matrix(*euler)
                # Combine into one matrix
                double=pose_rot.dot(rand_rot)
                # Recover combined angles
                euler=eulerFromMatrix(double)

                tr = [int(m.groups()[3]), int(m.groups()[4]), int(m.groups()[5])]
                tr=tuple([ t - dockdataObj.nCells if t > dockdataObj.nCells / 2 else t for t in tr ])
                dockdataObj.push(x, euler, tr)
                x+=1
                if len(dockdataObj.pList) == maxPose:
                    return dockdataObj
        return dockdataObj

def mpCcmap(datum):
    zPack, dist = datum
    for z in zPack:
            z.ccmap(dist=dist)
    return zPack
