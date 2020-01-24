#! /usr/bin/env python3

import sys, os, re, pickle, json,math, numpy as np

import pyproteinsExt.structure.coordinates as PDB
import pyproteinsExt.structure.operations as PDBop
from DockingPP.core_scores import Scores, multiPlot3D, countNative
from DockingPP.core_stats import ResStats, ContactStats , CmapRes, writeScores
from DockingPP.core_clustering import BSAS, birchCluster, wardCluster, herarCluster, ClusterColl
from DockingPP.rotation_utils import trans_matrix, eulerFromMatrix
import ccmap

from multiprocessing import Pool

parserPDB = PDB.Parser()

# self.columns={"original_rank":0,"r_size":1,"res_fr_sum":2,"res_mean_fr": 3, "res_log_sum": 4, "res_sq_sum": 5, "c_size": 6, "con_fr_sum": 7,"con_mean_fr" : 8,"con_log_sum" : 9,"con_sq_sum" : 10, "rmsd": 11}

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
        ## Check Receptor and ligand have been set
        self.has_rec()
        self.has_lig()

        if not self._ccmap:
            pdbObjRec = self.belongsTo.pdbObjReceptor
            pdbObjLig = self.belongsTo.pdbObjLigand
            #tmp_ccmap = ccmap.zmap(dist, self.euler)
            ligResCount=len(pdbObjLig.getResID)
            self.dictorizedReceptor = pdbObjRec.atomDictorize
            self.dictorizedLigand = pdbObjLig.atomDictorize
            #print dist
            pccmap = ccmap.zmap( (self.dictorizedReceptor, self.dictorizedLigand), dist,
                                   self.euler, self.translate, self.recOffset,  self.ligOffset)

            # Turn ccmap values into sparse matrix indexes : [(i,j) , ... ] euclidian division
            indexes=[(int(i/ligResCount), i% ligResCount)   for i in pccmap]
            self._ccmap = indexes
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
    def translatedCcmap(self) :
        """This function turns ccmap indexes into Contacts list [(recResidue, ligResidue), ... ]"""
        self.has_ccmap()
        def index(ResID):
            items=ResID.replace(":","").split()
            return ":".join(items)

        rec_residues=self.belongsTo.pdbObjReceptor.getResID
        lig_residues=self.belongsTo.pdbObjLigand.getResID
        tr_ccmap= [(index(rec_residues[i[0]]), index(lig_residues[i[1]])) for i in self._ccmap]
        return tr_ccmap


    @property
    def resMapList(self):
        """Returns a list of all residues implied in the contact between ligand and receptor in this Pose"""
        self.has_ccmap()
        residues=[]
        for contact in self.translatedCcmap:
            residues.append(CmapRes(contact[0], role='Rec').index)
            residues.append(CmapRes(contact[1],role='Lig').index)
        return list(set(residues))

    @property
    def contactMapList(self):
        """Returns a list of all contacts between ligand and receptor in this Pose"""
        self.has_ccmap ()
        contacts=[]
        for contact in self.translatedCcmap:
            contacts.append((CmapRes(contact[0], role='Rec').index,CmapRes(contact[1],role='Lig').index))
        return contacts

    @property
    def resSize(self):
        """Returns the number of residues in the contact map between ligand and receptor in this Pose """
        self.has_ccmap()
        return len(self.resMapList)

    @property
    def conSize(self):
        """Returns the number of contacts between ligand and receptor in this Pose"""
        self.has_ccmap()
        consize=0
        for contact in self.translatedCcmap:
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
            for contact in self.translatedCcmap:
                    score+=conStats.get(CmapRes(contact[0], role='Rec').index,CmapRes(contact[1],role='Lig').index)
            return score
        elif method == 'freq' :
            freqs=conStats.contactFreq
            for contact in self.translatedCcmap:
                score+=freqs.get(CmapRes(contact[0], role='Rec').index,CmapRes(contact[1],role='Lig').index)
            return score
        elif method == 'log' :
            freqs=conStats.contactFreq
            for contact in self.translatedCcmap:
                score+=math.log(freqs.get(CmapRes(contact[0], role='Rec').index,CmapRes(contact[1],role='Lig').index))
            return score
        else :
            raise Exception("Unknown '" + method + "' method value")

    def cmapMeanScore(self,conStats, method= 'plain'):
        """ This function returns the sum of the contacts statistics,
         method can take three values : 'plain' , 'freq'"""
        self.has_ccmap()
        score=0
        if method == 'plain' :
            for contact in self.translatedCcmap:
                score+=conStats.get(CmapRes(contact[0], role='Rec').index,CmapRes(contact[1],role='Lig').index)
            return score/self.conSize
        elif method == 'freq' :
            for contact in self.translatedCcmap :
                score+=conStats.contactFreq.get(CmapRes(contact[0], role='Rec').index,CmapRes(contact[1],role='Lig').index)
            return score/self.conSize
        else :
            raise Exception("Unknown '" + method + "' method value")

    def cmapSquareSumScore(self,conStats, method= 'plain') :
        """ This function returns the square sum of the contacts statistics,
        method can take three values : 'plain' , 'freq' or 'log' """
        self.has_ccmap()
        score=0
        if method == 'plain' :
            for contact in self.translatedCcmap:
                score+=conStats.get(CmapRes(contact[0], role='Rec').index,CmapRes(contact[1],role='Lig').index)**2
            return score
        elif method == 'freq' :
            for contact in self.translatedCcmap:
                score+=conStats.contactFreq.get(CmapRes(contact[0], role='Rec').index,CmapRes(contact[1],role='Lig').index)**2
            return score
        elif method == 'log' :
            for contact in self.translatedCcmap:
                score+=math.log(conStats.get(CmapRes(contact[0], role='Rec').index,CmapRes(contact[1],role='Lig').index))**2
            return score
        else :
            raise Exception("Unknown '" + method + "' method value")

    def has_ccmap(self, error=True):
        """ Check ccmap has been calculated for this pose """
        if self._ccmap != None :
            return True
        else :
            if error:
                raise Exception("You must calculate ccmap before using this function, use : DockData.ccmap(start=...,stop=...,dist=...)")
            else:
                return False

    def has_rec(self):
        """ Check Receptor has been set for this pose """
        if self.belongsTo.fileREC != None :
            return True
        else :
            raise Exception("You must set Receptor PDB before using this function, use : DockData.setReceptor(PDB_FILE)")

    def has_lig(self):
        """ Check Receptor has been set for this pose """
        if self.belongsTo.fileLIG != None :
            return True
        else :
            raise Exception("You must set Ligand PDB before using this function, use : DockData.setLigand(PDB_FILE)")

    @property
    def scores(self):
        if self.belongsTo:
            if self.belongsTo.scores:
                return self.belongsTo.scores[self.id]

class DockData(object):

    def dump(self):
        return pickle.dumps(self)


    def __init__(self, **kwargs):
        self._name="complex"
        self.step = None
        self.nCells = None
        self.eulerRand = None
        self.fileREC = None
        self.fileLIG = None
        self.baryREC = None
        self.baryLIG = None
        self.pList = []
        self.pdbObjLigand = None
        self.pdbObjReceptor = None
        self.truePos=[]
        self.scores=None

    def loadRMSD(self,filename=None):
        # Can be used for ZD or to use different RMSDS
        if not filename:
            raise Exception("You must set an rmsd file to use")
        RMSDS=[]
        with open(filename,'r') as f:
            for line in f.readlines():
                RMSDS.append(line.split('\t')[1].strip('\n'))
        for i in range(min(len(self.pList), len(RMSDS))):
            self.pList[i].set_RMSD(float(RMSDS[i]))
        return True


    def setScores(self, scores=None, filename=None):
        if filename:
            self.scores=Scores(filename=filename, data=None)
        elif scores:
            self.scores=Scores(data=scores, filename=None)
        else :
            try :
                self.scores=Scores(data=self.all_scores())
            except :
                return False
        self.scores.setPoses(self.pList)
        self.scores.belongsTo=self
        return True

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
            if not self[i].has_ccmap(error=False):
                print(f"{i} poses analysed")
                break
            residues=[]
            try:
                for contact in self[i].translatedCcmap:
                    rootRes = CmapRes(contact[0], role='Rec')
                    if rootRes.index not in residues:
                        res_stats.addRes(rootRes)
                        residues.append(rootRes.index)
                    res_stats[rootRes.index].increase_count(count='pond')

                    partnerRes = CmapRes(contact[1],role='Lig')
                    if partnerRes.index not in residues:
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
            if not self[i].has_ccmap(error=False):
                print(f'Warning : only {i} poses could be analysed ' )
                con_stats.setSize(int(i))
                break
            try:
                for contact in self[i].translatedCcmap:
                    rootRes = CmapRes(contact[0], role='Rec')
                    partnerRes = CmapRes(contact[1],role='Lig')
                    con_stats.incrMdTree(rootRes.index, partnerRes.index)

            except TypeError:
                pass
            #     con_stats.setSize(int(i))
            #     break
        return con_stats

    @property
    def resStats(self):
        """Return ResStats object for the experiment """
        n=len(self.pList)
        res_stats=ResStats(n, name=self.complexName)

        for i in range(n):
            if not self[i].has_ccmap(error=False):
                print(f'Warning : only {i} poses could be analysed ' )
                res_stats.setSize(int(i))
                break
            residues=[]
            try:
                for contact in self[i].translatedCcmap:
                    rootRes = CmapRes(contact[0], role='Rec')
                    if rootRes.index not in residues:
                        res_stats.addRes(rootRes)
                        residues.append(rootRes.index)
                    res_stats[rootRes.index].increase_count(count='pond')

                    partnerRes = CmapRes(contact[1],role='Lig')
                    if partnerRes.index not in residues:
                        residues.append(partnerRes.index)
                        res_stats.addRes(partnerRes)
                    res_stats[partnerRes.index].increase_count(count='pond')

                for i in list(set(residues)):
                    res_stats[i].increase_count(count='plain')

            except TypeError:
                pass
            #     print(f'Warning : only {i} poses could be analysed ' )
            #     break

        return res_stats

    def write_all_scores(self, size=1 , filename="scores", title='Exp1', header=None, F=False,resStats=None,conStats=None,maxPose=None) :
        header = ["Surface size", "Residue freq sum", "Residue mean freq", "Residue log sum", "Residue square sum", "Number of contacts", "Contact freq sum", "Contact mean freq", "Contact log sum", "Contact square sum"]

        if resStats and conStats:
            resS,conS=resStats,conStats
        else:
            resS , conS=self.getStats
        scores=self.all_scores(resStats=resS, conStats=conS,maxPose=maxPose)

        if maxPose:
            maxP=maxPose
        else:
            maxP=len(scores)

        resS.write(filename+"_resstats.tab", F=F)
        conS.write(filename+"_constats.tab",F=F)
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

            for pose in range(maxP):
                f.write("\t".join([str(i) for i in scores[pose]]) + "\n")
        self.setScores(scores=scores)
        return score_file

    def all_scores(self, resStats=None, conStats=None,maxPose=None) :
        header = ["Surface size", "Residue freq sum", "Residue mean freq", "Residue log sum", "Residue square sum", "Number of contacts", "Contact freq sum", "Contact mean freq", "Contact log sum", "Contact square sum"]
        scores= []
        if resStats and conStats :
            resS , conS = resStats, conStats
        else :
            resS , conS = self.getStats

        if maxPose:
            size=maxPose
        else:
            size=resS.expSize

        rfreqs=resS.resFreq
        cfreqs=conS.contactFreq
        for i,p in enumerate(self.pList[:size]) :
            if not p.has_ccmap(error=False):
                print(f"Warning : only {i} poses could be evaluated" )
                break
            complex=[None,0,0,0,0,0,0,0,0,0,0]
            complex[0]=p.id
            for res in p.resMapList :
                complex[2] += rfreqs[res] if res in rfreqs else 1/size
                complex[4] += math.log(rfreqs[res]) if res in rfreqs else math.log(1/size)
                complex[5] += rfreqs[res]**2 if res in rfreqs else 1/size**2
            complex[1] = p.resSize
            try :
                complex[3] = complex[2]/complex[1]  # mean freq
            except ZeroDivisionError:
                if complex[2]==0:
                    complex[3]==0
                else :
                    raise Exception("weird pose : " + str(p.id) + " with size 0 resSize and "+ str(p.ccmap))
            for contact in p.contactMapList:
                # print(contact)
                complex[7] += cfreqs.get(contact[0],contact[1])
                complex[9] += math.log(cfreqs.get(contact[0],contact[1]))
                complex[10] += (cfreqs.get(contact[0],contact[1]))**2
            complex[6] = p.conSize
            complex[8] = complex[7]/p.conSize if p.conSize !=0 else 0# mean contact freq
            scores.append(complex)
        self.setScores(scores=scores)
        return scores

    def has_scores(self):
        if self.scores :
            return True
        else :
            return False

    def has_rmsd(self):
        rmsd=True
        for p in self.pList :
            if not p.rmsd:
                rmsd=False
        return rmsd

    def __str__(self):
        return str({ 'step' : self.step, 'nCells' : self.nCells , 'eulerRand' : self.eulerRand,
                       'fileREC' : self.fileREC, 'baryREC' : self.baryREC,
                       'fileLIG' : self.fileLIG, 'baryLIG' : self.baryLIG,
                        'pList'  : self.pList
                   })

    def __getitem__(self, key):
        return self.pList[key]

    def __iter__(self):
        for Pose in self.pList:
            yield Pose

                ###############################################
                ##                                           ##
                ##  High level functions from scores object  ##
                ##                                           ##
                ###############################################

    def rankedPoses(self, element="original_rank", start=0, stop=None):
        if self.has_scores() or element=="original_rank":
            try :
                return self.scores.rankedPoses(element=element, start=start, stop=stop)
            except KeyError :
                raise Exception("pick a sorting element from \n \
'original_rank', 'r_size', 'res_fr_sum', 'res_mean_fr', 'res_log_sum', 'res_sq_sum', \n \
'c_size', 'con_fr_sum', 'con_mean_fr', 'con_log_sum', 'con_sq_sum'")
        else :
            raise Exception("You must set scores using self.setScores() or compute them with self.all_scores() before using rescoring functions")

    def rankedIDs(self, element="original_rank", start=0, stop=None):
        try :
            return [p.id for p in self.rankedPoses(element=element, start=start, stop=stop)]
        except KeyError :
            raise Exception("pick a sorting element from \n \
            'original_rank', 'r_size', 'res_fr_sum', 'res_mean_fr', 'res_log_sum', 'res_sq_sum', \
            'c_size', 'con_fr_sum', 'con_mean_fr', 'con_log_sum', 'con_sq_sum'")
        return [p.id for p in self.rankedPoses(element=element, start=start, stop=stop)]

    def ranks(self,element="original_rank"):
        if self.has_scores() or element=="original_rank":
            try :
                return self.scores.ranks(element=element)
            except KeyError :
                raise Exception("pick a sorting element from \n \
                'original_rank', 'r_size', 'res_fr_sum', 'res_mean_fr', 'res_log_sum', 'res_sq_sum', \
                'c_size', 'con_fr_sum', 'con_mean_fr', 'con_log_sum', 'con_sq_sum'")
        else :
            raise Exception("You must set scores before using rescoring functions")

    def rmsds(self):
        return self.rankedRmsds(element="original_rank", start=0, stop=None)

    def rankedRmsds(self,element="original_rank", start=0, stop=None):
        if (self.has_scores() or element=="original_rank") and self.has_rmsd() :
            try :
                rankedPoses=self.scores.rankedPoses(element=element, start=start, stop=stop)
                return self.scores.rankedRmsds(rankedPoses)
            except KeyError :
                raise Exception("pick a sorting element from \n \
                'original_rank', 'r_size', 'res_fr_sum', 'res_mean_fr', 'res_log_sum', 'res_sq_sum', \
                'c_size', 'con_fr_sum', 'con_mean_fr', 'con_log_sum', 'con_sq_sum'")
        else :
            raise Exception("You must compute scores ( use self.all_scores()) and set RMSDs (use loadRMSD(filename=FILENAME)) before using this function")

    def countNatives(self, element="original_rank", cutoff=5):
        if (self.has_scores() or element=="original_rank") and self.has_rmsd() :
            try :
                rankedPoses=self.scores.rankedPoses(element=element)
                rankedRmsds=self.scores.rankedRmsds(rankedPoses)
                return countNative(rankedRmsds, cutoff=cutoff)
            except KeyError :
                raise Exception("pick a sorting element from \n \
                'original_rank', 'r_size', 'res_fr_sum', 'res_mean_fr', 'res_log_sum', 'res_sq_sum', \
                'c_size', 'con_fr_sum', 'con_mean_fr', 'con_log_sum', 'con_sq_sum'")
        else :
            raise Exception("You must set scores and RMSDs before using this function")

    def plotFromPoses(self, rankedPoses, name='My complex', title='Docking decoys', size = (600,400)):
        """ Notebook version """
        self.scores.plot3D(rankedPoses,name=name,title=title, size = size)


    def plot3D(self, element="res_fr_sum", name='My complex', title='Docking decoys'):
        """ Notebook version """
        if self.has_scores() or element=="original_rank":
            try :
                rankedPoses=self.scores.rankedPoses(element=element)
                self.scores.plot3D(rankedPoses,name=name,title=title)
            except KeyError :
                raise Exception("pick a sorting element from \n \
                'original_rank', 'r_size', 'res_fr_sum', 'res_mean_fr', 'res_log_sum', 'res_sq_sum', \
                'c_size', 'con_fr_sum', 'con_mean_fr', 'con_log_sum', 'con_sq_sum'")
        else :
            raise Exception("You must set scores before using rescoring functions")

    def rmsdPlot(self, element="original_rank", start=0, stop=None, plot=None, title=None ):
        if self.has_scores() or element=="original_rank":
            try :
                rankedPoses= self.scores.rankedPoses(element=element)
                self.scores.rmsdGraphGenerator(rankedPoses, start=start, stop=stop, plot=plot, title=title )
            except KeyError :
                raise Exception("pick a sorting element from \n \
                'original_rank', 'r_size', 'res_fr_sum', 'res_mean_fr', 'res_log_sum', 'res_sq_sum', \
                'c_size', 'con_fr_sum', 'con_mean_fr', 'con_log_sum', 'con_sq_sum'")
        else :
            raise Exception("You must set scores before using rescoring functions")

    def multiPlot3D(self, wanted_scores,  title='Docking decoys', size=(600,400)):
        try :
            ranks=[self.scores.ranks(element=score) for score in wanted_scores]
            multiPlot3D([self.scores for i in wanted_scores], ranks, wanted_scores, title=title,size=size)
        except KeyError :
            raise Exception("pick a sorting element from \n \
            'original_rank', 'r_size', 'res_fr_sum', 'res_mean_fr', 'res_log_sum', 'res_sq_sum', \
            'c_size', 'con_fr_sum', 'con_mean_fr', 'con_log_sum', 'con_sq_sum'")


                ###############################################
                ##                                           ##
                ## High level functions from core_clustering ##
                ##                                           ##
                ###############################################

    def BSAS(self, maxd , element="original_rank", out="dict", start=0,stop=None):
        try :
            if out=='dict':
                clusters=ClusterColl(BSAS(self.rankedPoses(element=element), maxd, out=out, start=start,stop=stop ), DDObj=self)
            if out=='list' :
                clusters=BSAS(self.rankedPoses(element=element), maxd, out=out, start=start,stop=stop )
        except KeyError :
            raise Exception("pick a sorting element from \n \
            'original_rank', 'r_size', 'res_fr_sum', 'res_mean_fr', 'res_log_sum', 'res_sq_sum', \
            'c_size', 'con_fr_sum', 'con_mean_fr', 'con_log_sum', 'con_sq_sum', 'rmsd'")

        return clusters

    def wardCluster(self,maxd,  start=0, stop=None):
        clusters=ClusterColl(wardCluster(self, maxd, start=start, stop=stop ), DDObj=self)
        return clusters

    def birchCluster(self, maxd, out='dict', N=None):
        clusters=ClusterColl(birchCluster(self, maxd, out=out, N=N), DDObj=self)
        return clusters

    def herarCluster(self, maxc=None, linkage='complete', start=0, stop=None ):
        clusters=ClusterColl(herarCluster(self, maxc=maxc, linkage=linkage, start=start, stop=stop), DDObj=self)
        return clusters



        ###############################################################
        ##                                                           ##
        ##                                                           ##
        ##                      Parsing Functions                    ##
        ##                                                           ##
        ##                                                           ##
        ###############################################################

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
                dockdataObj.eulerRand = ( float(m.groups()[0]), float(m.groups()[1]),float(m.groups()[2]) )
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
                    if dockdataObject.eulerRand != (0,0,0) :
                        rand_rot=trans_matrix(*dockdataObj.eulerRand)
                        pose_rot=trans_matrix(*euler)
                        # Combine into one matrix
                        double=pose_rot.dot(rand_rot)
                        # Recover combined angles
                        euler=eulerFromMatrix(double)

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
                dockdataObj.eulerRand = ( float(m.groups()[0]), float(m.groups()[1]),float(m.groups()[2]) )
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
                # If two rotations are applied
                if dockdataObj.eulerRand != (0,0,0) :
                    # Make rotation matrices
                    rand_rot=trans_matrix(*dockdataObj.eulerRand)
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
