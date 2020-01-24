#! /usr/bin/env python3

import matplotlib.pyplot as plt
import math
import numpy as npy
import pickle
from math import sqrt
# Import dependencies
import plotly
import plotly.graph_objs as go
from mpl_toolkits.mplot3d import Axes3D
import mpld3


class Scores(object):
    """This class allows to read a rescoring file"""
    def __init__(self, filename=None, data=None):
        self.data=None
        if filename :
            self.loadFile(filename)
        elif data:
            self.loadScores(data)
        self.columns={"original_rank":0,"r_size":1,"res_fr_sum":2,"res_mean_fr": 3, "res_log_sum": 4, "res_sq_sum": 5, "c_size": 6, "con_fr_sum": 7,"con_mean_fr" : 8,"con_log_sum" : 9,"con_sq_sum" : 10}
        self.belongsTo=None
        self.poses=None

    def __repr__(self):
        return str(self.data)

    def __str__(self):
        return str(self.data)

    def loadFile(self,file):
        data=[]
        with open(file,'r') as f:
            for line in f.readlines()[3:]:
                pose=[]
                for info in line.strip("\n").split("\t"):
                    pose.append(float(info))
                data.append(pose)
        self.data=data
        return True

    def loadScores(self,scores):
        self.data=scores

    def setPoses(self, pList):
        self.poses=pList

    def __getitem__(self,poseID):
        dict_for_pose={}
        for key in self.columns:
            try :
                dict_for_pose[key] = self.data[poseID-1][self.columns[key]]
            except IndexError :
                break
        return dict_for_pose

    @property
    def getPoses(self):
        # Poses defined with setPoses will be preferentially used
        if self.poses:
            return self.poses
        if self.belongsTo :
            return self.belongsTo.pList
        else :
            raise Exception("Please define poses using setPoses()")

    def getScore(self,poseid, score):
        return self.scoresDict[int(poseid)-1][self.columns[score]]

    @property
    def scoresDict(self):
        scoresdict={}
        for info in self.data:
            scoresdict[int(info[0])]=info[1:]
        return scoresdict

    @property
    def scoresList(self):
        return self.data

    @property
    def rmsds(self):
        """ return Rmsds in the original order """
        if self.getPoses :
            return [pose.rmsd for pose in self.getPoses]
        else :
            raise Exception("Please define poses using setPoses")


    def coordDict(self, start=0, stop=None):
        """dockingPP poses position dictorizer """
        coord={'x':[],'y':[],'z':[],'a1':[],'a2':[],'a3':[]}
        if not stop:
            stop=len(self.getPoses)
        for pose in self.getPoses[int(start):int(stop)]:
            coord['x'].append(pose.translate[0])
            coord['y'] .append(pose.translate[1])
            coord['z'] .append(pose.translate[2])
            coord['a1'] .append(pose.euler[0])
            coord['a2'] .append(pose.euler[1])
            coord['a3'] .append(pose.euler[2])
        return coord



    ###############################################
    ##                                           ##
    ##            Rescoring Functions            ##
    ##                                           ##
    ###############################################


    def rankedPoses(self, element="res_mean_fr", start=0, stop=None):
        """Returns poses' ids,  sorted according to rescoring element : ]
        "original_rank","s_size","res_fr_sum","res_mean_fr","res_log_sum","res_sq_sum",
        "c_size","con_fr_sum","con_mean_fr","con_log_sum","con_sq_sum"
        Means the pose classified 1st is 1st position of the list.   """
        try :
            p=self.getPoses
        except :
            raise Exception("Please define poses using setPoses() function")
        if not stop:
            stop=len(self.data)
        col=self.columns[element]
        if element == "original_rank" :
            r=sorted(self.data[int(start):int(stop)],key=lambda o:float(o[col]))
        else:
            r=sorted(self.data[int(start):int(stop)],key=lambda o:float(o[col]),reverse=True)

        sorted_i=[self.getPoses[int(pose[0])-1] for pose in r] # get pose from ID with getPoses function
        return sorted_i



    def rankedRmsds(self, rankedPoses):
        rmsds=[]
        for i in rankedPoses:
            rmsds.append(i.rmsd)
        return rmsds



    def ranks(self, element="res_mean_fr", start=0, stop=None):
        """Returns ranks for poses in the original order
        At position 0 we have pose 0 with rank 37 : rank[0] = 37
        Means the pose in original rank 0 is in 37th position according to rescoring
        element:
        "original_rank","s_size","res_fr_sum","res_mean_fr","res_log_sum","res_sq_sum",
        "c_size","con_fr_sum","con_mean_fr","con_log_sum","con_sq_sum"  """

        return self.ranksFromRankedPoses(self.rankedPoses(element=element, start=start, stop=stop))

    def ranksFromRankedPoses(self,rankedPoses):
        rank=[t[0]+1 for t in sorted([(u,v.id) for u,v in enumerate(rankedPoses)], key=lambda o:o[1])]
        return rank


    ###############################################
    ##                                           ##
    ##         Visualization Functions          ##
    ##                                           ##
    ###############################################

    def rmsdGraphGenerator(self, rankedPoses, start=0, stop=None, plot=None, title=None ):
        pl=plot if plot else plt
        if title:
            try :
                pl.set_title(title)
            except AttributeError :
                pl.title(title)
        try:
            pl.xlabel("Rank")
            pl.ylabel("RMSD")
        except AttributeError :
            pl.set_xlabel("Rank")
            pl.set_ylabel("RMSD")
        x=0
        if not stop:
            stop=len(rankedPoses)
        rmsds=self.rankedRmsds(rankedPoses)
        colors=colorsFromRmsd(rmsds)[int(start):int(stop)]
        x=[i if i!=0 else 0 for i in range(len(colors))]
        y=rmsds[int(start):int(stop)]
        # print(len(x))
        # print(len(y))
        pl.scatter(x, rmsds[int(start):int(stop)], c=colors)

    def histRmsd(self,plot=None, stop=None):
        pl= plot if plot else plt
        s=stop if stop else len(self.data)
        pl.hist(sorted([i[1] for i in self.data[:stop]]))

    def trace(self, rankedPoses, name):
        colorscale=make_colorScale()

        x = [ p.translate[0] for p in rankedPoses]
        y = [ p.translate[1] for p in rankedPoses]
        z = [ p.translate[2] for p in rankedPoses]
        # pos=self.coordDict(start=0, stop=len(rank)-1)
        C=[i+1 for i in range(len(rankedPoses))]
        # print(C)
        S=[ 15 if p.rmsd > 5 else 30 for p in rankedPoses ]
        # S=[ 15]
        # Configure the trace.
        trace = go.Scatter3d(
            x=x, y=y, z=z,
            mode = 'markers',
            marker={
                'size': S,
                'opacity': 1,
                'color' : C,
                'colorscale' : colorscale,
                'colorbar': dict(title = 'Ranks',
            lenmode='fraction',
            len=0.75,
            thickness=20,
            titleside = 'top',
            tickmode = 'array',
            x=0.99)},
            name=name,
            text=['id=' + str(p.id) + ', rank=' + str(u+1) + ', rmsd='+ str(p.rmsd) for u,p in enumerate(rankedPoses)])

        return trace

    def plot3D(self, rankedPoses, name=' ', title='Docking decoys',  size = (530,480)):
        # ranks=self.ranksFromRankedPoses(rankedPoses)
        width,height=size
        """ Suited for notebook display """
        """ Warning : ranks must be in the original order, since positions and rmsds are in the original order """
        # Configure Plotly to be rendered inline in the notebook.
        title_text=title
        plotly.offline.init_notebook_mode()
        trace=self.trace(rankedPoses, name)
        layout = go.Layout(autosize=False,
        width=width,
        height=height,
        title=go.layout.Title(text=title_text, y=0.98 ),
        hovermode= 'closest',
        hoverlabel=dict(bgcolor='white', bordercolor='black'),
        margin={'l': 0, 'r': 0, 'b': 0, 't': 10})
        data = [trace]
        plot_figure = go.Figure(data=data, layout=layout)
        # Render the plot.
        plotly.offline.iplot(plot_figure,filename='Docking Decoys')
        return plot_figure

    # NB : This function should take different sets of points and not sc object. --> put to core_visuals ?

def multiPlot3D(sc, ranks, names, title='Docking decoys', size = (530,480)):
    """ Warning : ranks must be in the original order, since positions and rmsds are in the original order
    ranks is a list of lists of ranks"""
    # Configure Plotly to be rendered inline in the notebook.
    title_text=title
    width,height=size
    plotly.offline.init_notebook_mode()
    traces=[sc[i].trace(ranks[i],names[i]) for i in range(len(ranks))]
    layout = go.Layout(
    # margin=go.layout.Margin(
    #     l=50,
    #     r=50,
    #     b=100,
    #     t=100,
    #     pad=4 ),
    autosize=False,
    width=width,
    height=height,
    title=go.layout.Title(text=title_text, y=0.98 ),
    hovermode= 'closest'),
    margin={'l': 0, 'r': 15, 'b': 0, 't': 10}
    data = traces
    plot_figure = go.Figure(data=data, layout=layout)
    # Render the plot.
    plotly.offline.iplot(plot_figure,filename='Docking Decoys')

    return plot_figure


def make_colorScale():
    colorscale= [
    # Let first 10% (0.1) of the values have color rgb(0, 0, 0)
    # [0, 'rgb(400, 180, 180)'],
    # [0.1, 'rgb(400, 180, 180)'],
    [0, 'rgb(400, 0, 0)'],
    [0.05, 'rgb(400, 0, 0)'],
    [0.05, 'rgb(400, 100, 40)'],
    [0.1, 'rgb(400, 100, 40)'],
    # Let values between 10-20% of the min and max of z
    # have color rgb(20, 20, 20)
    [0.1, 'rgb(400, 120, 60)'],
    [0.2, 'rgb(400, 120, 60)'],

    # Values between 20-30% of the min and max of z
    # have color rgb(40, 40, 40)
    [0.2, 'rgb(400, 100, 80)'],
    [0.3, 'rgb(400, 100, 80)'],

    [0.3, 'rgb(400, 80, 100)'],
    [0.4, 'rgb(400, 80, 100)'],

    [0.4, 'rgb(300, 80, 120)'],
    [0.5, 'rgb(300, 80, 120)'],

    [0.5, 'rgb(250, 80, 200)'],
    [0.6, 'rgb(250, 80, 200)'],

    [0.6, 'rgb(200, 50, 250)'],
    [0.7, 'rgb(200, 50, 250)'],

    [0.7, 'rgb(150, 50, 300)'],
    [0.8, 'rgb(150, 50, 300)'],

    [0.8, 'rgb(100, 0, 350)'],
    [0.9, 'rgb(100, 0, 350)'],

    [0.9, 'rgb(0, 0, 400)'],
    [1.0, 'rgb(0, 0, 400)'] ]
    return colorscale

def true3D(pos, complex):
    """Make sure the first 7 positions of your dictionnary pos are Native Like Solutions"""
    # Configure Plotly to be rendered inline in the notebook.
    plotly.offline.init_notebook_mode()
    # Configure the trace.
    decoy = go.Scatter3d(
        x=pos['x'][7:], y=pos['y'][7:], z=pos['z'][7:],
        mode = 'markers',
        marker={
            'size': 3,
            'opacity': 1,
            'color' : "lightblue",
        },
        text=[str(i+7) for i in range(len(pos['x'][7:])) ])
    truepos = go.Scatter3d(
        x=pos['x'][:7], y=pos['y'][:7], z=pos['z'][:7],
        mode = 'markers',
        marker={
            'size': 5,
            'opacity': 0.8,
            'color' : "blue",
        },
        text=[str(i) for i in range(7) ])
    layout = go.Layout(title='Docking decoys and True position of '+complex +"'s ligand'",
    hovermode= 'closest',)
    margin={'l': 0, 'r': 0, 'b': 0, 't': 10}
    data = [decoy,truepos]
    plot_figure = go.Figure(data=data, layout=layout)
    # Render the plot.
    plotly.offline.iplot(plot_figure,filename='hover-chart-basic')

def colorsFromRmsd(rmsds):
    colors=[]
    for u in rmsds:
        if u<5:
            colors.extend(['green'])
            continue
        elif 5<u<10 :
            colors.extend(['red'])
            continue
        else :
            colors.extend(['blue'])
    return colors


def countNative(rmsds, cutoff=5):
    counts={5:0, 10:0, 20:0, 100:0, 200:0, "out":0}
    x=0
    for i in rmsds:
        if i<=cutoff:
            if x<200:
                counts[200]+=1
                if x<100:
                    counts[100]+=1
                    if x<20:
                        counts[20]+=1
                        if x<10:
                            counts[10]+=1
                            if x<5 :
                                counts[5]+=1
            else:
                counts["out"]+=1
        x+=1
    return counts

def eval_natives(natives,n):
    """Takes a dictionnary {complexe : countNative(rmsd)} and a top limit for good complexes (prediction in top n)
    and returns a list of succeded complexes and failed ones """
    positives=0
    good=[]
    bad=[]
    for c in natives:
    #     print(natives[c])
        if natives[c][n]>0:
            positives+=1
            good.append(c)
        if natives[c][200]==0 and natives[c]["out"]==0:
            bad.append(c)
    return (good, bad)



def multiplePlots(num, size=(20,10), ylim=None, xlim=None, title=None):
    fig=plt.figure(figsize=size)
    if title:
        fig.suptitle(title)
    axlist=[]
    for i in range(1,num+1):
        ax=plt.subplot(1, num, i)
        if ylim:
            ax.set_ylim(*ylim)
        if xlim:
            ax.set_xlim(*xlim)
        axlist.append(ax)
    return axlist

def scatter3d(data,groups, title=None):
    fig = plt.figure(figsize=(20,20))
    ax = fig.gca(projection='3d')

    sc=ax.scatter(data['x'], data['y'], data['z'],'o',c=groups, s=[i*10 for i in groups], cmap='autumn', picker=True)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    # Rotate it
    ax.view_init(45, 45)

    ### NEED TO FIND A WAY TO HANDLE HOVERING OVER POINST
    ### NEED TO ADD COLOR BAR
    ### NEED TO TAKE IN POSES


    plt.show()
