#! /usr/bin/env python3

import matplotlib.pyplot as plt
import math

from src.core_clustering import rankCluster
import pickle
from math import sqrt
# Import dependencies
import plotly
import plotly.graph_objs as go


class Scores(object):
    """This class allows to read a rescoring file"""
    def __init__(self, filename=None, data=None):
        self.data=None
        if filename :
            self.loadFile(filename)
        elif data:
            self.loadScores(data)
        self.columns={"pose":0,"s_size":1,"res_fr_sum":2,"res_mean_fr": 3, "res_log_sum": 4, "res_sq_sum": 5, "c_size": 6, "con_fr_sum": 7,"con_mean_fr" : 8,"con_log_sum" : 9,"con_sq_sum" : 10, "rmsd": 11}
        self.poses=None

    @property
    def originalRank(self):
        if not self.poses:
            raise Exception("Please define poses using setPoses")
        return [i for i in range(len(self.poses))]

    def loadFile(self,file):
        with open(file,'r') as f:
            self.data=[line.strip("\n").split("\t") for line in f.readlines()[3:]]
        return True

    def loadScores(self,scores):
        self.data=scores

    def setPoses(self, pList):
        self.poses=pList

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
        if self.poses :
            return [pose.rmsd for pose in self.poses]
        else :
            try:
                return [float(info[self.columns["rmsd"]].strip('\n')) for info in self.data]
            except IndexError:
                raise Exception("Please define poses using setPoses")


    def getScore(self,poseid, score):
        return self.scoresDict[int(poseid)-1][self.columns[score]]


    def rankedRmsds(self, rankedPoses):
        rmsds=[]
        for i in rankedPoses:
            rmsds.append(self.rmsds[i])
        return rmsds

    def coordDict(self, start=0, stop=None):
        """dockingPP poses position dictorizer """
        coord={'x':[],'y':[],'z':[],'a1':[],'a2':[],'a3':[]}
        if not stop:
            stop=len(self.poses)-1
        for pose in self.poses[int(start):int(stop)]:
            coord['x'].append(pose.translate[0])
            coord['y'] .append(pose.translate[1])
            coord['z'] .append(pose.translate[2])
            coord['a1'] .append(pose.euler[0])
            coord['a2'] .append(pose.euler[1])
            coord['a3'] .append(pose.euler[2])
        return coord

    def rmsdGraphGenerator(self, rankedPoses, start=0, stop=None, plot=None, title=None ):
        pl=plot if plot else plt
        if title:
            pl.set_title(title)
        x=0
        if not stop:
            stop=len(rankedPoses)-1
        rmsds=self.rankedRmsds(rankedPoses)
        colors=colorsFromRmsd(rmsds)[int(start):int(stop)]
        x=[i if i!=0 else 0 for i in range(len(colors))]
        y=rmsds[int(start):int(stop)]
        # print(len(x))
        # print(len(y))
        pl.scatter(x, rmsds[int(start):int(stop)], c=colors)



    def ranks(self, element="res_mean_fr", start=0, stop=None):
        """Returns ranks for poses in the original order
        At position 0 we have pose 0 with rank 37 : rank[0] = 37
        Means the pose in original rank 0 is in 37th position according to rescoring
        element:
        "s_size","res_fr_sum","res_mean_fr","res_log_sum","res_sq_sum",
        "c_size","con_fr_sum","con_mean_fr","con_log_sum","con_sq_sum"  """

        sorted_i=self.rankedPoses(element=element, start=start, stop=stop)
        rank=[u[0] for u in sorted([(u,v) for u,v in enumerate(sorted_i)], key=lambda o:o[1])]
        return rank

    def rankedPoses(self, element="res_mean_fr", start=0, stop=None):
        """Returns poses' ids,  sorted according to rescoring element : ]
        "s_size","res_fr_sum","res_mean_fr","res_log_sum","res_sq_sum",
        "c_size","con_fr_sum","con_mean_fr","con_log_sum","con_sq_sum"
        Means the pose classified 1st is 1st position of the list.   """

        try :
            p=self.poses
        except :
            raise Exception("Please define poses using setPoses() function")
        if not stop:
            stop=len(self.data)
        col=self.columns[element]
        r=sorted(self.data[int(start):int(stop)],key=lambda o:float(o[col]),reverse=True)
        sorted_i=[int(pose[0]) for pose in r]
        return sorted_i

    def histRmsd(self,plot=None, stop=None):
        pl= plot if plot else plt
        s=stop if stop else len(self.data)
        pl.hist(sorted([i[1] for i in self.data[:stop]]))

    def trace(self, rank, name):
        colorscale=make_colorScale()

        rmsds=self.rmsds
        pos=self.coordDict(start=0, stop=len(rank)-1)
        C=[r+1 for r in rank]
        # print(C)
        S=[ 15 if rmsd > 5 else 30 for rmsd in rmsds ]
        # S=[ 15]
        # Configure the trace.
        trace = go.Scatter3d(
            x=pos['x'], y=pos['y'], z=pos['z'],
            mode = 'markers',
            marker={
                'size': S,
                'opacity': 1,
                'color' : C,
                'colorscale' : colorscale,
                'colorbar': dict(title = 'Ranks',
            titleside = 'top',
            tickmode = 'array',
            x=-0.25)
            },
            name=name,
            text=['rank=' + str(rank[i]+1 ) + ', rmsd='+ str(self.rmsds[i]) for i in range(len(rank))])

        return trace

    def plot3D(self, ranks, name=' ', title='Docking decoys'):
        """ Warning : ranks must be in the original order, since positions and rmsds are in the original order """
        # Configure Plotly to be rendered inline in the notebook.
        plotly.offline.init_notebook_mode()
        trace=self.trace(ranks, name)
        layout = go.Layout(autosize=False,
        width=650,
        height=400,
        title=title,
        hovermode= 'x',
        hoverlabel=dict(bgcolor='white', bordercolor='black'),
        margin={'l': 0, 'r': 0, 'b': 0, 't': 10})
        data = [trace]
        plot_figure = go.Figure(data=data, layout=layout)
        # Render the plot.
        plotly.offline.iplot(plot_figure,filename='Docking Decoys')
        return plot_figure

    # NB : This function should take different sets of points and not sc object. --> put to core_visuals ?

def multiPlot3D(sc, ranks, names, title='Docking decoys'):
    """ Warning : rank must be in the original order, since positions and rmsds are in the original order
    ranks is a list of lists of ranks"""
    # Configure Plotly to be rendered inline in the notebook.
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
    width=600,
    height=400,
    title=title, hovermode= 'closest',)
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



def multiplePlots(num, size=(20,10), ylim=None):
    fig=plt.figure(figsize=size)
    axlist=[]
    for i in range(1,num+1):
        ax=plt.subplot(1, num, i)
        if ylim:
            ax.set_ylim(*ylim)
        axlist.append(ax)
    return axlist

def scatter3d(data,groups):
    fig = plt.figure(figsize=(20,20))
    ax = fig.gca(projection='3d')
    ax.scatter(data['x'][7:], data['y'][7:], data['z'][7:],'o',c=groups, cmap='rainbow')
    # Rotate it
    ax.view_init(45, 45)
    plt.set_cmap('rainbow')
    plt.show()
