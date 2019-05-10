#! /usr/bin/env python3

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from mpl_toolkits.mplot3d import Axes3D
from core_clustering import rankCluster
import pickle
from dockingPP import parse

import plotly.plotly as py
# Import dependencies
import plotly
import plotly.graph_objs as go

# #labels = range(1, 11)
# plt.figure(figsize=(10, 7))
# plt.subplots_adjust(bottom=0.1)
# plt.scatter('x','y', data=madata, label='True Position')
# plt.plot(madata['x'][:6],madata['y'][:6],'rs')

#3D viewing
# Make the plot

# zD = pickle.load(open('/Users/jprieto/Docking/scripts/1BJ1_4000_full.pickle','rb'))
# zD=parse('/Users/jprieto/Docking/data/subsamples/1BJ1_r-1BJ1_l.detail')
# madata=zD.dictPos
# groups=rankCluster(zD.pList,5)
# Données a definir
# S=[(sqrt(10000-i)/10) for i in rank]
# D=[(100-rmsd)/20 for rmsd in rmsds]

class Scores(object):
    def __init__(self, file):
        self.data=None
        with open(file,'r') as f:
            self.data=[line.split("\t") for line in f.readlines()[3:]]
        self.columns={"pose":0,"s_size":1,"res_fr_sum":2,"res_mean_fr": 3, \
        "res_log_sum": 4, "res_sq_sum": 5, "c_size": 6, "con_fr_sum": 7,"con_mean_fr" : 8,"con_log_sum" : 9,"con_sq_sum" : 10}
        self.poses=None

    def setPoses(self, pList):
        self.poses=pList

    @property
    def scoresDict(self):
        scoresdict={}
        for info in data:
            scoresdict[int(info[0])]=info[1:]
        return scoresDict

    @property
    def scoresList(self):
        return self.data

    @property
    def rmsds(self):
        if not self.poses :
            raise Exception("Please define poses using setPoses")
        return [poses.rmsd for pose in self.poses]


    def rankedRmsds(self, ranked):
        rmsds=[]
        for i in ranked:
            rmsds.append(self.poses[i].rmsd)
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

    def rmsdGraphGenerator(self, ranks, start=0, stop=None, plot=None ):
        pl=plot if plot else plt
        x=0
        if not stop:
            stop=len(ranks)-1
        rmsds=self.rankedRmsds(ranks)
        colors=colorsFromRmsd(rmsds)[int(start):int(stop)]
        x=[i for i in range(len(colors))]
        y=rmsds[int(start):int(stop)]
        # print(len(x))
        # print(len(y))
        pl.scatter(x, rmsds[int(start):int(stop)], c=colors)
        # for i in rmsds[int(start):int(stop)]:
        #     if self.poses[i].rmsd<5:
        #         pl.plot(x,self.poses[i].rmsd, 'go')
        #         pl.annotate('%s' % self.poses[i].id, xy=(x,self.poses[i].rmsd), textcoords='data')
        #
        #     elif 5<self.poses[i].rmsd<10 :
        #         pl.plot(x,self.poses[i].rmsd, 'ro')
        #     # elif  self.poses[i].rmsd<20 :
        #     else :
        #         pl.plot(x,self.poses[i].rmsd, 'bo')
        #     x+=1


    def rankPoses(self, element="res_mean_fr", start=0, stop=None):
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

        if not self.poses :
            raise Exception("Please define poses using setPoses() funciton")
        if not stop:
            stop=len(self.data)-1
        col=self.columns[element]
        r=sorted(self.data[int(start):int(stop)],key=lambda o:float(o[col]),reverse=True)
        sorted_i=[int(pose[0]) for pose in r]
        return sorted_i

    def histRmsd(self,plot=None, stop=None):
        pl= plot if plot else plt
        s=stop if stop else len(self.data)-1
        pl.hist(sorted([i[1] for i in self.data[:stop]]))


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

def countNative(rmsds):
    firstC=0
    firstK=0
    firstKK=0
    out=0
    x=0
    for i in rmsds:

        if i<5:
            if x<100:
                firstC+=1
            if x<1000:
                firstK+=1
            if x<2000:
                firstKK+=1
            if x>2000:
                out+=1
    return [(100,firstC),(1000,firstK),(2000,firstKK),("out",out)]

def countValid(pList):
    for cluster in clusters:
        for pose in cluster:
            print(pose)


def multiplePlots(num, size=(20,10)):
    fig=plt.figure(figsize=size)
    axlist=[]
    for i in range(1,num+1):
        ax=plt.subplot(1, num, i)
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

def plotly3D(pos, S, D):
    # Configure Plotly to be rendered inline in the notebook.
    plotly.offline.init_notebook_mode()
    # Configure the trace.
    trace = go.Scatter3d(
        x=pos['x'],  # <-- Put your data instead
        y=pos['y'],  # <-- Put your data instead
        z=pos['z'],  # <-- Put your data instead
        mode = 'markers',
        marker={
            'size': S,
            'opacity': 1,
            'color' : D,
            'colorscale' : 'Jet'
        }
    )

    # Configure the layout.
    layout = go.Layout(
        margin={'l': 0, 'r': 0, 'b': 0, 't': 0}
    )

    data = [trace]

    plot_figure = go.Figure(data=data, layout=layout)

    # Render the plot.
    plotly.offline.iplot(plot_figure)
