#Here we will read a model and produce an escher map so that all reactions
#are placed on the canvas.

class REACTION:
    def __init__(self,name,bigg_id,reversibility,gene_R,GENES,metabolites):
        self.name=name
        self.bigg=bigg_id
        self.rev=reversibility
        self.gene_rule=gene_R
        self.genes=GENES
        self.metabolites=metabolites


class MAP_REACTION:
    def __init__(self,name,bigg_id,reversibility,gene_R,genes,center):


class HUH:
    def __init__(self,center):
        self.c=center
    @property
    def c(self):
        return self._c
    @c.setter
    def c(self,c):
        self._c=c
        self.x=self._c+50

    


                 

#Note, if you do vars(REACTION) it will tell you all the .name and stuff on that object

import json

NAME='iJN678'
FileIN=NAME+'.json'
OutName=NAME+'_Escher_Test_Map'
FileOUT=OutName+'.json'
print('loading COBRA model '+FileIN)
with open(FileIN,'r') as f:
    M=json.load(f)
React=M.get('reactions')

GE=M.get('genes')
GEID=list()
for i in GE:
    GEID.append(i['id'])
GENES=dict(zip(GEID,GE))

Reactions=list()
for i in React:
    if i['lower_bound']==0:
        rev=False
    else:
        rev=True
    G=list()
    for j in i['gene_reaction_rule'].replace(' or ','').replace(' and ','').replace('(','').lstrip().split(')'):
        if j!='':
            G.append(GENES.get(j))
    Reactions.append(REACTION(i['name'],i['id'],rev,i['gene_reaction_rule'],G,i['metabolites']))
print('Have '+str(len(Reactions))+' reactions')

#Now build header, Header can be anything
#at the end we will assemble a string:
# [ {HEADER},{REST OF Stuff}]
HEADER='"map_name":"'+OutName+'","map_id":"OffLineMark1","map_description":"TestMap_Made_in_Python_to_ape_Escher_Map","homepage":"NONE","schema":"https://escher.github.io/escher/jsonschema/1-0-0#"'


#Now build Escher map, reactions and nodes with connections with gross position to be decided

#It appears that the logic for choosing primary metabolites is as follows:
#In the cobra model, the fist listed positive and negative valued metabolite are the primary metabolite
#We can't do the exact same thing here because the collection of metabolites is unordered. But this and
#the documentation for Escher seems to suggest that Primary/Secondary is arbitrary and to be decided by
#human intervention.  Program just picks two so the auto system has something to go on.

#Note, Primary nodes are arbitrary. They are to be chosen by the user at map draw time. when a reaction is
#selected they use the C key to rotate the metabolites, can toggle primary/secondary status with the p key.


#Maybe a class of reaction on map, contains center, relative position and name of metabolite nodes, and the arrows associated
#Think to the future, how to find and condense all the nodes.


#Now place the reactions on the canvas, define the canvas
