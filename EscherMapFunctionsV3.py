SegID = 1
def next_id():
    global SegID
    res = SegID
    SegID += 1
    return str(res)



class Nodes:
    #This will be a dictionary with ID incrementer
    #On initalization, you feed it the starting ID
    #We will give it the Add property which will add a node and increment the ID
    def __init__(self,StartingID=0):
        self.C = str(StartingID)
        self.N=dict()
    def Add(self,node):
        self.N.update({str(self.C):node})
        self.C=str(int(self.C)+1)
class Node:
    def __init__(self,NType,x,y,bigid=None,name=None,labx=None,laby=None,prim=None):
        self.N=dict()
        self.node_type=NType
        self.x=x
        self.y=y
        if bigid != None:
            self.bigid=bigid
            self.name=name
            self.labx=labx
            self.laby=laby
            self.prim=prim
    def Print(self): #must make this call before getting self.N
        self.N.update({'node_type':self.node_type})
        self.N.update({'x':self.x})
        self.N.update({'y':self.y})
        try:
            self.N.update({'bigg_id':self.bigid})
            self.N.update({'name':self.name})
            self.N.update({'label_x':self.labx})
            self.N.update({'label_y':self.laby})
            self.N.update({'node_is_primary':self.prim})
        finally:
            return(self.N)
class Reaction:
    def __init__(self,name,biggid,reversibility,labx,laby,genereactionrule,genes,metabolites,segments):
        self.R=dict()
        self.name=name
        self.labx=labx
        self.laby=laby
        self.biggid=biggid
        self.rev=reversibility
        self.genrule=genereactionrule
        self.gen=genes
        self.met=metabolites
        self.seg=segments
    def Print(self): #must make this call before getting self.N
        self.R.update({'name':self.name})
        self.R.update({'bigg_id':self.biggid})
        self.R.update({'reversibility':self.rev})
        self.R.update({'label_x':self.labx})
        self.R.update({'label_y':self.laby})
        self.R.update({'gene_reaction_rule':self.genrule})
        self.R.update({'genes':self.gen})
        self.R.update({'metabolites':self.met})
        self.R.update({'segments':self.seg})
        return(self.R)

    
class Reactions:
    #This will be a dictionary with ID incrementer
    #On initalization, you feed it the starting ID
    #We will give it the Add property which will add a node and increment the ID
    def __init__(self,StartingID=0):
        self.C = StartingID
        self.R=dict()
    def Add(self,reaction):
        self.R.update({str(self.C):reaction})
        self.C=str(int(self.C)+1)


def MAKE_ESCHER_MAP_FILE_FROM_CLASSES(Nodes, Reactions, Name):
    #Makes an escher map file called Name.json where Name is variable Name.
    #Nodes is a object of class Nodes, Reactions is of class Reactions.
    import json
    HEADER={'map_name':Name,'map_id':'MarkProgEscher2MakeEscherMapFromCLasses','map_description':'Map was build using a python  script written by Mark Layer on 8/4/2015','homepage':"https://escher.github.io",'schema':'https://escher.github.io/escher/jsonschema/1-0-0#'}

    MaxX=0
    MaxY=0
    MinX=0
    MinY=0
    for n in Nodes.N:
        MaxX=max(Nodes.N[n]['x'],MaxX)
        MaxY=max(Nodes.N[n]['y'],MaxY)
        MinX=min(Nodes.N[n]['x'],MinX)
        MinY=min(Nodes.N[n]['y'],MinY)
    Canvas={'x':MinX-50,'y':MinY-50,'width':(MaxX-MinX)+100,'height':(MaxY-MinY)+100}
    
    REST={'reactions':Reactions.R,'text_labels':{},'canvas':Canvas,'nodes':Nodes.N}
    PrintString=str('['+json.dumps(HEADER)+','+json.dumps(REST)+']')
    PrintString=PrintString.replace('"false"','false').replace('"true"','true').replace('"null"','null')
    open(Name+'.json','w').write(PrintString)


def PLOT_REACTION_RELATIVE_TO_MIDMARKER(Nodes,Reaction,MetaboliteIdNameDict):
    #Feed this a reaction and the global node class, it will generate nodes, place them relative to midmarker which is centered at (0,0)
    #The reaction is at default oriented along the X axis, with poroduct marker at (100,0) and reactant marker at (-100,0) 
    #at zero zero, and make segments which point to correct stuff
    #idea is later functions will rotate/translate a plotted reaction
    #node coordinates will change by following segment node IDs
    import math
    
    MidID=Nodes.C
    Nodes.Add(Node('midmarker',0,0).Print())
    ProdID=Nodes.C
    Nodes.Add(Node('multimarker',100,0).Print())
    ReactID=Nodes.C
    
    Nodes.Add(Node('multimarker',-100,0).Print())
    NEWNODES=[MidID,ProdID,ReactID]
    Prods=list()
    Reacts=list()#contains Node ID of Prods/Reacts for segment construct
    ProdCount=0
    ReactCount=0
    for Met in Reaction.met: #Reactions.met is a list of dictionary entries
        Coeff=Met['coefficient']
        if Coeff>0:
            PMY=((ProdCount%2)-.5)*(2)
            OffSet=(ProdCount//2)*50
            Prods.append(Nodes.C)
            NEWNODES.append(Nodes.C)
            Nodes.Add(Node('metabolite',100+OffSet,PMY*(100+OffSet),Met['bigg_id'],MetaboliteIdNameDict[Met['bigg_id']],100+OffSet-20,PMY*(100+OffSet+20),'false').Print())
            ProdCount+=1
        else:
            PMY=((ReactCount%2)-.5)*(2)
            OffSet=(ReactCount//2)*50
            Reacts.append(Nodes.C)
            NEWNODES.append(Nodes.C)
            Nodes.Add(Node('metabolite',(-100)-OffSet,PMY*(100+OffSet),Met['bigg_id'],MetaboliteIdNameDict[Met['bigg_id']],(-100)-OffSet-20,PMY*(100+OffSet+20),'false').Print())
            ReactCount+=1
    segments = dict()
    segments.update({next_id():{'from_node_id':MidID,'to_node_id':ProdID,'b1':'null','b2':'null'}})
    segments.update({next_id():{'from_node_id':MidID,'to_node_id':ReactID,'b1':'null','b2':'null'}})
    for NodeID in Prods:
        segments.update({next_id():{'from_node_id':ProdID,'to_node_id':NodeID,'b1':{'x':140,'y':0},'b2':{'x':Nodes.N[NodeID]['x'],'y':Nodes.N[NodeID]['y']-math.copysign(40,Nodes.N[NodeID]['y'])}}})
    for NodeID in Reacts:
        segments.update({next_id():{'from_node_id':ReactID,'to_node_id':NodeID,'b1':{'x':-140,'y':0},'b2':{'x':Nodes.N[NodeID]['x'],'y':Nodes.N[NodeID]['y']-math.copysign(40,Nodes.N[NodeID]['y'])}}})      
    Reaction.seg=segments
    return NEWNODES

    
def GROW_FROM_EXTERNAL_LIKE_CRYSTALS(Cobra_Model):
    import json
    with open(Cobra_Model,'r') as f:
        M=json.load(f)
    MetaboliteID_To_Name=dict()
    AllNodes=Nodes(0)
    AllReactions=Reactions(0)
    for met in M['metabolites']:
        MetaboliteID_To_Name.update({met['id']:met['name']})
    #print(list(M['metabolites'])[0])
    #print(list(MetaboliteID_To_Name)[0]+' '+MetaboliteID_To_Name[list(MetaboliteID_To_Name)[0]])
    GeneID_To_Name=dict()
    for gene in M['genes']:
        GeneID_To_Name.update({gene['id']:gene['name']})
    for React in M['reactions']:
        if (React['lower_bound']<0):
            Rev='true'
        else:
            Rev='false'
        GeneList=list()
        for UpG in React['gene_reaction_rule'].split(' or '):
            for DG in UpG.split(' and '):
                gene=DG.replace('(','').replace(')','').replace(' ','')
                if gene!='':
                    GeneList.append({gene:GeneID_To_Name[gene]})
        MetList=list()
        for Met in React['metabolites']:
            MetList.append({'bigg_id':Met,'coefficient':React['metabolites'][Met]})
        AllReactions.Add(Reaction(React['name'],React['id'],Rev,0,50,React['gene_reaction_rule'],GeneList,MetList,'').Print())
    NodesCollections=dict()
    CountNodes=0
    for RR in AllReactions.R:
        ReactionsNode=Nodes(CountNodes)
        PLOT_REACTION_RELATIVE_TO_MIDMARKER(ReactionsNode,AllReactions.R[RR],MetaboliteID_To_Name)
        NodesCollections.update({RR:ReactionsNode})
        if int(ReactionsNode.C)>20:
            break

    H=Reactions(0)
    H.Add(AllReactions.R[RR])
    MAKE_ESCHER_MAP_FILE_FROM_CLASSES(NodesCollections[RR],H,'Tester111')
    #print(str(AllNodes.N).replace("'",'"').replace('"false"','false').replace('"true"','true'))
    #print('*********************')
    #print('"'+RR+'":'+str(AllReactions.R[RR]).replace("'",'"').replace('"null"','null').replace('"false"','false').replace('"true"','true'))



def GROWFROMSINGLEREACTION(Cobra_Model,ReactionID,ComplexCutoff=3,FirstNumber=3,ReactionsToIgnore=[],OutName='',KeepGoing=None,OriginModel='',ManyMaps=''):
    #Give it a reaction ID and a cobra model.

    #Ultimate plan: It will plot that reaction in center, plot all reactions coming off of it
    #grow these into a single network
    #has some cutoff beyond which not to consider a connection
    #link to that reaction
    import json
    from collections import Counter
    MetCount=Counter()
    PrefixCount=Counter()
    with open(Cobra_Model,'r') as f:
        M=json.load(f)

    if OriginModel!='':
        with open(OriginModel,'r') as f:
            OLDModel=json.load(f)
    else:
        OLDModel=M

    MetaboliteID_To_Name=dict()
    AllNodes=Nodes(1000)
    AllReactions=Reactions(0)
    CoreReaction=dict()
    for met in M['metabolites']:
        MetaboliteID_To_Name.update({met['id']:met['name']})
    GeneID_To_Name=dict()
    for gene in M['genes']:
        GeneID_To_Name.update({gene['id']:gene['name']})
    for React in M['reactions']:
        if (React['lower_bound']<0):
            Rev='true'
        else:
            Rev='false'
        GeneList=list()
        for UpG in React['gene_reaction_rule'].split(' or '):
            for DG in UpG.split(' and '):
                gene=DG.replace('(','').replace(')','').replace(' ','')
                if gene!='':
                    GG=dict()
                    GG.update({"name":GeneID_To_Name[gene]})
                    GG.update({"bigg_id":gene})
  
                    GeneList.append(GG)
        MetList=list()
        for Met in React['metabolites']:
            MetList.append({'bigg_id':Met,'coefficient':React['metabolites'][Met]})
            
        if React['id'] not in ReactionsToIgnore:
            AllReactions.Add(Reaction(React['name'],React['id'],Rev,0,50,React['gene_reaction_rule'],GeneList,MetList,''))
            if React['id']==ReactionID:
                CoreReaction=[int(AllReactions.C)-1,AllReactions.R.pop(str(int(AllReactions.C)-1))]

    for React in OLDModel['reactions']:
        for Met in React['metabolites']:
            MetCount.update([Met])
            PrefixCount.update([Met[:-2]])
############################           
    #HARDCODE=['o2_x','nadp_x','nadp_m','nadp_h','nadp_c','nadph_x','nadph_m','nadph_h','nadph_c','atp_c','atp_h','atp_x','atp_m','o2_c','o2_m','o2_x','dtdp_c','dttp_c','ppi_c','ppi_m','ppi_h','ppi_x','amp_c','amp_h','amp_m','amp_x']    
    NeverCombine=set ()# (HARDCODE)
    #PREFIXLIST=set(['o2','nadp','nadph','atp','adp','dtdp','amp','h','ppi','dttp'])

    CountCutoff=5
    for i in PrefixCount.items():
        if i[1]>min(max(ComplexCutoff,FirstNumber)+5,CountCutoff):
            PREFIXLIST.add(i[0])
        #experiment, check prefix minus the final two characters, Maybe prefix list to catch all compartments
##########################################
            
    for i in MetCount.items():
        if i[1]>min(max(ComplexCutoff,FirstNumber)+5,CountCutoff):
            NeverCombine.add(i[0])
        #experiment, check prefix minus the final two characters, Maybe prefix list to catch all compartments
        if i[0][:-2] in PREFIXLIST:
            NeverCombine.add(i[0])


    ContinueFlag=1
    GlobolX=0
    GloH=0
    TotR=list()
    TotN=list()
    CCCCC=0
    AllReacts=Reactions(100000)
    while len(AllReactions.R)>0 and ContinueFlag==1:


        if KeepGoing==None or KeepGoing==0:
            ContinueFlag=0
        

        if len(CoreReaction)==0:
            print('new '+str(len(AllReactions.R)))
            CoreReaction=[list(AllReactions.R)[0],AllReactions.R.pop(list(AllReactions.R)[0])]
        NumChange=1
        KeptReacts=[CoreReaction]
        
        Links=dict()

        MtoReact=dict()
        for m in M['metabolites']:
            if m['id'] not in NeverCombine:
                ReactList=[]
                for RR in AllReactions.R:
                    PP=AllReactions.R[RR]
                    for om in PP.met:
                        if m['id']==om['bigg_id']:
                            ReactList.append(RR)
                            break
                MtoReact.update({m['id']:ReactList})
        while(NumChange>0):
            print(NumChange)
            NumChange=0
            #MtoReact=dict()
            for m in MtoReact:
                for R in MtoReact[m]:
                    if R not in AllReactions.R:
                        MtoReact[m].pop(MtoReact[m].index(R))
                        
           # for m in M['metabolites']:
            #    ReactList=[]
             #   for RR in AllReactions.R:
              #      PP=AllReactions.R[RR]
               #     for om in PP.met:
                #        if m['id']==om['bigg_id']:
                 #           ReactList.append(RR)
                  #          break
               # MtoReact.update({m['id']:ReactList})


             #oneway links from Key to List of children
            
            ProcessingReactions=list(KeptReacts)
            KeptReacts=list()
            
            while len(ProcessingReactions)>0:
                CCC=ProcessingReactions.pop(0)
                CurrentR=CCC[1]
                LinkList=[]
                for L in CurrentR.met:
                    ID=L['bigg_id']
                    if ID not in NeverCombine:
                        AdjList=MtoReact[ID]
                        if CurrentR==CoreReaction[1]:
                            ComCut=FirstNumber
                        else:
                            ComCut=ComplexCutoff
                        if not(len(AdjList)>ComCut):
                            for RItem in AdjList:
                                if RItem in AllReactions.R:
                                    ProcessingReactions.append([RItem,AllReactions.R.pop(RItem)])
                                    NumChange=NumChange+1
                                    LinkList.append(RItem)
                KeptReacts.append(CCC)
                OlList=list()
                if(CCC[0] in Links):
                    OlList=Links[CCC[0]]
                for k in LinkList:
                    OlList.append(k)
                Links.update({CCC[0]:OlList})
        Root=KeptReacts.pop(0)
        ConsideredReactions=dict()
        for r in KeptReacts:
            ConsideredReactions.update({r[0]:r[1]})
        BigList=TREEPLOT(ConsideredReactions,Links,AllNodes,MetaboliteID_To_Name,Root,NeverCombine)
        TranslateNodes(BigList[1],AllNodes,GlobolX+BigList[0][0]/2,0)
        TranslateReactions(BigList[2],GlobolX+BigList[0][0]/2,0)
        GlobolX=GlobolX+BigList[0][0]+150
        GloH=max(BigList[0][1],GloH)
        CoreReaction=dict()
        if ManyMaps!='':
            COMBINE_NODES(AllNodes,BigList[2],NeverCombine)
            AllReacts=Reactions(100000)
            RESET_ALL_BEVELS(AllNodes,BigList[2])
            for R in BigList[2]:
                AllReacts.Add(R.Print())
            MAKE_ESCHER_MAP_FILE_FROM_CLASSES(AllNodes,AllReacts,OutName+str(CCCCC))
            AllNodes=Nodes(1000)
            CCCCC=CCCCC+1
            GlobolX=0
            GloH=0
        else:
            for R in TotR:
                AllReacts.Add(R.Print())
            for R in BigList[2]:
                TotR.append(R)
    if ManyMaps=='':
        COMBINE_NODES(AllNodes,TotR,NeverCombine)
        AllReacts=Reactions(100000)
        RESET_ALL_BEVELS(AllNodes,TotR)
        for R in TotR:
            AllReacts.Add(R.Print())
        
        MAKE_ESCHER_MAP_FILE_FROM_CLASSES(AllNodes,AllReacts,OutName)

def RESET_ALL_BEVELS(AllNodes,AllReacts):
    for R in AllReacts:
        
        for S in R.seg:
            SEG=R.seg[S]
            #Remember, b1 is based at the from node, b2 is based at the to node
            if SEG['b1']!='null' and SEG['b1']!=None:
                ID=R.seg[S]['to_node_id']
                OID=R.seg[S]['from_node_id']
                TONODE=AllNodes.N[ID]
                FROMNODE=AllNodes.N[OID]
                FromX=FROMNODE['x']
                FromY=FROMNODE['y']
                ToX=TONODE['x']
                ToY=TONODE['y']
                BTWO=GET_FORTY_AWAY(FromX,FromY,ToX,ToY)
                SEG['b2']['x']=BTWO[0]
                SEG['b2']['y']=BTWO[1]
                BONE=GET_FORTY_AWAY(ToX,ToY,FromX,FromY)
                SEG['b1']['x']=BONE[0]
                SEG['b1']['y']=BONE[1]
                


def COMBINE_NODES(AllNodes,ReactionList,NeverCombine):
    KeepNodes=dict()
    Remap=dict()
    for N in AllNodes.N:
        if AllNodes.N[N]['node_type']=='metabolite':
            if AllNodes.N[N]['bigg_id'] not in NeverCombine:
                if AllNodes.N[N]['bigg_id'] not in KeepNodes: #Only keep the node with the lowest Y coordinate
                    KeepNodes.update({AllNodes.N[N]['bigg_id']:N})
                else:
                    OldY=AllNodes.N[KeepNodes[AllNodes.N[N]['bigg_id']]]['y']
                    NewY=AllNodes.N[N]['y']
                    if NewY<OldY:
                        KeepNodes.update({AllNodes.N[N]['bigg_id']:N})

    for N in AllNodes.N:
        if AllNodes.N[N]['node_type']=='metabolite':
            if AllNodes.N[N]['bigg_id'] not in NeverCombine:
                if N != KeepNodes[AllNodes.N[N]['bigg_id']]:
                    Remap.update({N:KeepNodes[AllNodes.N[N]['bigg_id']]})
    for N in Remap:
        del AllNodes.N[N]
        AllNodes.N[Remap[N]]['node_is_primary']='true'
    for R in ReactionList:
        for S in R.seg:
            if R.seg[S]['to_node_id'] in  Remap:
                ID=R.seg[S]['to_node_id']
                OID=R.seg[S]['from_node_id']
#Remember, b1 is based at the from node, b2 is based at the to node
                #here we change the to node, so b2 needs to change
                #Got from node, got new to node, can place b2 40 units from to node pointing at from node
                FromX=AllNodes.N[OID]['x']
                FromY=AllNodes.N[OID]['y']
                ToX=AllNodes.N[Remap[ID]]['x']
                ToY=AllNodes.N[Remap[ID]]['y']
                POINT=GET_FORTY_AWAY(FromX,FromY,ToX,ToY)
                R.seg[S]['to_node_id']=Remap[ID]
                R.seg[S]['b2']['x']=POINT[0]
                R.seg[S]['b2']['y']=POINT[1]              


            if R.seg[S]['from_node_id'] in  Remap:
                OID=R.seg[S]['to_node_id']
                ID=R.seg[S]['from_node_id']
                R.seg[S]['from_node_id']=Remap[ID]
                FromX=AllNodes.N[ID]['x']
                FromY=AllNodes.N[ID]['y']
                ToX=AllNodes.N[Remap[OID]]['x']
                ToY=AllNodes.N[Remap[OID]]['y']
                POINT=GET_FORTY_AWAY(ToX,ToY,FromX,FromY)
                R.seg[S]['b1']['x']=POINT[0]
                R.seg[S]['b1']['y']=POINT[1] 
def GET_FORTY_AWAY(FromX,FromY,ToX,ToY):
    import math
    VABx=FromX-ToX
    VABy=FromY-ToY
    VABL=math.sqrt(VABx*VABx + VABy*VABy)
    if VABL!=0:
        Ux=VABx/VABL
        Uy=VABy/VABL
    else:
        Ux=0
        Uy=1
    X=Ux*40
    Y=Uy*40
    Px=ToX+X
    Py=ToY+Y
    return([Px,Py])
def ROTATE_POINT(Point,angle,center):
    import math
    s=math.sin(angle)
    c=math.cos(angle)
    PX=Point[0]-center[0]
    PY=Point[1]-center[1]
    return((PX*c-PY*s,PX*s+PY*c))

def ROTATE_BY_DEG_NODE(NODELIST,AllNodes,Reaction,ROT):
    import math
    angle=math.radians(ROT)
    NODE=-1
    maxx=0
    minx=1000
    maxy=0
    miny=1000
    for N in NODELIST:
        maxx=max(maxx,AllNodes.N[N]['x'])
        maxy=max(maxy,AllNodes.N[N]['y'])
        minx=min(minx,AllNodes.N[N]['x'])
        miny=min(miny,AllNodes.N[N]['y'])
        if AllNodes.N[N]['node_type']=='midmarker':
            NODE=N
            break
    if NODE!=-1:
        CenterX=AllNodes.N[NODE]['x']
        CenterY=AllNodes.N[NODE]['y']
    else:
        CenterX=(minx+maxx)/2
        CenterY=(miny+maxy)/2
    #TranslateNodes(NODEIDLIST,AllNodes,-CenterX,-CenterY)
    for N in NODELIST:
        News=ROTATE_POINT((AllNodes.N[N]['x'],AllNodes.N[N]['y']),angle,(CenterX,CenterY))
        AllNodes.N[N]['x']=News[0]
        AllNodes.N[N]['y']=News[1]
        if AllNodes.N[N]['node_type']=='metabolite':
            Labs=ROTATE_POINT((AllNodes.N[N]['label_x'],AllNodes.N[N]['label_y']),angle,(CenterX,CenterY))
            AllNodes.N[N]['label_x']=Labs[0]
            AllNodes.N[N]['label_y']=Labs[1]
    RL=ROTATE_POINT((Reaction.labx,Reaction.laby),angle,(CenterX,CenterY))
    Reaction.labx=RL[0]
    Reaction.laby=RL[1]
    for S in Reaction.seg:
        if Reaction.seg[S]['b1']!=None and Reaction.seg[S]['b1']!='null':
            B=ROTATE_POINT((Reaction.seg[S]['b1']['x'],Reaction.seg[S]['b1']['y']),angle,(CenterX,CenterY))
            Reaction.seg[S]['b1']['x']=B[0]
            Reaction.seg[S]['b1']['y']=B[1]
            B=ROTATE_POINT((Reaction.seg[S]['b2']['x'],Reaction.seg[S]['b2']['y']),angle,(CenterX,CenterY))
            Reaction.seg[S]['b2']['x']=B[0]
            Reaction.seg[S]['b2']['y']=B[1]

    #TranslateNodes(NODEIDLIST,AllNodes,CenterX,CenterY)
    
def TREEPLOT(ConsideredReactions,Links,AllNodes,MetaboliteID_To_Name,Reaction,NeverCombine):
    #returns [Dimensions,NodeIDs,ReactionIDs]
    #Dimensions are [width,heighth]
    REACTLIST=[Reaction[1]]    
    NODELIST=PLOT_REACTION_RELATIVE_TO_MIDMARKER(AllNodes,REACTLIST[0],MetaboliteID_To_Name)
    H=GETHEIGHT(NODELIST,AllNodes)
    W=GETWIDTH(NODELIST,AllNodes)
    #translate so that no node has negative Y coordinate
    TranslateNodes(NODELIST,AllNodes,-W[0],-H[0])
    TranslateReactions(REACTLIST,-W[0],-H[0])
    Height=H[1]-H[0]
    WIDTH=W[1]-W[0]
    TranslateNodes(NODELIST,AllNodes,-WIDTH/2,0)
    TranslateReactions(REACTLIST,-WIDTH/2,0) #now root is centered with center x=0
    ChildrenNodeLists=list()
    ChildrenReactionLists=list()
    ChildrenDimensions=list()

    FProdC=dict()
    FReaC=dict()
    PL=REACTLIST[0].met
    Pro=set()
    Rea=set()
    for M in PL:
        if M['bigg_id'] not in NeverCombine:
            if M['coefficient']>0:
                Pro.add(M['bigg_id'])
            if M['coefficient']<0:
                Rea.add(M['bigg_id'])
    MyMet=dict()




        
    for CID in Links[Reaction[0]]:
        Child=[CID,ConsideredReactions[CID]]
        BigList=TREEPLOT(ConsideredReactions,Links,AllNodes,MetaboliteID_To_Name,Child,NeverCombine)
        ChildrenNodeLists.append(BigList[1])
        ChildrenReactionLists.append(BigList[2])
        ChildrenDimensions.append(BigList[0])
        Cmet=ConsideredReactions[CID].met    
        ChList=set()
        for ChildMet in Cmet:
            ChList.add(ChildMet['bigg_id'])
        MyMet.update({CID:ChList})
    for P in Pro:
        PList=list()
        for R in dict(MyMet):
            if P in MyMet[R]:
                PList.append(R)
                del MyMet[R]
        if len(PList)>0:
            FProdC.update({P:PList})
    for P in Rea:
        PList=list()
        for R in dict(MyMet):
            if P in MyMet[R]:
                PList.append(R)
                del MyMet[R]
        if len(PList)>0:
            FReaC.update({P:PList})
    #so FProdC is a dictionary, contains products and the children associated with
    #FReaC is a dict, contains reactants and the children associated with

    MetToNodeID=dict()
    for M in FProdC:
        for N in NODELIST:
            if(AllNodes.N[N]['node_type']=='metabolite'):
            #print(AllNodes.N[N])
                if AllNodes.N[N]['bigg_id']==M:
                    MetToNodeID.update({M:N})
                    break
    for M in FReaC:
        for N in NODELIST:
            if(AllNodes.N[N]['node_type']=='metabolite'):
                if AllNodes.N[N]['bigg_id']==M:
                    MetToNodeID.update({M:N})
                    break

    #want to do, when placing children, sort them by parent metabolite with reactants first then products
    #then move the associated metabolite to be closer to child cluster

    if len(Links[Reaction[0]])==1: #Same rules as before apply to things with one child, but already did legwork
        if len(FProdC)>0: #the child is of a product
            ROT=90
            ROTATE_BY_DEG_NODE(NODELIST,AllNodes,REACTLIST[0],ROT)
            JJ=WIDTH
            WIDHT=Height
            Height=WIDTH
        else: #the child is of a reactant
            ROT=-90
            ROTATE_BY_DEG_NODE(NODELIST,AllNodes,REACTLIST[0],ROT)
            JJ=WIDTH
            WIDHT=Height
            Height=WIDTH
            
    #OK, keep a dict, Met in Father with Child
    #arrange children grouping them by Met in Father
##    #move Met in Father to be near the group of children
##    if len(Links[Reaction[0]])==1:  #has one or no children
##        #check if met in common is product, reactant or both
##        PL=REACTLIST[0].met
##        MM=set()
##        for M in PL:
##            MM.add(M['bigg_id'])
##            
##        CL=ConsideredReactions[Links[Reaction[0]][0]].met
##        CC=set() 
##        for C in CL:
##            CC.add(C['bigg_id'])
##        Shist=set(NeverCombine)
##        Common=MM.intersection(CC)
##        Common=Common-Shist
##        Prod=0
##        React=0
##        for X in Common:
##            for M in PL:
##                if X==M['bigg_id']:
##                    if M['coefficient']>0:
##                        Prod=1
##                    if M['coefficient']<0:
##                        React=1
##        if Prod==1 and React==0:
##            ROT=90
##            ROTATE_BY_DEG_NODE(NODELIST,AllNodes,REACTLIST[0],ROT)
##            JJ=WIDTH
##            WIDHT=Height
##            Height=WIDTH
##        if React==1 and Prod==0:
##            ROT=-90
##            ROTATE_BY_DEG_NODE(NODELIST,AllNodes,REACTLIST[0],ROT)
##            JJ=WIDTH
##            WIDHT=Height
##            Height=WIDTH
        
        #ROTATE_BY_DEG_REACT(REACTLIST,ROT)
                    
            

    
    

    ChildWidth=0 #width of all children 
    for D in ChildrenDimensions:
        w=D[0]
        if ChildWidth==0:
            ChildWidth=ChildWidth+w
        else:
            ChildWidth=ChildWidth+w+250
    ChildWidth=max(ChildWidth,WIDTH)
    XdispRoot=ChildWidth/2
    TranslateNodes(NODELIST,AllNodes,XdispRoot,0)
    TranslateReactions(REACTLIST,XdispRoot,0) #center the root at middle of children
    count=0
    drawX=0
    
    MaxY=0


    #Links[Reaction[0]] list of children root reaction IDs
    #ChildrenDimensions list ordered same as above
    #ChildrenNodeLists nods ordered same as above
    #ChildrenReactionLists contains all reactions in child
    #FReaC dict containing reactants and associated children roots
    #FProdC dict containing products and associated children roots

    for M in FReaC:
        LeftMostPoint=drawX
        for Children in FReaC[M]:
            Index=Links[Reaction[0]].index(Children)
            Xdisp=ChildrenDimensions[Index][0]/2+drawX
            MaxY=max(MaxY,ChildrenDimensions[Index][1])
            drawX=drawX+ChildrenDimensions[Index][0]+250
            TranslateNodes(ChildrenNodeLists[Index],AllNodes,Xdisp,Height+350)
            NODELIST=NODELIST+ChildrenNodeLists[Index]
            TranslateReactions(ChildrenReactionLists[Index],Xdisp,Height+350)
            REACTLIST=REACTLIST+ChildrenReactionLists[Index]
        RightMostPoint=drawX-250
        #now translate metabolite M to the center of this cluster
        CenterCoord=(RightMostPoint+LeftMostPoint)/2
        OldCoord=AllNodes.N[MetToNodeID[M]]['x']
        TransX=CenterCoord-OldCoord
        TranslateNodes([MetToNodeID[M]],AllNodes,TransX,150)
    for M in FProdC:
        LeftMostPoint=drawX
        for Children in FProdC[M]:
            Index=Links[Reaction[0]].index(Children)
            Xdisp=ChildrenDimensions[Index][0]/2+drawX
            MaxY=max(MaxY,ChildrenDimensions[Index][1])
            drawX=drawX+ChildrenDimensions[Index][0]+250
            TranslateNodes(ChildrenNodeLists[Index],AllNodes,Xdisp,Height+350)
            NODELIST=NODELIST+ChildrenNodeLists[Index]
            TranslateReactions(ChildrenReactionLists[Index],Xdisp,Height+350)
            REACTLIST=REACTLIST+ChildrenReactionLists[Index]
        RightMostPoint=drawX-250
        #now translate metabolite M to the center of this cluster
        CenterCoord=(RightMostPoint+LeftMostPoint)/2
        OldCoord=AllNodes.N[MetToNodeID[M]]['x']
        TransX=CenterCoord-OldCoord
        TranslateNodes([MetToNodeID[M]],AllNodes,TransX,150)
    
##    while count<len(ChildrenNodeLists):#start left most point at 0, stack children left to right
##        #children are assumed to be their width and centered at 0, so
##        Xdisp=ChildrenDimensions[count][0]/2+drawX
##        if MaxY<ChildrenDimensions[count][1]:
##            MaxY=ChildrenDimensions[count][1]
##        #Ydisp needs work of course
##        drawX=drawX+ChildrenDimensions[count][0]+250
##        TranslateNodes(ChildrenNodeLists[count],AllNodes,Xdisp,Height+250)
##        NODELIST=NODELIST+ChildrenNodeLists[count]
##        TranslateReactions(ChildrenReactionLists[count],Xdisp,Height+250)
##        REACTLIST=REACTLIST+ChildrenReactionLists[count]
##        count=count+1
    TotalH=Height+MaxY
    if ChildWidth>WIDTH:
        K=ChildWidth
    else:
        K=WIDTH
    TranslateNodes(NODELIST,AllNodes,-K/2,0)
    TranslateReactions(REACTLIST,-K/2,0)
    return [[K,TotalH],NODELIST,REACTLIST]
def GETWIDTH(NODELIST,AllNodes):
    xmax=0
    xmin=1000000
    for ID in NODELIST:
        if xmax<AllNodes.N[ID]['x']:
            xmax=AllNodes.N[ID]['x']
        if xmin>AllNodes.N[ID]['x']:
            xmin=AllNodes.N[ID]['x']
    return([xmin,xmax])

def GETHEIGHT(NODELIST,AllNodes):
    ymax=0
    ymin=1000000
    for ID in NODELIST:
        if ymax<AllNodes.N[ID]['y']:
            ymax=AllNodes.N[ID]['y']
        if ymin>AllNodes.N[ID]['y']:
            ymin=AllNodes.N[ID]['y']
    return([ymin,ymax])
        

def TranslateNodes(NODEIDLIST,NODES,Xdisp,Ydisp):
    for ID in NODEIDLIST:
        NODES.N[ID]['x']=NODES.N[ID]['x']+Xdisp
        NODES.N[ID]['y']=NODES.N[ID]['y']+Ydisp
        if NODES.N[ID]['node_type']=='metabolite':
            NODES.N[ID]['label_x']=NODES.N[ID]['label_x']+Xdisp
            NODES.N[ID]['label_y']=NODES.N[ID]['label_y']+Ydisp
def TranslateReactions(Reactions,Xdisp,Ydisp):
    for R in Reactions:
        R.labx=R.labx+Xdisp
        R.laby=R.laby+Ydisp
        for S in R.seg:
            if R.seg[S]['b1']!='null' and R.seg[S]['b1']!=None:
                R.seg[S]['b1']['x']=R.seg[S]['b1']['x']+Xdisp
                R.seg[S]['b2']['x']=R.seg[S]['b2']['x']+Xdisp
                R.seg[S]['b1']['y']=R.seg[S]['b1']['y']+Ydisp
                R.seg[S]['b2']['y']=R.seg[S]['b2']['y']+Ydisp

def COMBINE_TWO_MAPS(Map1,Map2):
    import json
    
    with open(Map1,'r') as f:
        First=json.load(f)
    with open(Map2,'r') as f:
        Second=json.load(f)   
    Head=First[0]

    Head['map_name']=Head['map_name']+'_'+Second[0]['map_name']
    FirstNodeID=0
    FirstRID=0
    FirstSegID=0
    MaxX=0
    MinX=0
    for n in First[1]['nodes']:
        FirstNodeID=max(FirstNodeID,int(n))
        MaxX=max(float(First[1]['nodes'][n]['x']),MaxX)
        MinX=min(float(First[1]['nodes'][n]['x']),MinX)
    for r in First[1]['reactions']:
        FirstRID=max(FirstRID,int(r))
        for s in First[1]['reactions'][r]['segments']:
            FirstSegID=max(FirstSegID,int(s))
    XDisp=MaxX-MinX+500
    NewN=FirstNodeID+1
    NewR=FirstRID+1
    NewS=FirstSegID+1
    NewNodes=dict()
    for n in Second[1]['nodes']:
        NODE=Second[1]['nodes'][n]
        NODE['x']=float(NODE['x'])+XDisp
        if NODE['node_type']=='metabolite':
            NODE['label_x']=float(NODE['label_x'])+XDisp
        for r in Second[1]['reactions']:
            for s in Second[1]['reactions'][r]['segments']:
                if str(Second[1]['reactions'][r]['segments'][s]['from_node_id'])==str(n):
                    Second[1]['reactions'][r]['segments'][s]['from_node_id']=NewN
                if str(Second[1]['reactions'][r]['segments'][s]['to_node_id'])==str(n):
                    Second[1]['reactions'][r]['segments'][s]['to_node_id']=NewN                
        NewNodes.update({NewN:NODE})
        NewN=NewN+1
    NewReactions=dict()
    for r in Second[1]['reactions']:
        Reaction=Second[1]['reactions'][r]
        Reaction['label_x']=float(Reaction['label_x'])+XDisp
        NewSegments=dict()
        for s in Reaction['segments']:
            Seg=Reaction['segments'][s]
            if(Seg['b1']!=None)and(Seg['b1']!='null'):
                Seg['b1']['x']=float(Seg['b1']['x'])+XDisp
                Seg['b2']['x']=float(Seg['b2']['x'])+XDisp
            NewSegments.update({NewS:Seg})
            NewS=NewS+1
        Reaction.update({'segments':NewSegments})
        NewReactions.update({NewR:Reaction})
        NewR=NewR+1
    for i in NewNodes:
        First[1]['nodes'].update({i:NewNodes[i]})
    for r in NewReactions:
        First[1]['reactions'].update({r:NewReactions[r]})
    First[1]['canvas']['width']=First[1]['canvas']['width']+Second[1]['canvas']['width']
    First[1]['canvas']['height']=max(First[1]['canvas']['height'],Second[1]['canvas']['height'])
    PrintString=json.dumps(First)
    open(Map1.replace('.json','+')+Map2,'w').write(PrintString)


def REPORTSUBSYSTEM(Cobra_Model):
    import json
    with open(Cobra_Model,'r') as f:
        M=json.load(f)
    SUBSYSTEMS=dict()
    for R in M['reactions']:
        Sub=R['subsystem']
        Col=set()
        if Sub.split(':')[0] in SUBSYSTEMS:
            Col=SUBSYSTEMS[Sub.split(':')[0]]
        if Sub!=Sub.split(':')[0]:
            Col.add(Sub)
        SUBSYSTEMS.update({Sub.split(':')[0]:Col})
    PrintList=list()
    for S in SUBSYSTEMS:
        PrintList.append(S)
        for H in SUBSYSTEMS[S]:
            PrintList.append('\t'+H)
    open('CurrentProgressPathways.txt','w').write('\n'.join(PrintList))

def SPLIT_MODEL_INTO_MODELS_BY_SUBSYSTEM(Cobra_Model):
    import json
    import os

    d = Cobra_Model.replace('.json','_Subsystems/')
    if not os.path.exists(d):
        os.makedirs(d)
    with open(Cobra_Model,'r') as f:
        M=json.load(f)
    SUBSYSTEMS=set()
    for R in M['reactions']:
        SUBSYSTEMS.add(R['subsystem'].split(':')[0].split(',')[0])
    FILE=dict(M)
    for Sub in SUBSYSTEMS:
        LocalReacts=list()
        for R in M['reactions']:
            if R['subsystem'].split(':')[0].split(',')[0]==Sub:
                LocalReacts.append(R)
        FILE.update({'reactions':LocalReacts})
        open(d+Sub.replace('/','').replace('.','').replace(' ','').replace(':','')+'.json','w').write(json.dumps(FILE))

    #open('SubsystemReport.txt','w').write('\n'.join(list(SUBSYSTEMS)))

def MAKE_ALL_FROM_FOLDER(InFolder,OutFolder,Oldmodel=''):
    import os
    if not os.path.exists(OutFolder):
        os.makedirs(OutFolder)
    for File in os.listdir(InFolder):
        print(File)
        OutName=OutFolder+'/'+File
        GROWFROMSINGLEREACTION(InFolder+'/'+File,'',8,12,[],OutName.replace('.json',''),'All',Oldmodel)


def GEN_SPLITTED_REPORTS(Cobra_Model):
    import json
    import os
    import csv

    EndTrans=dict()
    with open('pt3topt2anno_BBH.csv', 'rt') as csvfile:
        RReater = csv.reader(csvfile, delimiter=',')
        count=0
        for row in RReater:
            if count==0:
                count=1
            else:
                pt3=row[0]
                KEGG=row[8]
                EndTrans.update({pt3:KEGG})

#OK, so each subssystem is a folder
    #folder contains Json model of reactions only
    #list of reactions
    #csv of present gene IDs and the KEGGID of that gene
    #list of metabolites
    #list of metabolite BiGG IDS
    d = Cobra_Model.replace('.json','_Reports/')
    if not os.path.exists(d):
        os.makedirs(d)
    with open(Cobra_Model,'r') as f:
        M=json.load(f)
    MetCodeName=dict()
    for m in M['metabolites']:
        MetCodeName.update({m['id']:m['name']})
    SUBSYSTEMS=set()
    for R in M['reactions']:
        #SUBSYSTEMS.add(R['subsystem'].split(':')[0].split(',')[0])
        SUBSYSTEMS.add(R['subsystem'])
    FILE=dict(M)


    
    for Sub in SUBSYSTEMS:
        LocFoldLocation=d+Sub.replace('/','').replace('.','').replace(' ','').replace(':','')+'/'
        if not os.path.exists(LocFoldLocation):
            os.makedirs(LocFoldLocation)
        LocalReacts=list()

        REACTNAMES=list();
        MetCodes=set();
        MetNames=set();
        GENRULES=set();
        REACTNAMES.append('"Reaction BiGG ID","Reaction Name"')
        for R in M['reactions']:
            #if R['subsystem'].split(':')[0].split(',')[0]==Sub:
            if R['subsystem']==Sub:
                #print(R)
                LocalReacts.append(R)
                REACTNAMES.append('"'+R['id']+'","'+R['name']+'"')
                for m in R['metabolites']:
                    MetCodes.add(m)
                    MetNames.add(MetCodeName[m])
                GRULE=R['gene_reaction_rule'].replace('(','').replace(')','').replace(' and ',';').replace(' or ',';').split(';')
                GENRULES.update(GRULE)
        MetReportDict=dict()
        for m in MetCodes:
            if MetCodeName[m] not in MetReportDict:
                MetReportDict.update({MetCodeName[m]:[m]})
            else:
                OlList=MetReportDict[MetCodeName[m]]
                OlList.append(m)
                MetReportDict.update({MetCodeName[m]:OlList})
        MetReport=['"Metabolite Name","Metabolite ID Codes"']
        for m in MetReportDict:
            MetReport.append('"'+m+'","'+','.join(MetReportDict[m])+'"')
        GenReport=['"pt3_ID","KeGG_ID","LINK"']
        for g in GENRULES:
            if g!='':
                if g in EndTrans:
                    if EndTrans[g]!='':
                        GenReport.append('"'+g+'","'+EndTrans[g]+'",http://www.genome.jp/dbget-bin/www_bget?'+EndTrans[g])
                    else:
                        GenReport.append('"'+g+'","KEGG NOT FOUND"')
                elif 'Phatr_' in g:
                    num=g.replace('Phatr_','')
                    GenReport.append('"'+g+'","pti:PHATRDRAFT_'+num+'",http://www.genome.jp/dbget-bin/www_bget?pti:PHATRDRAFT_'+num)

                else:
                    GenReport.append('"'+g+'","KEGG NOT FOUND"')
                
        FILE.update({'reactions':LocalReacts})
        open(LocFoldLocation+Sub.replace('/','').replace('.','').replace(' ','').replace(':','').replace(',','')+'Model.json','w').write(json.dumps(FILE))
        open(LocFoldLocation+Sub.replace('/','').replace('.','').replace(' ','').replace(':','').replace(',','')+'_ReactionReport.csv','w').write('\n'.join(REACTNAMES))
        open(LocFoldLocation+Sub.replace('/','').replace('.','').replace(' ','').replace(':','').replace(',','')+'_MetaboliteReport.csv','w').write('\n'.join(MetReport))
        open(LocFoldLocation+Sub.replace('/','').replace('.','').replace(' ','').replace(':','').replace(',','')+'_GeneReport.csv','w').write('\n'.join(GenReport))
        GROWFROMSINGLEREACTION(LocFoldLocation+Sub.replace('/','').replace('.','').replace(' ','').replace(',','').replace(':','')+'Model.json','',8,12,[],LocFoldLocation+Sub.replace('/','').replace('.','').replace(' ','').replace(':','')+'SampleMap.json','All','Pti_manual_curated_balanced_09092015.json')

        print(Sub)

        
def GEN_SPLITTED_REPORTS_TARGET(Cobra_Model,OnlyThis):
    import json
    import os
    import csv

    EndTrans=dict()
    with open('pt3topt2anno_BBH.csv', 'rt') as csvfile:
        RReater = csv.reader(csvfile, delimiter=',')
        count=0
        for row in RReater:
            if count==0:
                count=1
            else:
                pt3=row[0]
                KEGG=row[8]
                EndTrans.update({pt3:KEGG})

#OK, so each subssystem is a folder
    #folder contains Json model of reactions only
    #list of reactions
    #csv of present gene IDs and the KEGGID of that gene
    #list of metabolites
    #list of metabolite BiGG IDS
    d = Cobra_Model.replace('.json','_Reports/')
    if not os.path.exists(d):
        os.makedirs(d)
    with open(Cobra_Model,'r') as f:
        M=json.load(f)
    MetCodeName=dict()
    for m in M['metabolites']:
        MetCodeName.update({m['id']:m['name']})
    SUBSYSTEMS=set()
    for R in M['reactions']:
        SUBSYSTEMS.add(R['subsystem'])
    FILE=dict(M)


    
    for Sub in SUBSYSTEMS:
        #print(Sub)
        if Sub.split(':')[0].split(',')[0].replace(' ','')==OnlyThis:
            LocFoldLocation=d+Sub.replace('/','').replace('.','').replace(' ','').replace(':','')+'/'
            if not os.path.exists(LocFoldLocation):
                os.makedirs(LocFoldLocation)
            LocalReacts=list()

            REACTNAMES=list();
            MetCodes=set();
            MetNames=set();
            GENRULES=set();
            REACTNAMES.append('"Reaction BiGG ID","Reaction Name"')
            for R in M['reactions']:
                if R['subsystem']==Sub:
                    #print(R)
                    LocalReacts.append(R)
                    REACTNAMES.append('"'+R['id']+'","'+R['name']+'"')
                    for m in R['metabolites']:
                        MetCodes.add(m)
                        MetNames.add(MetCodeName[m])
                    GRULE=R['gene_reaction_rule'].replace('(','').replace(')','').replace(' and ',';').replace(' or ',';').split(';')
                    GENRULES.update(GRULE)
            MetReportDict=dict()
            for m in MetCodes:
                if MetCodeName[m] not in MetReportDict:
                    MetReportDict.update({MetCodeName[m]:[m]})
                else:
                    OlList=MetReportDict[MetCodeName[m]]
                    OlList.append(m)
                    MetReportDict.update({MetCodeName[m]:OlList})
            MetReport=['"Metabolite Name","Metabolite ID Codes"']
            for m in MetReportDict:
                MetReport.append('"'+m+'","'+','.join(MetReportDict[m])+'"')
            GenReport=['"pt3_ID","KeGG_ID","LINK"']
            for g in GENRULES:
                if g!='':
                    if g in EndTrans:
                        if EndTrans[g]!='':
                            GenReport.append('"'+g+'","'+EndTrans[g]+'",http://www.genome.jp/dbget-bin/www_bget?'+EndTrans[g])
                        else:
                            GenReport.append('"'+g+'","KEGG NOT FOUND"')
                    elif 'Phatr_' in g:
                        num=g.replace('Phatr_','')
                        GenReport.append('"'+g+'","pti:PHATRDRAFT_'+num+'",http://www.genome.jp/dbget-bin/www_bget?pti:PHATRDRAFT_'+num)

                    else:
                        GenReport.append('"'+g+'","KEGG NOT FOUND"')
                    
            FILE.update({'reactions':LocalReacts})
            open(LocFoldLocation+Sub.replace('/','').replace('.','').replace(' ','').replace(':','')+'Model.json','w').write(json.dumps(FILE))
            open(LocFoldLocation+Sub.replace('/','').replace('.','').replace(' ','').replace(':','')+'_ReactionReport.csv','w').write('\n'.join(REACTNAMES))
            open(LocFoldLocation+Sub.replace('/','').replace('.','').replace(' ','').replace(':','')+'_MetaboliteReport.csv','w').write('\n'.join(MetReport))
            open(LocFoldLocation+Sub.replace('/','').replace('.','').replace(' ','').replace(':','')+'_GeneReport.csv','w').write('\n'.join(GenReport))
            GROWFROMSINGLEREACTION(LocFoldLocation+Sub.replace('/','').replace('.','').replace(' ','').replace(':','')+'Model.json','',2,6,[],LocFoldLocation+Sub.replace('/','').replace('.','').replace(' ','').replace(':','')+'SampleMap.json','All','Pti_manual_curated_balanced_09092015.json')
            print(Sub)



    
    #open('SubsystemReport.txt','w').write('\n'.join(list(SUBSYSTEMS)))
    
#GROWFROMSINGLEREACTION('Pti_manual_curated_balanced_09092015.json','ECH121_2E_m',2)
#GROWFROMSINGLEEACTION('Pti_manual_curated_balanced_09092015.json','biomass_carb_c',2)
#GROWFROMSINGLEREACTION('Pti_manual_curated_balanced_09092015.json','EX_photon_e',6)
#GROWFROMSINGLEREACTION('Pti_manual_curated_balanced_09092015.json','PHYPQOX_h',3)
#GROWFROMSINGLEREACTION('Pti_manual_curated_balanced_09092015.json','ACCOAC_c',5)
#GROWFROMSINGLEREACTION('Pti_manual_curated_balanced_09092015.json','BCTAL_c',2)
#GROWFROMSINGLEREACTION('Pti_manual_curated_balanced_09092015.json','G6PDH_c',3)
#GROWFROMSINGLEREACTION('Pti_manual_curated_balanced_09092015.json','biomass_pro_c',2)
#GROWFROMSINGLEREACTION('Pti_manual_curated_balanced_09092015.json','2MBCOAD_m',2)
#GROWFROMSINGLEREACTION('Pti_manual_curated_balanced_09092015.json','MALOAAt_m',2,6)
#GROWFROMSINGLEREACTION('Pti_manual_curated_balanced_09092015.json','ACCOAC_c',2,5)

    
BIOMASSIGNORE=['biomass_pro_c','biomass_pigm_h','biomass_DNA_c','biomass_mem_lipids_c','biomass_TAG_c','bof_c']
#BIOMASSIGNORE=list()

#GROWFROMSINGLEREACTION('Pti_manual_curated_balanced_09092015.json','G6PDH_c',3,6,BIOMASSIGNORE)
#GROWFROMSINGLEREACTION('Pti_manual_curated_balanced_09092015.json','NOR_c',3,6,BIOMASSIGNORE)
#GROWFROMSINGLEREACTION('Pti_manual_curated_balanced_09092015.json','ACCOAC_c',2,6,BIOMASSIGNORE)
#GROWFROMSINGLEREACTION('NucleotidemetabolismModel.json','',2,6,'','Nucleotidemetabolism_LessConnected','All','Pti_manual_curated_balanced_09092015.json','yes')
#GROWFROMSINGLEREACTION('N-GlycanbiosynthesisModel.json','',2,6,'DOLK_c','N-Glycanbiosynthesis_MAP_SIMPLE','','Pti_manual_curated_balanced_09092015.json')



#GROWFROMSINGLEREACTION('TriacylglycerolbiosynthesisModel.json','',3,6,'','TriacylglycerolbiosynthesisModel_MAP_SIMPLE','','Pti_manual_curated_balanced_09092015.json') #Do lipid
GROWFROMSINGLEREACTION('FattyaciddegradationModel.json','',3,6,'','Fattyaciddegradation_map_2_For_Alex','','Pti_manual_curated_balanced_09092015.json') #Do lipid


#REPORTSUBSYSTEM('Pti_manual_curated_balanced_09092015.json')
#SPLIT_MODEL_INTO_MODELS_BY_SUBSYSTEM('Pti_manual_curated_balanced_09092015.json')
#MAKE_ALL_FROM_FOLDER('Pti_manual_curated_balanced_09092015_Subsystems','Pti_manual_curated_balanced_09092015_MAPS','Pti_manual_curated_balanced_09092015.json')
#GEN_SPLITTED_REPORTS('Pti_manual_curated_balanced_09092015.json');
#GEN_SPLITTED_REPORTS_TARGET('Pti_manual_curated_balanced_09092015.json','Nucleotidemetabolism')

#GROW_FROM_EXTERNAL_LIKE_CRYSTALS('iJN678_Fixed_Names.json')
#COMBINE_TWO_MAPS('Tester111.json','Tester222.json')
