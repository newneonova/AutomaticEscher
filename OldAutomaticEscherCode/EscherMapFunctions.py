#This file contains three functions which read Escher map json files and works
#with them.




#

def GATHER_BIGG_DESCRIPTIVE_NAMES(Cobra_Model):
    import json
    with open(Cobra_Model,'r') as f:
        M=json.load(f)
    MetList=M['metabolites']
    import time #the BiGG database does not want more than 10 requests per second, so slow it down
    #also for each model, save a ID/Name reference file
    #when running, look for file, if file doesn't exist make a new one
    #that level of logic in upper function, this one makes the file
    import urllib.request
    ROOT_ADDRESS='http://bigg.ucsd.edu/api/v2/universal/metabolites/'
    SaveMets=dict()
    print('Fetching from BiGG Database, limited to 10 per second')
    count=0
    for i in MetList:
        BiGG_ID=i['id'][0:len(i['id'])-2]
        HUH=json.loads(urllib.request.urlopen(ROOT_ADDRESS+BiGG_ID).read().decode("utf-8"))
        SaveMets.update({BiGG_ID+'_c':HUH['name']})
        count+=1
        time.sleep(.1)
        if count%10==0:
            print(str(count)+'/'+str(len(MetList)))
    print('Saving BiGG ID/Name reference exchange in "'+Cobra_Model.replace('.json','_BiGG_NAME_Metabolites.json')+'"')
    open(Cobra_Model.replace('.json','_BiGG_NAME_Metabolites.json'), 'w').write(json.dumps(SaveMets))





def COMBINE_LIKE_NODES_WITH_EXCLUSIONS(Escher_map,Exclusion_list=None):
    #COMBINE_LIKE_NODES_WITH_EXCLUSIONS takes two arguments, the address of the
    #Escher_map as name.json if in the same direcotry as the running script, and
    #the list of excluded metabolites. The exclusion set is a list of BiGG IDs
    #for the metabolites that we do not want to combine
    if Exclusion_list is None:
        Exclusion_list = []
    Exclusion=set(Exclusion_list)
    import json
    print('loading Escher Map model '+Escher_map)
    with open(Escher_map,'r') as f:
        M=json.load(f)
    OutName=str(M[0].get('map_name')).replace(' ','_')+'_Combined_Metabolites.json'
    Header=M[0]
    text_labels=M[1].get('text_labels')
    canvas=M[1].get('canvas')
    StartReaction=M[1].get('reactions')
    StartNodes=M[1].get('nodes')
    ReplaceList=list() #of the form ((keepid,replaceid),...)
    MetNodes=dict()
    KeepNodes=dict()
    print('Keeping unique nodes')
    for i in StartNodes:
        ThisNode=StartNodes.get(i)
        if (ThisNode.get('node_type')!='metabolite') or (StartNodes.get(i).get('bigg_id') in Exclusion) :
            KeepNodes.update({i:ThisNode})
        elif StartNodes.get(i).get('bigg_id') not in MetNodes:
            MetNodes.update({StartNodes.get(i).get('bigg_id'):i})
            KeepNodes.update({i:ThisNode})
        else:
            ReplaceList.append(str(MetNodes.get(StartNodes.get(i).get('bigg_id')))+','+str(i))
            #print(str(MetNodes.get(StartNodes.get(i).get('bigg_id')))+','+str(i))
    ReactString=json.dumps(StartReaction)
    print('Redirecting old paths to new nodes')
    for i in ReplaceList:
        Keep=i.split(',')[0]
        Toss=i.split(',')[1]
        ReactString=ReactString.replace('id": "'+Toss,'id": "'+Keep)
        ReactString=ReactString.replace('id": '+Toss,'id": '+Keep)
    EndReaction=json.loads(ReactString)
    Rest=dict({'reactions':EndReaction})
    Rest.update({'canvas':canvas})
    Rest.update({'text_labels':text_labels})
    Rest.update({'nodes':KeepNodes})
    Whole=dict(Header)
    Whole.update(Rest)
    PrintString='['+json.dumps(Header)+','+json.dumps(Rest)+']'
    open(OutName,'w').write(PrintString)
    print('Done')

def CHECK_FOR_MISSING_METABOLITES(Escher_map,Cobra_file):     
    #This program examines the reactions present in the cobra file and checks to see if the reactions
    #in the map have the same metabolites, This is to make sure that in the mapping process no metabolie
    #was accidently deleted.
    #this results in a text file listing all reactions which are missing metabolites.
    import json
    print('Loading Escher Map')
    with open(Escher_map,'r') as f:
        M=json.load(f)
    OutName=str(M[0].get('map_name')).replace(' ','_')+'_Bad_metabolites.txt'
    print('Loading COBRA file')
    with open(Cobra_file,'r') as f:
        C=json.load(f)

    CR=C.get('reactions')
    CobraDict=dict()
    for i in CR:
        ID=i.get('id')
        CobraDict.update({ID:i.get('metabolites')})
    MR=M[1].get('reactions')
    MN=M[1].get('nodes')
    BadList=list()
    for i in MR:
        REACT=MR.get(i)
        ID=REACT.get('bigg_id')
        SEG=REACT.get('segments')
        nodeset=set()
        for s in SEG:
            nodeset.add(SEG.get(s).get('from_node_id'))
            nodeset.add(SEG.get(s).get('to_node_id'))
        MetSet=set()
        for n in list(nodeset):
            NO=MN.get(n)
            if NO.get('node_type')=='metabolite':
                MetSet.add(NO.get('bigg_id'))
        for Metabol in CobraDict.get(ID):
            if Metabol not in MetSet:
                BadList.append(ID)
                break
    print('There are '+str(len(BadList))+' reactions with incorrect metabolites')
    if len(BadList)>0:
        open(OutName,'w').write('\n'.join(BadList))
    print('Done')



#The following pair of functions, PLACE_SUB_REACTION and PLACE_CORE_REACTION are used in
#AUTO_MAP_LINKED to group bunches of reactions into simple trees for plotting on to Escher maps
#To make conversion from this grouping function to actual Escher notation, we need to generate more information
#and pass the functions more information.
#Namely, we will also pass the functions COBRA model reaction information so we can add all the secondary
#metabolite nodes.  Then on each level we need to generate a Dict labeled NODES and one SEGMENTS
#NODES contains all the node info for the reaction on that level, primary, secondary, mid and multi marker
#Segments contian all the information on connections from node to node

def MMAAXX(MaxAndMin,Test):
    import math
    return((max(MaxAndMin[0],Test[0]),max(MaxAndMin[1],Test[1]),min(MaxAndMin[2],Test[0]),min(MaxAndMin[3],Test[1])))



def MAKE_NODES(MaxAndMin,reaction,center,angle,primary_metabolites,U,NameKey,ReactMet=''):
    import math

    LevelInfo=dict()
    NodeDic=dict()
    G=0
    NodeDic.update({0:{'node_type':'midmarker','x':center[0],'y':center[1]}})
    G+=1
    
    if len(primary_metabolites)==2:
        J=list(primary_metabolites)
        if abs(reaction.get('metabolites').get(J[0])+reaction.get('metabolites').get(J[1]))==2:
            angle=angle+math.pi/2
 
    if ReactMet!='': #then was called as root reaction
        Source_Coeff=reaction.get('metabolites').get(ReactMet)
    elif len(primary_metabolites)==1:
        Source_Coeff=-reaction.get('metabolites').get(list(primary_metabolites)[0]) #pretend we came from opposite direction of only child metabolite
        
    elif len(primary_metabolites)==2:
        #assign one of the metabolites as the source coefficient, assign the one with negative y coord as source (because this happens only when we are center, with two its verticle)
        if primary_metabolites[list(primary_metabolites)[0]][1]>1:
            Source_Coeff=reaction.get('metabolites').get(list(primary_metabolites)[1])
        else:
            Source_Coeff=reaction.get('metabolites').get(list(primary_metabolites)[0])
        
    else:
        Source_Coeff=1
    if Source_Coeff > 0: #then linking metabolite is a reactant
        Dir=-1
        NodeDic.update({1:{'node_type':'multimarker','x':(center[0]+Dir*U/10*math.cos(angle)),'y':(center[1]+Dir*U/10*math.sin(angle))}}) #Node 1 is toward products/links to products
        G+=1
        NodeDic.update({2:{'node_type':'multimarker','x':(center[0]-Dir*U/10*math.cos(angle)),'y':(center[1]-Dir*U/10*math.sin(angle))}}) #Node 2 is toward reactants
        
    else: #linking metabolite is a product
        Dir=1
        NodeDic.update({1:{'node_type':'multimarker','x':(center[0]+Dir*U/10*math.cos(angle)),'y':(center[1]+Dir*U/10*math.sin(angle))}})#Node 1 is toward products/links to products
        G+=1
        NodeDic.update({2:{'node_type':'multimarker','x':(center[0]-Dir*U/10*math.cos(angle)),'y':(center[1]-Dir*U/10*math.sin(angle))}})#Node 2 is toward reactants
        

                       




    countPos=2
    countNeg=2
    #print(reaction)
    for Guy in list(reaction.get('metabolites')):
        #print(Guy)
        G+=1
        COEFF=float(reaction.get('metabolites').get(Guy))
        if Guy in primary_metabolites:
            COR=primary_metabolites.get(Guy)
            Prim='true'
        else:

            if COEFF<0:
                Side=(((countNeg%2)-.5)/.5)
                Slot=math.floor(countNeg/2)
                D=Side*Slot*100
                APPLE=angle+math.pi/2+Dir*Side*Slot*math.pi/10
                CorBaseNode=(NodeDic[2]['x'],NodeDic[2]['y'])
                COR=(CorBaseNode[0]+D*math.cos(APPLE),CorBaseNode[1]+D*math.sin(APPLE))
                
                
                countNeg+=1
            elif COEFF>0:
                Side=(((countPos%2)-.5)/.5)
                Slot=math.floor(countPos/2)
                D=Side*Slot*100
                APPLE=angle+math.pi/2-Dir*Side*Slot*math.pi/10
                CorBaseNode=(NodeDic[1]['x'],NodeDic[1]['y'])
                COR=(CorBaseNode[0]+D*math.cos(APPLE),CorBaseNode[1]+D*math.sin(APPLE))
                countPos+=1
            Prim='false'
        N=Guy
        try:
            Name=NameKey[Guy]
        except KeyError:
            Name=N
        if COEFF<0:
            NodeNum=2
        else:
            NodeNum=1
        NX=NodeDic.get(NodeNum).get('x')
        NY=NodeDic.get(NodeNum).get('y')
        OY=NodeDic.get(0).get('y')
        OX=NodeDic.get(0).get('x')
        LCor=(COR[0]+(NX-OX),COR[1]+(NY-OY))
        NodeDic.update({G:{'node_type':'metabolite','x':COR[0],'y':COR[1],'label_x':LCor[0],'label_y':LCor[1],'name':Name,'bigg_id':N,'node_is_primary':Prim}})


    LevelInfo.update({'nodes':NodeDic})
    LevelInfo.update({'label_x':(center[0]-U/5*math.cos(angle+math.pi/2))})
    LevelInfo.update({'label_y':(center[1]-U/5*math.sin(angle+math.pi/2))})
    SegDic=dict()
    H=0
    SegDic.update({H:{'from_node_id':1,'to_node_id':0,'b1':'null','b2':'null'}})
    H+=1
    SegDic.update({H:{'from_node_id':2,'to_node_id':0,'b1':'null','b2':'null'}})
    H+=1
    #rather than go from - to + and work on the multi marker, connect directly to midmarker
    #b1 is close to the from node in the direction of the to node
    #b2 is close to the to node in the direction of the from node

    #make node 1 the +, make node 2 the - node. to and from do not seem to matter.
    for j in range(3,G+1):
        N_T=NodeDic.get(j)
        name=N_T.get('bigg_id')
        COEFF=reaction.get('metabolites').get(name)
        if COEFF > 0:
            NodeNum=1
        else:
            NodeNum=2
        NX=NodeDic.get(NodeNum).get('x')
        NY=NodeDic.get(NodeNum).get('y')
        OY=NodeDic.get(0).get('y')
        OX=NodeDic.get(0).get('x')
        ManX=NodeDic.get(j).get('x')
        ManY=NodeDic.get(j).get('y')
        D=math.sqrt(pow((N_T.get('x')-NX),2)+pow((N_T.get('y')-NY),2))
        #KAPPAONE={'x':(NX+D/10*(N_T.get('x')-NX)),'y':(NY+D/10*(N_T.get('y')-NY))}
        #KAPPATWO={'x':(N_T.get('x')+D/10*(NX-N_T.get('x'))),'y':(N_T.get('y')+D/10*(NY-N_T.get('y')))}
        #KAPPAONE={'x':NX+1,'y':NY+1}
        #KAPPATWO={'x':N_T.get('x')-1,'y':N_T.get('y')-1}
        KAPPAONE={'x':NX-2*(OX-NX),'y':NY-2*(OY-NY)}
        KAPPATWO={'x':ManX-(ManX-(NX-2*(OX-NX)))/4,'y':ManY-(ManY-(NY-2*(OY-NY)))/4}        
        
        SegDic.update({H:{'from_node_id':NodeNum,'to_node_id':j,'b1':KAPPAONE,'b2':KAPPATWO}})
        H+=1
    LevelInfo.update({'segments':SegDic})
    for i in NodeDic:
        MaxAndMin=MMAAXX(MaxAndMin,(NodeDic.get(i).get('x'),NodeDic.get(i).get('y')))
    if False and reaction['id']=='PGSA160':
        print(primary_metabolites)
        print(Dir)
        print(Source_Coeff)
        print(NodeDic)
        print(ReactMet=='')
        import sys
        #sys.exit(0)
    return([LevelInfo,MaxAndMin])
       
#Now have a dict, reactions and their neighbors by connection

    #Now use this and the adj counter information to start grouping them. Also start positional
    #Layout by graph theory

    #Start a Sub_Map
    #Choose the reaction with the most connections, delete it from counter
    #place it in the center of virtural space. Have growths sticking off of it.
    #Angle of the growths is 360 degrees divided by number of connections
    #first growth at 0, nth at (n-1)360/n.

    #Alg works as follows.
        #Load Dic entry for Reaction, Delete from DIc (WIll ignore loops)
        #for each growth direction:
            #Pop tuple from Dic Entry
            #Stick node for metabolite at ANGLE and distance 1 from reaction
            #if New reaction in Dic, proceed, else is done with this growth.
            #stick new reaction at ANGLE and distance 2 from Core reaction
            #Resolve new reaction given CoreName and ANGLE:
                #New reaction like core reaction except connection 1 already done
                #remove entry from count, load dic entry and delete it from Dic
                #Delete from DIC list CoreName tuple. if no tuples in list then done
                #else array growths as n=NumInTupleList+1
                #First growth is skipped
                #2 to n growths at ANGLENEW=(n-1)360/n+ANGLE-180 MOD 360
                #For each growth, Metabolite at ANGLENEW distance 1
                #save metabolite in list
                # if NEWREACTION in DIC, NEWREACT at ANGLENEW distance 2
                #Call this function with REACTNAME, ANGLENEW
                #If no more growth nodes,populate secondary metabolites, distance 1/2
                    #if metabolite not in saved list, place node at 10 deg increments from coreangle                      

def PLACE_SUB_REACTION(MaxAndMin,Upname,U,Relations,Cou,InAngle,NewCenter,ReactPair,OMetCoord,REACTIONSOURCE,Metabolites):
    import math
    ReactName=ReactPair[0]
    #print(Upname)
    #print(ReactName)
    ReactMet=ReactPair[1]
    #print(Upname+', '+ReactName)
    BabyInfo=dict()
    #BabyInfo.update({'React_Name':ReactName})
    #BabyInfo.update({'Uplink_Metabolite':ReactMet})
    #BabyInfo.update({'center':NewCenter})
    del Cou[ReactName]
    ReaRel=Relations.pop(ReactName)
    AdjCount=len(ReaRel)
    del ReaRel[ReaRel.index((Upname,ReactPair[1]))]
    PrimaryMetabolites=dict()
    SubReactions=list()
    PrimaryMetabolites.update({ReactPair[1]:OMetCoord})
    for NUM in range(1,AdjCount): #NUM =0 reserved for parent link
        NewAngle=math.radians((((NUM*360)/AdjCount)+math.degrees(InAngle)-180)%360)
        Current=ReaRel.pop()
        if Current[0] in Relations:
            MetCoord=(NewCenter[0]+(U*math.cos(NewAngle)),NewCenter[1]+(U*math.sin(NewAngle)))
            MaxAndMin=MMAAXX(MaxAndMin,MetCoord)
            NewRCent=(NewCenter[0]+(2*U*math.cos(NewAngle)),NewCenter[1]+(2*U*math.sin(NewAngle)))
            MaxAndMin=MMAAXX(MaxAndMin,NewRCent)
            PrimaryMetabolites.update({Current[1]:MetCoord})
            Loona=PLACE_SUB_REACTION(MaxAndMin,ReactName,U,Relations,Cou,NewAngle,NewRCent,Current,MetCoord,REACTIONSOURCE,Metabolites)
            SubReactions.append(Loona[0])
            MaxAndMin=Loona[1]
    #BabyInfo.update({'primay_metabolites':PrimaryMetabolites})
    #BabyInfo.update({'child_reactions':SubReactions})
    #BabyInfo.update({'LevelInfo':MAKE_NODES(REACTIONSOURCE.get(ReactName),NewCenter,InAngle,PrimaryMetabolites,U)})
    Kale=MAKE_NODES(MaxAndMin,REACTIONSOURCE.get(ReactName),NewCenter,InAngle,PrimaryMetabolites,U,Metabolites,ReactMet)
    BabyInfo.update({ReactName:Kale[0]})
    MaxAndMin=Kale[1]
    for Lev in SubReactions:
        BabyInfo.update(Lev)
    return([BabyInfo,MaxAndMin])


def UPDATE_COORD(Child,ChildOffset):
    OldCent=Child['center']
    Child.update({'center':(OldCent[0]+ChildOffset[0],OldCent[1]+ChildOffset[1])})
    MM=Child['maxandmins']
    Child.update({'maxandmins':(MM[0]+ChildOffset[0],MM[1]+ChildOffset[1],MM[2]+ChildOffset[0],MM[3]+ChildOffset[1])})
    for n in Child['nodes']:
        Oldx=Child['nodes'][n]['x']
        Child['nodes'][n].update({'x':Oldx+ChildOffset[0]})
        Oldy=Child['nodes'][n]['y']
        Child['nodes'][n].update({'y':Oldy+ChildOffset[1]})
        if 'label_x' in Child['nodes'][n]:
            OldLx=Child['nodes'][n]['label_x']
            Child['nodes'][n].update({'label_x':OldLx+ChildOffset[0]})
            OldLy=Child['nodes'][n]['label_y']
            Child['nodes'][n].update({'label_y':OldLy+ChildOffset[1]})
    for R in Child['reactions']:
        OldLy=Child['reactions'][R]['label_y']
        Child['reactions'][R].update({'label_y':OldLy+ChildOffset[1]})
        OldLx=Child['reactions'][R]['label_x']
        Child['reactions'][R].update({'label_x':OldLx+ChildOffset[0]})
        for Seg in Child['reactions'][R]['segments']:
            if 'x' in Child['reactions'][R]['segments'][Seg]['b1']:
                oldx=Child['reactions'][R]['segments'][Seg]['b1']['x']
                Child['reactions'][R]['segments'][Seg]['b1'].update({'x':oldx+ChildOffset[0]})
                oldy=Child['reactions'][R]['segments'][Seg]['b1']['y']
                Child['reactions'][R]['segments'][Seg]['b1'].update({'y':oldy+ChildOffset[1]})
                oldx=Child['reactions'][R]['segments'][Seg]['b2']['x']
                Child['reactions'][R]['segments'][Seg]['b2'].update({'x':oldx+ChildOffset[0]})
                oldy=Child['reactions'][R]['segments'][Seg]['b2']['y']
                Child['reactions'][R]['segments'][Seg]['b2'].update({'y':oldy+ChildOffset[1]})

def UPDATE_MAX_AND_MIN(Child):
    MAXANDMIN=(Child['center'][0],Child['center'][1],Child['center'][0],Child['center'][1])
    for n in Child['nodes']:
        MAXANDMIN=MMAAXX(MAXANDMIN,(Child['nodes'][n]['x'],Child['nodes'][n]['y']))
        if 'label_x' in Child['nodes'][n]:
            MAXANDMIN=MMAAXX(MAXANDMIN,(Child['nodes'][n]['label_x'],Child['nodes'][n]['label_y']))
    Child.update({'maxandmins':MAXANDMIN})  

def PLACE_CORE_REACTION(Relations,Cou,REACTIONSOURCE,Metabolites): #returns dic with react name : node and connection information
    import math
    MaxAndMin=(0,0,0,0)
    SubMapInfo=dict()
    #Reactions are points, Connecting metabolites are points with names
    #will work to convert to Escher map things on upper level
    ReactName=Cou.most_common(1)[0][0]
    AdjCount=Cou.most_common(1)[0][1]
    CENTER=(0,0)
    SubMapInfo.update({'core_name':ReactName})
    SubMapInfo.update({'center':CENTER})
    U=300 #How far each unit is
    SubMapInfo.update({'Unit_Length':U})
    del Cou[ReactName]
    ReaRel=Relations.pop(ReactName)
    PrimaryMetabolites=dict()
    SubReactions=list()
    for NUM in range(AdjCount):
        Angle=math.radians(NUM*360/(AdjCount)+90)
        Consider=ReaRel.pop()
        MetCoord=(U*math.cos(Angle),U*math.sin(Angle))
        MaxAndMin=MMAAXX(MaxAndMin,MetCoord)
        PrimaryMetabolites.update({Consider[1]:MetCoord})
        if Consider[0] in Relations:
            NewRCent=(2*U*math.cos(Angle),2*U*math.sin(Angle))
            MaxAndMin=MMAAXX(MaxAndMin,NewRCent)
            Loona=PLACE_SUB_REACTION(MaxAndMin,ReactName,U,Relations,Cou,Angle,NewRCent,Consider,MetCoord,REACTIONSOURCE,Metabolites)
            SubReactions.append(Loona[0])
            MaxAndMin=Loona[1]
    #SubMapInfo.update({'primay_metabolites':PrimaryMetabolites})
    #SubMapInfo.update({'child_reactions':SubReactions})
    Kale=MAKE_NODES(MaxAndMin,REACTIONSOURCE.get(ReactName),CENTER,math.radians(90),PrimaryMetabolites,U,Metabolites)
    BabyInfo=({ReactName:Kale[0]})
    MaxAndMin=Kale[1]
    for Lev in SubReactions:
        BabyInfo.update(Lev)
    SubMapInfo.update({'reactions':BabyInfo})
    SubMapInfo.update({'maxandmins':MaxAndMin})
    return(SubMapInfo)

def GEN_CLUSTERS(GoodMet,ReactGroup,RIndexedByID,Metabolites):
        from collections import Counter
        #Format as follows:
        #Dict of reactions
        #data is tuples, (adj-reaction , metabolite link)
        ReactionLink=dict()
        ReactionAdjCount=Counter()
        for MLite in GoodMet:
            REARS=ReactGroup.get(MLite).split(',')
            ReactionAdjCount.update([REARS[0]])
            ReactionAdjCount.update([REARS[1]])
            if REARS[0] in ReactionLink:
                KillerApp=list(ReactionLink.get(REARS[0]))
                KillerApp.append((REARS[1],MLite))
                ReactionLink.update({REARS[0]:KillerApp})
            else:
                ReactionLink.update({REARS[0]:[(REARS[1],MLite)]})
            if REARS[1] in ReactionLink:
                KillerApp=list(ReactionLink.get(REARS[1]))
                KillerApp.append((REARS[0],MLite))
                ReactionLink.update({REARS[1]:KillerApp})
            else:
                ReactionLink.update({REARS[1]:[(REARS[0],MLite)]})
        MapInfo=dict()
        G=0
        while len(ReactionLink)>0:
            MapInfo.update({G:PLACE_CORE_REACTION(ReactionLink,ReactionAdjCount,RIndexedByID,Metabolites)})
            G=G+1
        return(MapInfo)


def COMBINE_TWO_METABOLITES(Clusters,NeverCombine): #New cluster format with nodes on top level  #Modified to also combine 3
    from collections import Counter
    for C in Clusters:
        MetCount=Counter()
        NODES=Clusters[C]['nodes']
        for n in list(NODES):
            NODE=NODES[n]
            if NODE['node_type']=='metabolite':
                MetCount.update([NODE.get('bigg_id')])
        CombList=set()
        for i in MetCount.items():
            if (i[1]==2 or i[1]==3) and i[0] not in NeverCombine:
                CombList.add(i[0])
        RepList=dict()
        HaveMet=dict()
        for n in list(NODES):
            NODE=NODES[n]
            if NODE['node_type']=='metabolite':
                if NODE.get('bigg_id') in HaveMet and NODE.get('bigg_id') in CombList:
                    RepList.update({n:HaveMet.get(NODE.get('bigg_id'))})
                    NODES.pop(n)
                    NODES[HaveMet.get(NODE.get('bigg_id'))].update({'node_is_primary':'true'})
                else:
                    HaveMet.update({NODE.get('bigg_id'):n})
        for R in Clusters[C]['reactions']:
            for S in Clusters[C]['reactions'][R]['segments']:
                Segment=Clusters[C]['reactions'][R]['segments'][S]
                OLDTO=Segment.get('to_node_id')
                OLDFROM=Segment.get('from_node_id')
                if OLDTO in RepList:
                    Segment.update({'to_node_id':RepList[OLDTO]})
                    #print(NODES[RepList[OLDTO]])

                    Segment['b2'].update({'x':(NODES[RepList[OLDTO]]['x']+50)})
                    Segment['b2'].update({'y':(NODES[RepList[OLDTO]]['y']+50)})
                    #Segment['b1'].update({'x':(NODES[RepList[OLDTO]]['x']+500000)})
                    #Segment['b1'].update({'y':(NODES[RepList[OLDTO]]['y']+500000)})

                if OLDFROM in RepList:
                    Segment.update({'from_node_id':RepList[OLDFROM]})
                    print('HAPPENED')
                    #Note, this seems to be something that does not happen, I guess we always go from
                    #midmarkers to nodes, never from nodes to midmarkers


def UPDATE_NODE_ID(Clusters,GloNodeID,GloSegID,NeverCombine):
    #Updates the node ID 
    from collections import Counter
    for C in Clusters:
        CLUS=Clusters.get(C)
        MetCount=Counter()
        for R in CLUS.get('reactions'):
            REACT=CLUS.get('reactions').get(R)
            TempNodeStore=dict()
            RepList=dict()
            for n in list(REACT.get('nodes')):
                NODE=REACT.get('nodes').pop(n)
                if NODE.get('node_type')=='metabolite':
                    MetCount.update([NODE.get('bigg_id')])
                TempNodeStore.update({GloNodeID:NODE})
                RepList.update({n:GloNodeID})
                GloNodeID+=1
            TempSegStore=dict()
            for s in list(REACT.get('segments')):
                SEGMENT=REACT.get('segments').pop(s)
                OLDTO=SEGMENT.get('to_node_id')
                OLDFROM=SEGMENT.get('from_node_id')
                SEGMENT.update({'to_node_id':RepList.get(OLDTO)})
                SEGMENT.update({'from_node_id':RepList.get(OLDFROM)})
                TempSegStore.update({GloSegID:SEGMENT})
                GloSegID+=1
            REACT.update({'nodes':TempNodeStore})
            REACT.update({'segments':TempSegStore})
    CombList=set()
    if False:
        for i in MetCount.items():
            if i[1]==2 and i[0] not in NeverCombine:
                CombList.add(i[0])
        RepList=dict()
        HaveMet=dict()
        for R in CLUS.get('reactions'):
            REACT=CLUS.get('reactions').get(R)
            for n in list(REACT.get('nodes')):
                NODE=REACT.get('nodes').get(n)
                if NODE.get('node_type')=='metabolite':
                    if NODE.get('bigg_id') in HaveMet and NODE.get('bigg_id') in CombList:
                        RepList.update({n:HaveMet.get(NODE.get('bigg_id'))})
                        REACT['nodes'][HaveMet.get(NODE.get('bigg_id'))].update({'node_is_primary':'true'})
                        del(REACT.get('nodes')[n])
                    else:
                        HaveMet.update({NODE.get('bigg_id'):n})
        for R in CLUS.get('reactions'):
            REACT=CLUS.get('reactions').get(R)
            for s in list(REACT.get('segments')):
                SEGMENT=REACT.get('segments').get(s)
                OLDTO=SEGMENT.get('to_node_id')
                OLDFROM=SEGMENT.get('from_node_id')
                if OLDTO in RepList:
                    SEGMENT.update({'to_node_id':RepList.get(OLDTO)})
                if OLDFROM in RepList:
                    SEGMENT.update({'from_node_id':RepList.get(OLDFROM)})
    return([GloNodeID,GloSegID])




def COMBINE_NEW_NODE(MetID,MC,F,Child):
    for n in list(Child['nodes']):
        NODE=Child['nodes'][n]
        if NODE['node_type']=='metabolite':
            if NODE['bigg_id']==MetID:
                NODE.update({'node_is_primary':'true'})
                ollx=NODE['label_x']
                olly=NODE['label_y']
                ox=NODE['x']
                oy=NODE['y']
                OFFX=ollx-ox
                OFFY=olly-oy
                NODE.update({'label_x':OFFX+MC[0]})
                NODE.update({'label_y':OFFY+MC[1]})
                NODE.update({'x':MC[0]})
                NODE.update({'y':MC[1]})
                RepID=n
    #Need to update the segments of F with the new b1 and b2
    #I believe b1 is associated with the from
    #and b2 the to
    #we will set the b1 and b2 (when appropriate) to be MC+1
    for R in Child['reactions']:
        for S in Child['reactions'][R]['segments']:
            if Child['reactions'][R]['segments'][S]['from_node_id']== RepID:
                Child['reactions'][R]['segments'][S]['b1'].update({'x':MC[0]+1})
                Child['reactions'][R]['segments'][S]['b1'].update({'y':MC[1]+1})
            if Child['reactions'][R]['segments'][S]['to_node_id']== RepID:
                Child['reactions'][R]['segments'][S]['b2'].update({'x':MC[0]+1})
                Child['reactions'][R]['segments'][S]['b2'].update({'y':MC[1]+1})
                
    for n in list(F['nodes']):
        NODE=F['nodes'][n]
        if NODE['node_type']=='metabolite':
            if NODE['bigg_id']==MetID:
                RepL=n
                F['nodes'].pop(n)
                
    for R in F['reactions']:
        for S in F['reactions'][R]['segments']:
            if F['reactions'][R]['segments'][S]['from_node_id']== RepL:
                F['reactions'][R]['segments'][S].update({'from_node_id':RepID})
                F['reactions'][R]['segments'][S]['b1'].update({'x':MC[0]+1})
                F['reactions'][R]['segments'][S]['b1'].update({'y':MC[1]+1})
                
            if F['reactions'][R]['segments'][S]['to_node_id']== RepL:
                F['reactions'][R]['segments'][S].update({'to_node_id':RepID})
                F['reactions'][R]['segments'][S]['b2'].update({'x':MC[0]+50})
                F['reactions'][R]['segments'][S]['b2'].update({'y':MC[1]+50})


def FIND_ANGLE(F,MetName):
    import math
    for NONUM in F['nodes']:
        if F['nodes'][NONUM]['node_type']=='metabolite':
            if F['nodes'][NONUM]['bigg_id']==MetName:
                MetX=F['nodes'][NONUM]['x']
                MetY=F['nodes'][NONUM]['y']
                break
    RelPoints=(MetX-F['center'][0],MetY-F['center'][1])
    Dist=math.hypot(F['center'][0] - MetX, F['center'][1] - MetY)
    if RelPoints[0]==0:
        if RelPoints[1]==0:
            Angle=0
        else:
            Angle=math.asin(RelPoints[1]/Dist)
    else:
        Angle=math.acos(RelPoints[0]/Dist)
    return(Angle)

def ROTATE_POINT(Point,angle,center):
    import math
    s=math.sin(angle)
    c=math.cos(angle)
    PX=Point[0]-center[0]
    PY=Point[1]-center[1]
    return((PX*c-PY*s,PX*s+PY*c))


def ROTATE_ABOUT_CENTER(Child,RotateBy):
    RotateCenter=Child['center']
    MM=Child['maxandmins']
    #We have a rectangle (MaxX,MaxY), (MaxX,MinY), (MinX,MaxY), (MixX,MinY)
    Cornerone=ROTATE_POINT((MM[0],MM[1]),RotateBy,RotateCenter)
    Cornertwo=ROTATE_POINT((MM[0],MM[3]),RotateBy,RotateCenter)
    Cornerthree=ROTATE_POINT((MM[2],MM[1]),RotateBy,RotateCenter)
    Cornerfour=ROTATE_POINT((MM[1],MM[3]),RotateBy,RotateCenter)
    MM=MMAAXX(MM,Cornerone)
    MM=MMAAXX(MM,Cornertwo)
    MM=MMAAXX(MM,Cornerthree)
    MM=MMAAXX(MM,Cornerfour)
    #Maxes=ROTATE_POINT((MM[0],MM[1]),RotateBy,RotateCenter)
    #Mins=ROTATE_POINT((MM[2],MM[3]),RotateBy,RotateCenter)
    Child.update({'maxandmins':MM})
    for n in Child['nodes']:
        Oldx=Child['nodes'][n]['x']
        Oldy=Child['nodes'][n]['y']
        News=ROTATE_POINT((Oldx,Oldy),RotateBy,RotateCenter)
        Child['nodes'][n].update({'x':News[0]})
        Child['nodes'][n].update({'y':News[1]})
        if 'label_x' in Child['nodes'][n]:
            OldLx=Child['nodes'][n]['label_x']
            OldLy=Child['nodes'][n]['label_y']
            News=ROTATE_POINT((OldLx,OldLy),RotateBy,RotateCenter)
            Child['nodes'][n].update({'label_x':News[0]})
            Child['nodes'][n].update({'label_y':News[1]})
    for R in Child['reactions']:
        OldLx=Child['reactions'][R]['label_x']
        OldLy=Child['reactions'][R]['label_y']
        News=ROTATE_POINT((OldLx,OldLy),RotateBy,RotateCenter)
        Child['reactions'][R].update({'label_y':News[1]})
        Child['reactions'][R].update({'label_x':News[0]})
        for Seg in Child['reactions'][R]['segments']:
            if 'x' in Child['reactions'][R]['segments'][Seg]['b1']:
                oldx=Child['reactions'][R]['segments'][Seg]['b1']['x']
                oldy=Child['reactions'][R]['segments'][Seg]['b1']['y']
                News=ROTATE_POINT((oldx,oldy),RotateBy,RotateCenter)
                Child['reactions'][R]['segments'][Seg]['b1'].update({'x':News[0]})
                Child['reactions'][R]['segments'][Seg]['b1'].update({'y':News[1]})
                oldx=Child['reactions'][R]['segments'][Seg]['b2']['x']
                oldy=Child['reactions'][R]['segments'][Seg]['b2']['y']
                News=ROTATE_POINT((oldx,oldy),RotateBy,RotateCenter)             
                Child['reactions'][R]['segments'][Seg]['b2'].update({'x':News[0]})
                Child['reactions'][R]['segments'][Seg]['b2'].update({'y':News[1]})    



def CORRECT_CHILD(Child,MetName,Angle):
    #Given a cluster called Child, and the name of a metabolite called MetName, and the incoming angle of connection
    #rotate child about it's center so that the node for MetName lies on the ray from center through -Angle
    OldAngle=FIND_ANGLE(Child,MetName)
    RotateBy=OldAngle+Angle
    ROTATE_ABOUT_CENTER(Child,RotateBy)



def FIND_RECTANGLE_INTERSECTION(M,Angle):
    import math
    OffX=(M[0]+M[2])/2
    OffY=(M[1]+M[3])/2
    AdjM=(M[0]-OffX,M[1]-OffY,M[2]-OffX,M[3]-OffY)
    a=AdjM[0]-AdjM[2]
    b=AdjM[1]-AdjM[3]

    R1=math.atan2(a/2,b/2)
    R2=math.atan2(-a/2,b/2)
    R3=math.atan2(-a/2,-b/2)
    R4=math.atan2(a/2,-b/2)

    if Angle%(2*math.pi)<=R1 or Angle%(2*math.pi)>=R4 :
        X=a/2
        Y=math.sin(Angle)
    elif (Angle%(2*math.pi)>=R2 and Angle%(2*math.pi)<=R3):
        X=-a/2
        Y=math.sin(Angle)
    elif (Angle%(2*math.pi)>R1 and Angle%(2*math.pi)<R2):
        X=math.cos(Angle)
        Y=b/2
    else:
        X=math.cos(Angle)
        Y=-b/2
    
    PointOfIntersection=(X+OffX,Y+OffY)
    return PointOfIntersection
    
    
def SUB_COMBINE(ChainLink,Clusters,Angle,ClusterCount,ChainLinkID):
    import math
    if ChainLinkID in ChainLink:
        if len(ChainLink[ChainLinkID])>0:
            ConnList=ChainLink.pop(ChainLinkID)
            F=Clusters[ChainLinkID]
            UPDATE_MAX_AND_MIN(F)
            NumDirChild=len(ConnList)
            for BUM in range(1,NumDirChild+1):
                NUM=BUM-1
                MetName=ConnList[NUM][1]
                #NewAngle=FIND_ANGLE(F,MetName)
                #NewAngle=NewAngle+math.radians(NUM*90/(NumDirChild))
                NewAngle=(Angle-math.pi)+math.radians(BUM*90/(NumDirChild+1))%(2*math.pi)

                RecF=FIND_RECTANGLE_INTERSECTION(F['maxandmins'],NewAngle)
                
                CID=ConnList[NUM][0]
                Child=Clusters[CID]
                CORRECT_CHILD(Child,MetName,NewAngle)
                if (ChainLinkID,ConnList[NUM][1]) in ChainLink[CID]:
                        ChainLink[CID].remove((ChainLinkID,ConnList[NUM][1]))
                UPDATE_MAX_AND_MIN(Child)
                RecC=FIND_RECTANGLE_INTERSECTION(Child['maxandmins'],(NewAngle+math.pi)%(2*math.pi))

                DistF=math.hypot(F['center'][0]-RecF[0],F['center'][1]-RecF[1])
                DistC=math.hypot(Child['center'][0]-RecF[0],Child['center'][1]-RecF[1])

                Dist=DistF+DistC+50

                ChildOffset=(F['center'][0]+Dist*math.cos(NewAngle),F['center'][1]+Dist*math.sin(NewAngle))
                UPDATE_COORD(Child,ChildOffset)
                MetID=ConnList[NUM][1]
                MetNodeCoord=(F['center'][0]+(Dist/2)*math.cos(NewAngle),F['center'][1]+(Dist/2)*math.sin(NewAngle))
                COMBINE_NEW_NODE(MetID,MetNodeCoord,F,Child)                
                #from newly placed child, do any lower level combinations
                SUB_COMBINE(ChainLink,Clusters,NewAngle,ClusterCount,CID)  
                
                #Now combine the child and the father
                C=Clusters.pop(CID)
                if CID in ChainLink:
                    Dump=ChainLink.pop(CID)
                #print(str(ConnList[NUM][0])+' was popped')
                
                F['nodes'].update(C['nodes'])
                F['reactions'].update(C['reactions'])

                #UPDATE_MAX_AND_MIN(F)
                #MaxAndMin=F['maxandmins']
                #F.update({'center':(((MaxAndMin[0]+MaxAndMin[2])/2.0),((MaxAndMin[1]+MaxAndMin[3])/2.0))})

                #UPDATE_COORD(F,(-((MaxAndMin[0]+MaxAndMin[2])/2.0),-((MaxAndMin[1]+MaxAndMin[3])/2.0)))#center each cluster at (0,0)
                #F.update({'center':(0,0)})

            
def COMBINE_CLUSTERS(SecMetLinkList,Clusters):
    from collections import Counter
    import math
    ClusterCount=Counter()
    ChainLink=dict()
    for i in SecMetLinkList:
        P=SecMetLinkList.get(i).split(',')
        REARS=[int(P[0]),int(P[1])]
        if REARS[0]!=REARS[1]:
            if REARS[0] in ChainLink:
                KillerApp=list(ChainLink.get(REARS[0]))
                Bad=1
                for Super in KillerApp:
                    if REARS[1] == Super[0]:
                        Bad=0
                if Bad==1:
                    KillerApp.append((REARS[1],i))
                    ClusterCount.update([REARS[0]])
                ChainLink.update({REARS[0]:KillerApp})
            else:
                ChainLink.update({REARS[0]:[(REARS[1],i)]})
                ClusterCount.update([REARS[0]])
            if REARS[1] in ChainLink:
                KillerApp=list(ChainLink.get(REARS[1]))
                Bad=1
                for Super in KillerApp:
                    if REARS[0]==Super[0]:
                        Bad=0
                if Bad==1:
                    KillerApp.append((REARS[0],i))
                    ClusterCount.update([REARS[1]])
                ChainLink.update({REARS[1]:KillerApp})
            else:
                ChainLink.update({REARS[1]:[(REARS[0],i)]})
                ClusterCount.update([REARS[1]])
    IterateChainLink=list(ChainLink)
    while len(IterateChainLink)>0:
        IterateChainLink=list(ChainLink) 
        for i in IterateChainLink:       
            if i in ChainLink:                             
                ConnList=ChainLink.pop(i)
                NumDirChild=ClusterCount[i]   
                F=Clusters[i]               
                for NUM in range(NumDirChild):
                    #Figure out the metabolite used in the combo
                    #get the metabolite's coordinates
                    #set angle to be on the line from center through metabolite
                    #Find met's coordinates in child relative to child center
                    #rotate child so met is closest to father
                    MetName=ConnList[NUM][1]
                    Angle=FIND_ANGLE(F,MetName)
                    #Now because it may be the case that a bunch of metabolites are off at the same ange
                    #we add an offset to the ange depending on the current NUM
                    if NumDirChild>1:
                        Angle=Angle+math.radians(NUM*90/(NumDirChild))
                    RecF=FIND_RECTANGLE_INTERSECTION(F['maxandmins'],Angle)
                    CID=ConnList[NUM][0]
                    Child=Clusters[ConnList[NUM][0]]
                    CORRECT_CHILD(Child,MetName,Angle) #Rotates the child so the metabolite is closest to the inangle
                    if (i,ConnList[NUM][1]) in ChainLink[CID]:
                        ChainLink[CID].remove((i,ConnList[NUM][1]))
                    UPDATE_MAX_AND_MIN(Child)
                    RecC=FIND_RECTANGLE_INTERSECTION(Child['maxandmins'],(Angle+math.pi)%(2*math.pi))
                    DistF=math.hypot(F['center'][0]-RecF[0],F['center'][1]-RecF[1])
                    DistC=math.hypot(Child['center'][0]-RecF[0],Child['center'][1]-RecF[1])
                    Dist=DistF+DistC+50

                    ChildOffset=(F['center'][0]+Dist*math.cos(Angle),F['center'][1]+Dist*math.sin(Angle))
                    UPDATE_COORD(Child,ChildOffset)
                    MetID=ConnList[NUM][1]
                    MetNodeCoord=(F['center'][0]+(Dist/2)*math.cos(Angle),F['center'][1]+(Dist/2)*math.sin(Angle))
                    COMBINE_NEW_NODE(MetID,MetNodeCoord,F,Child)                
                    #from newly placed child, do any lower level combinations
                    SUB_COMBINE(ChainLink,Clusters,Angle,ClusterCount,ConnList[NUM][0])

                    #Now combine the child and the father
                    C=Clusters.pop(CID)
                    if CID in ChainLink:
                        Dump=ChainLink.pop(CID)
                        #print('POPPED '+str(CID))
                    #print(str(ChainLink[i][NUM][0])+' was popped')
                    F['nodes'].update(C['nodes'])
                    F['reactions'].update(C['reactions'])
                    #UPDATE_MAX_AND_MIN(F)
                    #MaxAndMin=F['maxandmins']
                    #UPDATE_COORD(F,(-((MaxAndMin[0]+MaxAndMin[2])/2.0),-((MaxAndMin[1]+MaxAndMin[3])/2.0)))#center each cluster at (0,0)
                    #F.update({'center':(0,0)})
                UPDATE_MAX_AND_MIN(F)
                MaxAndMin=F['maxandmins']
                UPDATE_COORD(F,(-((MaxAndMin[0]+MaxAndMin[2])/2.0),-((MaxAndMin[1]+MaxAndMin[3])/2.0)))#center each cluster at (0,0)
                F.update({'center':(0,0)})

            
                
def COUNT_METABOLITES_IN_CLUSTERS(Clusters):
    from collections import Counter
    NewMetCount=Counter()
    for i in Clusters:
        TinyCount=Counter()
        for No in Clusters[i]['nodes']:
            if Clusters[i]['nodes'][No].get('node_type')=='metabolite':
                TinyCount.update([Clusters[i]['nodes'][No].get('bigg_id')])
        for Mets in TinyCount.items():
            if Mets[1]==1:
                NewMetCount.update([Mets[0]])
    return(NewMetCount)


def GET_TWO_LINK(REACT,ERMN):        
    from collections import Counter
    MetCount=Counter()
    MetLink=dict()    
    for i in REACT:
        MetabolList=list(i.get('metabolites'))
        if len(MetabolList) <= ERMN:
            for MET in MetabolList:
                MetCount.update([MET])
                if MET in MetLink:
                    MetLink.update({MET : MetLink.get(MET)+','+i.get('id')})
                else:
                    MetLink.update({MET : i.get('id')})
    GoodMet=list()
    ReactGroup=dict()
    for i in MetCount.items():
        if i[1] == 2:
            GoodMet.append(i[0])
            ReactGroup.update({i[0]:MetLink.get(i[0])})
    return([GoodMet,ReactGroup])

def COMBINE_ALL_CLUSTERS(Clusters,NeverCombine,IncludeReactions=0):
    #Include Reactions:
    #0 means do this on the entire cluster including clusters of single reaction only
    #1 means just work on the clusters of more than 1 reaction
    #2 means just work on the clusters of 1 reaction only
    from collections import Counter
    ExcludedClusters=set()
    if IncludeReactions==1:
        for i in Clusters:
            if len(Clusters[i]['reactions'])==1:
                ExcludedClusters.add(i)
    elif IncludeReactions==2:
        for i in Clusters:
            if len(Clusters[i]['reactions'])>1:
                ExcludedClusters.add(i)
    tempclus=dict()
    for Clus in Clusters:
        if Clus not in ExcludedClusters:
            tempclus.update({Clus:Clusters[Clus]})
    #print(len(Clusters))
    #print(len(tempclus))
    NewMetCount=COUNT_METABOLITES_IN_CLUSTERS(tempclus)
    ActiveList=list()
    for i in NewMetCount.items():
        if i[1]==2:
            ActiveList.append(i[0])
    while len(ActiveList)>0:        
        SecMetLinkList=dict()
        for Met in ActiveList:
            SecMetLinkList.update({Met:''})                
            for i in tempclus:
                TinyCount=Counter()
                for No in tempclus[i]['nodes']:
                    if tempclus[i]['nodes'][No].get('node_type')=='metabolite':
                        if tempclus[i]['nodes'][No].get('bigg_id') == Met:
                            TinyCount.update([Met])
                if TinyCount[Met]==1:
                    Old=SecMetLinkList.get(Met)
                    if Old =='':
                        New=str(i)
                    else:
                        New=Old+','+str(i)
                    SecMetLinkList.update({Met:New})
        COMBINE_CLUSTERS(SecMetLinkList,Clusters)
        COMBINE_TWO_METABOLITES(Clusters,NeverCombine)
        tempclus=dict()
        for Clus in Clusters:
            if Clus not in ExcludedClusters:
                tempclus.update({Clus:Clusters[Clus]})
        NewMetCount=COUNT_METABOLITES_IN_CLUSTERS(tempclus)
        ActiveList=list()
        for i in NewMetCount.items():
            if i[1]==2:
                ActiveList.append(i[0])
                
def AUTO_MAP_LINKED_TWO(Cobra_Model,With_Names=False):
    import json
    if With_Names==True:
        print('Using metabolite Names in map instead of BiGG ID')
        import os.path
        if os.path.isfile(Cobra_Model.replace('.json','_BiGG_NAME_Metabolites.json'))==True:
            print('found existing reference file '+str(Cobra_Model.replace('.json','_BiGG_NAME_Metabolites.json')))
        else:
            print('Could not find '+Cobra_Model.replace('.json','_BiGG_NAME_Metabolites.json')+'. Generating new replacement mapping from BiGG Database')
            GATHER_BIGG_DESCRIPTIVE_NAMES(Cobra_Model)
        with open(Cobra_Model.replace('.json','_BiGG_NAME_Metabolites.json'),'r') as f:
            NameKey=json.load(f)
    


    
    ERMN=15 #Max number of metabolites in a reaction before the reaction is ignored

    from collections import Counter
    
    print('Loading Cobra Model')
    with open(Cobra_Model,'r') as f:
        M=json.load(f)
    if With_Names==True:
        #for R in M['reactions']:
            #TMets=dict()
            #for i in R['metabolites']:
                #if i in NameKey:
                #    TMets.update({NameKey[i]:R['metabolites'][i]})
            #R.update({'metabolites':TMets})
        for mets in M['metabolites']:
            if mets['id'] in NameKey:
                OldID=mets['id']
                mets.update({'name':NameKey[OldID]})
    else:
        NameKey=dict()
        for i in M['metabolites']:
            NameKey.update({i['id']:i['name']})

    Metabolites=NameKey
    REACT=M.get('reactions')
    RIndexedByID=dict()
    for i in REACT:
        RIndexedByID.update({i.get('id'):i})

    NAME=M.get('description')+'_Auto_Built_Escher_Map_V2'

    H=GET_TWO_LINK(REACT,ERMN)
    GoodMet=H[0]
    ReactGroup=H[1]

    UsedReactions=set()
    Clusters=GEN_CLUSTERS(GoodMet,ReactGroup,RIndexedByID,Metabolites)
    while len(GoodMet)>0:
        for i in Clusters:
            UsedReactions=UsedReactions|set(Clusters[i]['reactions'])
        Again=list()
        for R in REACT:
            if R['id'] not in UsedReactions:
                Again.append(R)
        #print(len(GoodMet))
        H=GET_TWO_LINK(Again,ERMN)
        GoodMet=H[0]
        ReactGroup=H[1]
        ClustersTWO=GEN_CLUSTERS(GoodMet,ReactGroup,RIndexedByID,Metabolites)
        GLen=len(Clusters)
        for i in ClustersTWO:
            Clusters.update({GLen+i:ClustersTWO[i]})
    MetCount=Counter()
    VALID=dict()
    for R in REACT:
        if len(R['metabolites'])<ERMN:
            VALID.update({R['id']:R})
            for met in R['metabolites']:
                MetCount.update([met])
    NeverCombine=set()
    for i in MetCount.items():
        if i[1]>5:
            NeverCombine.add(i[0])
    #import sys
    #sys.exit(0)

    ###################################

    UnusedReaction=set(VALID)-UsedReactions
    UnusedReactionDict=dict()
    for RName in list(UnusedReaction):
        UnusedReactionDict.update({RName:VALID.get(RName)})
    
    #ExcludeList=list(set(MetCount)-(set(GoodMet)|set(ThreeFourFiveMet)))
    #now have three lists.  GoodMet contains all metabolites that lie in two reactions
    #ThreeFourFiveMet all metabolites lying in 3 4 or 5 reactions
    #ExcludeList all metabolites lying in 1 reaction or more that 5 reactions.





#These are clusters of more than 1

    #Now we want to count the metabolites in reactions which did not get clustered and use those to make more clusters (as a loop)


    #now adding the unused good reactions as clusters of 1
    EmptyLink=dict()
    TrivialCounter=Counter()
    for i in UnusedReactionDict:
        EmptyLink.update({i:''})
        TrivialCounter[i]=0
    G=len(Clusters)
    while len(EmptyLink)>0:
        Clusters.update({G:PLACE_CORE_REACTION(EmptyLink,TrivialCounter,RIndexedByID,Metabolites)})
        G+=1
    
 
    #Now have a dict of clusters
    #for each index 0,1,...
    #there is center, core_name, maxandmins, reactions, and unit_length

    #in 'reactions' is the names of each reaction
    #in each reaction are Nodes, x,y labels ('label_x) and segments

    #First, go through all clusters and update node id with a global node id
    GloNodeID=10000
    GloSegID=100000
    ID=UPDATE_NODE_ID(Clusters,GloNodeID,GloSegID,NeverCombine)#updates the node ids with global ids, also combines metabolites inside
    GloNodeID=ID[0]
    GloSegID=ID[1] #Now have Clusters, each cluster has center,core_name,maxandmins,reactions,unit_length
    #want to make it so we can do a next level of combination
    #Should make it so that segments linked with reaction. Nodes up on top
    
    for i in Clusters:
        Nodes=dict()
        for React in Clusters[i]['reactions']:
            Nodes.update(Clusters[i]['reactions'][React].pop('nodes'))
        Clusters[i].update({'nodes':Nodes})
        MaxAndMin=Clusters[i]['maxandmins']
        UPDATE_COORD(Clusters[i],(-((MaxAndMin[0]+MaxAndMin[2])/2.0),-((MaxAndMin[1]+MaxAndMin[3])/2.0)))#center each cluster at (0,0)
        Clusters[i].update({'center':(0,0)})

    
    #clusters[clusterID]['nodes'] contains all the node info for that cluster
    #now we want to query metabolite nodes
    #we want a list of metabolites which appear exactly twice, appear means:
        #either metabolite is on exactly one node in a cluster
        #or metabolite is present in an unused unexcluded reaction
   
    #count instances of metabolites in clusters
    Start=1
    End=0
    while Start-End!=0:
        Start=len(Clusters)
        print(Start)
        COMBINE_ALL_CLUSTERS(Clusters,NeverCombine,IncludeReactions=0)
        COMBINE_ALL_CLUSTERS(Clusters,NeverCombine,IncludeReactions=1)
        COMBINE_ALL_CLUSTERS(Clusters,NeverCombine,IncludeReactions=2)
        End=len(Clusters)
        print(End)

  

    #Now all the reactions (with fewer than ERMN metabolites) are on the board. All the metabolites
    #which are in either exactly two reactions or are in exactly two clusters have been combined.

    #Now print them, make the Escher map
    LBoundX=0
    UpBoundY=0
    MaxY=0
    MaxX=0
    #(0,0) is upper left corner,NOTE that each cluster is centered exactly at 0,0.  center of box
    #we have a function called UPDATE_COORD which we send the offset.
    for i in Clusters:
        if len(Clusters[i]['reactions'])!=1:
            Current=Clusters.get(i)
            MaxAndMin=Current['maxandmins']
            OFFSET=(LBoundX-MaxAndMin[2]+400,UpBoundY-MaxAndMin[3]+100)
            UPDATE_COORD(Current,OFFSET)
            MaxAndMin=Current['maxandmins']
            LBoundX=MaxAndMin[0]
            MaxY=max(MaxY,MaxAndMin[1])
            MaxX=max(MaxX,LBoundX+100)
            if LBoundX>20000:
                LBoundX=0
                UpBoundY=MaxY+100
    #Now prepare Clusters for printing.
    HEADER={'map_name':NAME,'map_id':'MARKPROGID','map_description':'Map was build using a python collection script written by Mark Layer on 7/10/2015','homepage':"https://escher.github.io",'schema':'https://escher.github.io/escher/jsonschema/1-0-0#'}
    #SegID and NodeID already universal
    #Need to assign ReactID
    #Need to collect Reactions, and Nodes
    
    ReactID=0
    Nodes=dict()
    Reactions=dict()
    Genes=dict()

    for GEN in M.get('genes'):
        Genes.update({GEN.get('id'):GEN})

    for i in Clusters:
        if len(Clusters[i]['reactions'])!=1:
            for ReactName in Clusters.get(i).get('reactions'):
                CR=Clusters.get(i).get('reactions').get(ReactName)
                CobraInfo=RIndexedByID.get(ReactName)
                genes=list()
                met=list()
                TempSegs=dict()
                for gg in ','.join(','.join(CobraInfo.get('gene_reaction_rule').split(' or ')).split(' and ')).split(','):
                    WORD=str(Genes.get(gg.replace('(','').replace(')',''))).replace('id','bigg_id')
                    genes.append(WORD)
                for mm in CobraInfo.get('metabolites'):
                    met.append({'bigg_id':mm,'coefficient':CobraInfo.get('metabolites').get(mm)})     
                if CobraInfo.get('lower_bound')<0:
                    Rev='true'
                else:
                    Rev='false'
                Reactions.update({ReactID:{'name':CobraInfo.get('name'),'bigg_id':CobraInfo.get('id'),'reversibility':Rev,'label_x':CR.get('label_x'),'label_y':CR.get('label_y'),'gene_reaction_rule':CobraInfo.get('gene_reaction_rule'),'genes':genes,'metabolites':met,'segments':CR['segments']}})
                ReactID+=1
            Nodes.update(Clusters[i]['nodes'])

    Canvas={'x':0,'y':0,'width':MaxX,'height':MaxY}
    print(str(MaxX)+','+str(MaxY))
    REST={'reactions':Reactions,'text_labels':{},'canvas':Canvas,'nodes':Nodes}
    PrintString=str('['+json.dumps(HEADER)+','+json.dumps(REST)+']')
    PrintString=PrintString.replace('"false"','false').replace('"true"','true').replace('"null"','null')
    #PrintString='['+json.dumps(TOP)+']'
    open(NAME+'.json','w').write(PrintString)
    print('Done')                
            


def FIX_NAMES_COBRA_MODEL(Cobra_Model):
    import json
  
    print('Modifying Cobra Model file so names are descriptive names')
    import os.path
    if os.path.isfile(Cobra_Model.replace('.json','_BiGG_NAME_Metabolites.json'))==True:
        print('found existing reference file '+str(Cobra_Model.replace('.json','_BiGG_NAME_Metabolites.json')))
    else:
        print('Could not find '+Cobra_Model.replace('.json','_BiGG_NAME_Metabolites.json')+'. Generating new replacement mapping from BiGG Database')
        GATHER_BIGG_DESCRIPTIVE_NAMES(Cobra_Model)
    with open(Cobra_Model.replace('.json','_BiGG_NAME_Metabolites.json'),'r') as f:
        NameKey=json.load(f)
    with open(Cobra_Model,'r') as f:
        M=json.load(f)
        for Met in M['metabolites']:
            try:
                Met.update({'name':NameKey[Met['id'][:len(Met['id'])-2]+'_c']})
            except KeyError:
                print('skipping')
                print(Met['id'])
                print(Met['id'][:len(Met['id'])-2]+'_c')
                print(Met['id'][:len(Met['id'])-2]+'_c' in NameKey)
    
    PrintString=json.dumps(M)
    open(Cobra_Model.replace('.json','_Fixed_Names.json'),'w').write(PrintString)
            

def METABOLITE_QUERY(Escher_map,Cobra_Model=''):
    #METABOLITE_QUERY takes an Escher Map file and reads all the nodes
    #it generates a list of all the metabolites present in the map and counts
    #the number of nodes per metabolite.  If a cobra model is specified it counts all instances of the metabolites in the
    #missing reactions and prints a file, first column metabolite node count and second, missing metabolite instances.


    #Escher_map should be the file address or NAME.json if it's in the same
    #directory
    import json
    from collections import Counter

    IgnoreNumber=15 #For column 3 if exists, ignores counting if in a reaction
    #containing more than this many metabolites (BioSyn and such)
    print('loading Escher Map model '+Escher_map)
    with open(Escher_map,'r') as f:
        M=json.load(f)
    OutName=str(M[0].get('map_name')).replace(' ','_')+'_Metabolite_Report.csv'
    Nodes=M[1].get('nodes')
    React=M[1].get('reactions')
    RSet=set()
    for R in React:
        RSet.add(React[R]['bigg_id'])
    MetCount=Counter()
    print('counting metabolite nodes')
    for i in Nodes:
        if Nodes.get(i).get('node_type')=='metabolite':
            MetCount.update([Nodes.get(i).get('bigg_id')])
    LineOne='BiGG_ID,Number_Of_Nodes_In_Map'
    if Cobra_Model!='':
        with open(Cobra_Model,'r') as f:
            L=json.load(f)
        MisMetCount=Counter()
        MisMetCountExclude=Counter()
        for Kapra in L['reactions']:
            if Kapra['id'] not in RSet:
                if len(Kapra['metabolites'])<IgnoreNumber:
                    for mm in Kapra['metabolites']:
                        MisMetCountExclude.update([mm])
                for mm in Kapra['metabolites']:
                    MisMetCount.update([mm])
        LineOne=LineOne+',Number_of_metabolites_not_on_map,Number_of_metabolites_not_on_map_excluding_big_reactions'
    printlist=[LineOne]
    if Cobra_Model!='':
        for met in list(MetCount):
            Shoe=MetCount.pop(met)
            try:
                Pull=MisMetCountExclude.pop(met)
            except KeyError:
                Pull=0
            try:
                Moo=MisMetCount.pop(met)
            except KeyError:
                Moo=0
            printlist.append(met+','+str(Shoe)+','+str(Moo)+','+str(Pull))
        for met in list(MisMetCount):
            Moo=MisMetCount.pop(met)
            try:
                Pull=MisMetCountExclude.pop(met)
            except KeyError:
                Pull=0
            printlist.append(met+','+str(0)+','+str(Moo)+','+str(Pull))
    else:
         for met in list(MetCount):
             Shoe=MetCount.pop(met)
             printlist.append(met+','+str(Shoe))
    print('outputting results to '+OutName)
    open(OutName,'w').write('\n'.join(printlist))
    print('Done')


#ExList=['h_c','nad_c']
#ExList=[]
#COMBINE_LIKE_NODES_WITH_EXCLUSIONS('RealAttemptAtiJN678.json',ExList)
#METABOLITE_QUERY('MarkRealAttemptatiJN678.json')
#CHECK_FOR_MISSING_METABOLITES('SuperTestFunction5.json','iJN678.json')
#AUTO_MAP_LINKED('iJN678.json')
#AUTO_MAP_LINKED_TWO('iJN678.json',True)
#METABOLITE_QUERY('iJN678_Auto_Built_Escher_Map_V2.json','iJN678_Fixed_Names.json')
#FIX_NAMES_COBRA_MODEL('iJN678.json')
METABOLITE_QUERY('iJN678_Auto_Map_Curated.json','iJN678_Fixed_Names.json')
