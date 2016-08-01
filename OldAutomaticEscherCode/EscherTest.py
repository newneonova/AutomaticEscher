#Here I will write a function to read an escher map into a series of variables






#top level has Header, rest of tree

#rest of tree contains:


def LISTELM(M):
    for i in M:
        print(i)


import json

FILE='Map1.json'
with open(FILE,'r') as f:
    M=json.load(f)

HEADER=M[0]
REACT=M[1].get('reactions')
NODES=M[1].get('nodes')
TEXT=M[1].get('text_labels')
CANVAS=M[1].get('canvas')
