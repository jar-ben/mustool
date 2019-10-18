import sys
sys.path.insert(0, "/home/xbendik/usr/lib/lib/python3.7/site-packages")
from pysat.card import *



def generateCard(ls, b):
    print("{} choose {}".format(len(ls), b))
    am1 = CardEnc.atmost(lits=ls, bound=(b - 1), encoding=8)
    return am1.clauses

def maxVar(cls):
    m = 0
    for cl in cls:
        m = max(m, max([abs(l) for l in cl]))
    return m

def generate(n):
    lits = [i for i in range(1, n+1)]
    card = generateCard(lits, int(n/2))

    m = maxVar(card)
    res = "p gcnf {} {} {}\n".format(maxVar(card), n + len(card), n)
    for l in lits:
        res += "{{{}}} {} 0\n".format(l,l)
    for cl in card:
        res += "{0} " + " ".join([str(l) for l in cl]) + " 0\n"
    return res

with open("generated.gcnf", "w") as f:
    f.write(generate(8))
