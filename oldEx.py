def negate(cl):
    return [-1 * l for l in cl]

def nvars(cls):
    vrs = []
    for cl in cls:
        for l in cl:
            if l < 0: l *= -1
            vrs += [l]
    return set(vrs)

## input: l: a list of lists of numbers
## output: a maximum number in l
def maxVar(l):
    maximum = 0
    for i in l:
        maximum = max(maximum, max([abs(j) for j in i]))
    return maximum


def tseitinCube(cube, current):
    current += 1
    res = []
    for l in cube:
        res.append([-current,l])
    res.append([current] + [-l for l in cube])
    return res, current

def tseitinClause(clause, current):
    current += 1
    res = []
    for l in clause:
        res.append([current,-l])
    res.append([-current] + [l for l in clause])
    return res, current

## input: cnf: a cnf formula, given as a list of lists of literals
## input: current: available variable index (i.e. not representing any variable yet)
## output: cnf formula logacically equivalent to the input formula, however transformed via the tseitin transformation
## that is, several tseitin activation variables are created (starting with current) and the original cnf formula is now
## true if an only if the main tseitin activation literal holds. All the added activation literals are returned in the list "acts",
## the main activation literal is the last one in the list. 
def tseitinOnCnf(cnf, current):
    activators = []
    acts = []
    res = []
    for cl in cnf:
        tscl, current = tseitinClause(cl, current)
        activators.append(current)
        res += tscl
    current += 1
    for a in activators:
        res.append([a,-current])
    res.append([current] + negate(activators))    
    return res, current

def tseitinOnDnf(dnf, current):
    activators = []
    acts = []
    res = []
    for cl in dnf:
        tscl, current = tseitinCube(cl, current)
        activators.append(current)
        res += tscl
    current += 1    
    res.append([-current] + activators)
    for a in activators:
        res.append([-a,current])
    return res, current

## input: cl: a XOR represented as a set of literals
## output: tseitin representation of cl
def tseitinXOR(cl, current):
    fst = cl[0]
    cnf = []
    for snd in cl[1:]:
        c, current = tseitinOnCnf([[fst, snd],[-fst, -snd]], current)
        cnf += c
        fst = current
    return cnf, current

## input: cnf: a cnf formula (not necessarily in a tseitin form
## output: negation of the input in a tseitin cnf form
def negateCNF(cnf, current):
    dnf = []
    for cl in cnf:
        dnf.append([-l for l in cl])
    return tseitinOnDnf(dnf, current)

## input: two tseintin cnf formula, cn1 and cnf2, with main activators act1 and act2
## output: tseintin OR of cn1 and cn2 with main activator current
def orTwoTseitinCnf(cnf1, cnf2, act1, act2, current):
    current += 1
    return cnf1 + cnf2 + [[-current, act1, act2], [-act1, current], [-act2, current]], current

## input: two tseintin cnf formula, cn1 and cnf2, with main activators act1 and act2
## output: tseintin AND of cn1 and cn2 with main activator current
def andTwoTseitinCnf(cnf1, cnf2, act1, act2, current):
    current += 1
    return cnf1 + cnf2 + [[current, -act1, -act2], [act1, -current], [act2, -current]], current

## input: two sets of variables A = {a1, ..., an}, B = {b1, ..., bn}
## output: tseitin cnf formula enforcing A to be a subset of B (in terms of valuations)
def properSubset(A, B, current):
    assert len(A) == len(B)
    #enforce subset
    sub = []
    for i in range(len(A)):
        sub.append([-A[i],B[i]])
    #enforce proper subset
    propDNF = []
    for i in range(len(A)):
        propDNF.append([-A[i],B[i]])

    tProp, current = tseitinOnDnf(propDNF, current)
    actProp = current
    tSub, current = tseitinOnCnf(sub, current)
    return andTwoTseitinCnf(tProp, tSub, actProp, current, current)

import sys
sys.path.insert(0, "/home/xbendik/usr/lib/lib/python3.7/site-packages")
from pysat.card import *
def offByOneSubset(A, B, current):
    assert len(A) == len(B)
    sub = []
    for i in range(len(A)):
        sub.append([-A[i],B[i]])

    Z = []
    ZActs = []
    for i in range(len(A)):
        cnf, current = tseitinOnCnf([[-A[i], -B[i]], [A[i], B[i]]], current)
        Z += cnf
        ZActs.append(current)
    print("Z {} {}".format(len(Z), len(nvars(Z))))

    am1 = CardEnc.atmost(lits=ZActs, bound=1, encoding=1)
    allC = sub + Z + am1.clauses + [ZActs]
    print("allC {} {}".format(len(allC), len(nvars(allC))))
    return tseitinOnCnf(allC, current)

## intput: C: set C = [c1, ..., cn] of constraints
## input: activators: set of activation literals [x1, ..., xn]
## output: tseitin cnf formula saying that C using the activators, i.e. \bigwedge_{i}(ci or -xi), is unsatisfiable
def unsat(C, activators, current):
    assert len(C) == len(activators)
    n = len(C)
    dnf = []
    for i in range(n):
        dnf.append([activators[i]] + [-l for l in C[i]])
    return tseitinOnDnf(dnf, current)

## intput: C: set C = [c1, ..., cn] of constraints
## input: activators: set of activation literals [x1, ..., xn]
## output: tseitin cnf formula saying that C using the activators, i.e. \bigwedge_{i}(ci or -xi), is satisfiable
def sat(C, activators, current):
    assert len(C) == len(activators)
    n = len(C)
    cnf = []
    for i in range(n):
        cnf.append([-activators[i]] + [l for l in C[i]])
    return tseitinOnCnf(cnf, current)


def exMUS(C):
    Vars = list(nvars(C)) ##these are in the 2nd quantificator level and are Exists
    n = len(C)
    m = len(Vars)

    primeC = []
    for cl in C:
        primeCl = []
        for l in cl:
            if l > 0: primeCl.append(l + m)
            else: primeCl.append(l - m)
        primeC.append(primeCl)


    primeVars = list(nvars(primeC)) ##these are in 3rd quantificator level and are ForAll
    assert primeVars == [max(Vars) + i + 1 for i in range(m)]
    activators = [max(primeVars) + i + 1 for i in range(n)] ##these are in 1st quantificator level and are ForAll
    primeActivators = [max(activators) + i + 1 for i in range(n)] ##these are in the 2nd quantificator level and are Exists
    current = primeActivators[-1]

    ## the chosen subset X has to be unsat
    unsatCNF, current = unsat(C, activators, current)
    unsatAct = current

    ## there is no unsat subset of X
    propSubCNF, current = properSubset(primeActivators, activators, current)
    propSubCNF, current = offByOneSubset(primeActivators, activators, current)
    negPropSubCNF, current = negateCNF(propSubCNF, current)
    negPropSubAct = current
    satCNF, current = sat(primeC, primeActivators, current)
    satAct = current
    noUnsatSubset, current = orTwoTseitinCnf(negPropSubCNF, satCNF, negPropSubAct, satAct, current)
    noUnsatSubsetAct = current
    eMUS, current = andTwoTseitinCnf(unsatCNF, noUnsatSubset, unsatAct, noUnsatSubsetAct, current)
    eMUS.append([current])
    
    currents = [i for i in range(primeActivators[-1] + 1, current + 1)]

    result = "p cnf {} {}\n".format(maxVar(eMUS), len(eMUS))
    result += "e " + " ".join([str(i) for i in activators]) + " 0 \n"
    result += "a " + " ".join([str(i) for i in Vars + primeActivators]) + " 0 \n"
    result += "e " + " ".join([str(i) for i in primeVars + currents]) + " 0 \n"
    for cl in eMUS:
        result += " ".join([str(l) for l in cl]) + " 0\n"
    
    return result



def parse(file):
    with open(file, "r") as f:
        lines = f.readlines()
        assert lines[0][0] == "p"
        C = []
        for line in lines[1:]:
            line = line.split(" ")
            C.append([int(i) for i in line[:-1]])
        return C
import subprocess as sp

def compute(filename, C):
    cmd = "/home/xbendik/bin/caqe/target/release/caqe --qdo {}".format(filename)
    proc = sp.Popen([cmd], stdout=sp.PIPE, shell=True)
    (out, err) = proc.communicate()
    out = out.decode("utf-8")

    print(out)
##    
##    assert "SAT" in out
##    if "UNSAT" in out:
##        cex = [int(v) for v in out.splitlines()[-1].split(" ")[1:]]
##        K = list(filter(lambda x: x in activators, cex))
##        res = []
##        for i in range(len(activators)):
##            if activators[i] in K:
##                res.append(S[i])
##        print ("UNSAT: {}".format(" ".join([str(c) for c in res])))
##    else:
##        print("SAT")



import sys
if __name__ == "__main__":
    assert len(sys.argv) > 1
    filename = sys.argv[1]
    C  = parse(filename)
    encoding = exMUS(C)
    qdimacs = "exencoded.qdimacs"
    with open(qdimacs, "w") as f:
        f.write(encoding)
    compute(qdimacs, C)
