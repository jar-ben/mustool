def negate(cl):
    return [-1 * l for l in cl]

## input: l: a list of lists of numbers
## output: a maximum number in l
def maxVar(l):
    maximum = 0
    for i in l:
        maximum = max(maximum, max([abs(j) for j in i]))
    return maximum

## input: cube: a cube, i.e. conjunction of literals, represented as list of ints
## input actLit: a literal
## output: cnf formula (written as a list of lists) representing (cube <==> actLit)
def CubeActLitEq(cube, actLit):
    res = []
    for l in cube:
        res.append([-actLit,l])
    res.append([actLit] + [-l for l in cube])
    return res

## input: clause: a clause, i.e. disjunction of literals, represented as list of ints
## input actLit: a literal
## output: cnf formula (written as a list of lists) representing (clause <==> actLit)
def ClauseActLitEq(clause, actLit):
    res = []
    for l in clause:
        res.append([actLit,-l])
    res.append([-actLit] + [l for l in clause])
    return res

def dnfToCnf(dnf):
    activators = []
    current = maxVar(dnf) + 1
    acts = []
    res = []
    for cl in dnf:
        activators.append(current)
        res += CubeActLitEq(cl,current)
        acts.append(current)
        current += 1
    res.append([-current] + activators)
    for a in activators:
        res.append([-a,current])
    acts.append(current)
    return res, acts

def cnfToCnf(cnf, current):
    activators = []
    acts = []
    res = []
    for cl in cnf:
        activators.append(current)
        res += ClauseActLitEq(cl,current)
        acts.append(current)
        current += 1

    for a in activators:
        res.append([a,-current])
    res.append([current] + negate(activators))    
    acts.append(current)
    return res, acts

def qdimacs(cnf, a, e):
    result = "p cnf {} {}\n".format(maxVar(cnf), len(cnf))
    result += "a " + " ".join([str(i) for i in a]) + "\n"
    result += "e " + " ".join([str(i) for i in e]) + "\n"
    for cl in cnf:
        result += " ".join([str(l) for l in cl]) + " 0\n"
    return result

def nvars(cls):
    vrs = []
    for cl in cls:
        for l in cl:
            if l < 0: l *= -1
            vrs += [l]
    return set(vrs)


def complement(S,C):
    cm = []
    for cl in C:
        if cl not in S: cm.append(cl)
    return cm

def prime(cl, n):
    primed = []
    for l in cl:
        if l > 0:
            primed.append(l + n)
        else:
            primed.append(l - n)
    return primed

def encode(C, S):
    Vars = list(nvars(C))
    n = len(Vars)
    primeVars = [v + n for v in Vars]
    activators = [2*n + i + 1 for i in range(len(S))]

    fst = []
    for i in range(len(S)):
        fst.append(prime(negate(S[i]),n) + [activators[i]])
    fst, actsFst = dnfToCnf(fst)
    snd = []
    for cl in complement(S,C):
        snd.append(cl)
    for i in range(len(S)):
        snd.append(S[i] + [-activators[i]])

    activationWiseSnd, actsSnd = cnfToCnf(snd,actsFst[-1] + 1)    

    actFst = actsFst[-1] #activator for the whole fst
    actSnd = actsSnd[-1] #activator for the whole snd
    mainAct = actSnd + 1 #activator for the disjunction of fst and snd
    completeCnf = fst + activationWiseSnd + [[-mainAct, actFst, actSnd], [-actFst, mainAct], [-actSnd, mainAct], [mainAct]]                           

    a = primeVars + activators
    e = Vars + actsSnd + actsFst + [mainAct]
    return qdimacs(completeCnf, a, e), activators


def parse(file):
    with open(file, "r") as f:
        lines = f.readlines()
        assert lines[0][0] == "p"
        C = []
        SIndices = []
        S = []
        for line in lines[1:]:
            line = line.split(" ")
            if line[0] != "c":
                line = [int(i) for i in line[:-1]]
                C.append(line)
            else:
                SIndices.append(int(line[1]))
                S.append(C[int(line[1])])
        return encode(C,S), SIndices

import subprocess as sp

def compute(filename, activators, S):
    cmd = "/home/xbendik/bin/cadet/cadet {}".format(filename)
    proc = sp.Popen([cmd], stdout=sp.PIPE, shell=True)
    (out, err) = proc.communicate()
    out = out.decode("utf-8")

    assert "SAT" in out
    if "UNSAT" in out:
        cex = [int(v) for v in out.splitlines()[-1].split(" ")[1:]]
        K = list(filter(lambda x: x in activators, cex))
        res = []
        for i in range(len(activators)):
            if activators[i] in K:
                res.append(S[i])
        print ("UNSAT: {}".format(" ".join([str(c) for c in res])))
    else:
        print("SAT")

import sys
if __name__ == "__main__":
    assert len(sys.argv) > 1
    filename = sys.argv[1]
    (encoded, activators), S  = parse(filename)
    qdimacs = "encoded.qdimacs"
    with open(qdimacs, "w") as f:
        f.write(encoded)
    compute(qdimacs, activators, S)
