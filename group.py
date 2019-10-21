
def variables(cls):
    vrs = []
    for cl in cls:
        for l in cl:
            if l < 0: l *= -1
            vrs += [l]
    return list(set(vrs))

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
    res.append([current] + [-1 * a for a in activators])    
    return res, current

## input: cl: a XOR represented as a set of literals
## output: tseitin representation of cl
def tseitinXOR(cl, current):
    if len(cl) == 1:
        return [cl], current

    fst = cl[0]
    cnf = []    
    for snd in cl[1:]:
        c, current = tseitinOnCnf([[fst, snd],[-fst, -snd]], current)
        cnf += c
        fst = current
    return cnf, current

## intput: C: set C = [c1, ..., cn] of constraints
## input: activators: set of activation literals [x1, ..., xn]
## output: tseitin cnf formula saying that C using the activators, i.e. \bigwedge_{i}(ci or -xi), is satisfiable
def sat(C, activators, current):
    assert len(C) == len(activators)
    n = len(C)
    cnf = []
    for i in range(n):
        cnf.append([-activators[i]] + [l for l in C[i]])
    print("\n\npresat\n", cnf, "\n\n")
    return tseitinOnCnf(cnf, current)

def sign(l):
    if l > 0: return 1
    else: return -1

def parseUnex(filename, activators, current):
    
    with open(filename, "r") as f:
        lines = f.readlines()
        assert lines[0][0] == "p"
        C = []
        Base = []
        for line in lines[1:]:
            line = line.split(" ")
            if line[0] == "x":
                lits = [int(i) for i in line[1:-1]]
                lits = [sign(l) * activators[abs(l)-1] for l in lits]
                cl, current = tseitinXOR(lits, current)
                Base += cl
                C += [[current]]
            else:
                cl = [int(i) for i in line[:-1]]
                if len(cl) > 0:
                    C.append([sign(l) * activators[abs(l)-1] for l in cl])
        return C, Base, current

def renderGcnf(G0, rest):
    nvariables = max(variables(G0 + rest))
    clauses = len(G0 + rest)
    groups = len(rest)

    result = "p gcnf {} {} {}\n".format(nvariables, clauses, groups)
    res2 = "p cnf {} {}\n".format(nvariables, clauses)
    for cl in G0:
        result += "{0} " + " ".join([str(l) for l in cl]) + " 0\n"
        res2 += " ".join([str(l) for l in cl]) + " 0\n"
    for i in range(len(rest)):
        result += "{{{}}} ".format(i + 1) + " ".join([str(l) for l in rest[i]]) + " 0\n"
        res2 += " ".join([str(l) for l in rest[i]]) + " 0\n"

    print(res2)

    return result

def negate(cls, base, current):
    t, tAct = tseitinOnCnf(cls, current)
    current = tAct
    

def exMUS(constraints, unex):
    C  = parse(constraints)
    Vars = variables(C)
    n = len(C)

    activators = [max(Vars) + i + 1 for i in range(n)]
    current = activators[-1]

    unexAndXor, unexBase, current = parseUnex(unex, activators, current)
    neg, negAct = tseitinOnCnf(unexAndXor, current)
    current = negAct
    G0 = unexBase + neg

    satT, satAct = sat(C, activators, current)
    current = satAct
    G0 += satT

    G0 += [[-negAct, satAct]]
    return renderGcnf(G0, [[l] for l in activators]), activators

def parse(filename):
    with open(filename, "r") as f:
        lines = f.readlines()
        assert lines[0][0] == "p"
        C = []
        for line in lines[1:]:
            line = line.split(" ")
            cl = [int(i) for i in line[:-1]]
            if len(cl) > 0:
                C.append(cl)
        return C

import subprocess as sp
def compute(gmus, acts):
    cmd = "./muser2-para -grp -comp {}".format(gmus)
    proc = sp.Popen([cmd], stdout=sp.PIPE, shell=True)
    (out, err) = proc.communicate()
    out = out.decode("utf-8")
    print(out)

import sys
if __name__ == "__main__":
    assert len(sys.argv) > 2
    filename = sys.argv[1]
    unex = sys.argv[2]
    encoding, activators = exMUS(filename, unex)
    gmus = "gmus.gcnf"
    with open(gmus, "w") as f:
        f.write(encoding)
    compute(gmus, activators)
