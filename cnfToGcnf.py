
def variables(cls):
    vrs = []
    for cl in cls:
        for l in cl:
            if l < 0: l *= -1
            vrs += [l]
    return list(set(vrs))

def renderCnf(cls):
    nvariables = max(variables(cls))
    clauses = len(cls)
    result = "p cnf {} {}\n".format(nvariables, clauses)
    for cl in cls:
        result += " ".join([str(l) for l in cl]) + " 0\n"
    return result


def renderGcnf(G0, rest):
    nvariables = max(variables(G0 + rest))
    clauses = len(G0 + rest)
    groups = len(rest)

    result = "p gcnf {} {} {}\n".format(nvariables, clauses, groups)
    for cl in G0:
        result += "{0} " + " ".join([str(l) for l in cl]) + " 0\n"
    for i in range(len(rest)):
        result += "{{{}}} ".format(i + 1) + " ".join([str(l) for l in rest[i]]) + " 0\n"

    return result

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
        

def checkSat(S):
    filepath = "checksat.cnf"
    with open(filepath, "w") as f:
        f.write(renderCnf(S))
    cmd = "z3 {}".format(filepath)
    proc = sp.Popen([cmd], stdout=sp.PIPE, shell=True)
    (out, err) = proc.communicate()
    out = out.decode("utf-8")
    return not "unsat" in out

def collectCrits(C, mus):
    crits = []
    for c in mus:
        toCheck = []
        for i in range(len(C)):
            if i != c:
                toCheck.append(C[i])
        if checkSat(toCheck):
            crits.append(c)
    return crits

def convert(filename, target):
    C = parse(filename)
    mus = getMUS(filename)
    crits = collectCrits(C, mus)
    with open(target, "w") as f:
        f.write(renderGcnf([C[i] for i in range(len(C)) if i in crits], [C[i] for i in range(len(C)) if i not in crits]))



import subprocess as sp
def getMUS(filepath):
    cmd = "./muser2-para -comp {}".format(filepath)
    proc = sp.Popen([cmd], stdout=sp.PIPE, shell=True)
    (out, err) = proc.communicate()
    out = out.decode("utf-8")
    assert "UNSATISFIABLE" in out
    reading = False
    for line in out.splitlines():
        if "UNSATISFIABLE" in line:
            reading = True
            continue
        if reading:
            line = line.split(" ")
            assert line[0] == "v"
            return [int(l) - 1 for l in line[1:-1]]


import sys
if __name__ == "__main__":
    assert len(sys.argv) > 1
    filename = sys.argv[1]
    target = filename + ".gcnf"
    if len(sys.argv) > 2:
        target = sys.argv[2]
    convert(filename, target)
