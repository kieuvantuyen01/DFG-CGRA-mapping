from pysat.formula import CNF
from pysat.solvers import Glucose4
from pysat.solvers import Cadical103
from pypblib import pblib
pb2 = pblib.Pb2cnf()
from typing import List
import math
import graph
import dfg
import maxsat
import time
import copy
try:
    from hierarchical_mapping import solve_mapping_hierarchical
    hierarchical_available = True
except ImportError:
    hierarchical_available = False
    print("Warning: Hierarchical mapping module not available. Install NetworkX to use this feature.")
nvars = 0
all_clauses = []
fincrement = True
solver = None
CGRA = None
DFG = None
ncycles = 0
tempnodes = False
finterpret = True
freduce = True
fmaxsat = False
var = {}

try:
    import config
except ImportError:
    # Giá trị mặc định nếu không có tập tin cấu hình
    config = type('', (), {})
    config.DEFAULT_NUM_CYCLES = 52
    config.DEFAULT_NREGS = 2
    config.DEFAULT_NOPES = 1
    config.USE_HIERARCHICAL = True
    config.NUM_PARTITIONS = 2
    config.FINCREMENT = True
    config.FINTERPRET = True
    config.FREDUCE = True
    config.FMAXSAT = False
    config.DFG_FILE = "f.txt"
    config.CGRA_FILE = "e.txt"

def add_clause(clause: list):
    global all_clauses,solver
    all_clauses.append(clause)
    # print(clause)
    solver.add_clause(clause)
    # solver.add_clause

def exact1(nvar: list):
    add_clause(nvar)
    atmost(1, nvar)

def naive_atmost1(nvar :list):
    for i in range(len(nvar)):
            for j in range(i+1,len(nvar)):
                add_clause([-nvar[i],-nvar[j]])


def bimander_atmost1(nvar: list, neach):
    global nvars
    n = len(nvar)
    m = int(n / neach) + (n % neach != 0)
    # print(n,m)
    if m==1:
        bi = 0
    else:
        bi = int(math.log2(m-1)+1)
    tmp = 0
    nvar2 = []
    for i in range(m):
        nvar2 = []
        tmp2 = tmp
        for j in range(neach):
            if(tmp >= n):
                break
            nvar2.append(nvar[tmp])
            tmp += 1
        if(len(nvar2) > 1): 
            naive_atmost1(nvar2)
        tmp2_ = tmp2
        for k in range(bi):
            tmp2 = tmp2_
            bit = (1<<k)
            if bit & i:
                for j in range(neach):
                    if(tmp2 >= n):
                        break
                    add_clause([-nvar[tmp2], nvars + k + 1])
                    # print([-nvar[tmp2], var + k + 1])
                    tmp2 += 1
            else:
                for j in range(neach):
                    if(tmp2 >= n):
                        break
                    add_clause([-nvar[tmp2], -(nvars + k + 1)])
                    # print([-nvar[tmp2], -(var + k + 1)])
                    tmp2 += 1
    nvars += bi

def atmost(k:int , nvar: list):
    global nvars, nvars
    if(k==1):
        bimander_atmost1(nvar,2)
        # naive_atmost1(nvar)
        return
    # formula = []
    # nvars = pb2.encode_at_most_k(nvar, k, formula, nvars + 1)
    # for clause in formula:
    #     add_clause(clause)
    
    x = nvar.copy()
    x.insert(0,0)
    r = [[0 for _ in range(k+1)] for _ in range(len(nvar)+1)]
    for i in range(1,k+1):
        for j in range(1,i+1):
            nvars += 1
            r[i][j] = nvars
    for i in range(k+1,len(nvar)):
        for j in range(1,k+1):
            nvars += 1
            r[i][j] = nvars
    
    # 1
    for i in range (1,len(nvar)):
        add_clause([-x[i],r[i][1]])

    # 2
    for i in range (2,len(nvar)):
        for j in range(1,min(i-1,k)+1):
            add_clause([-r[i-1][j],r[i][j]])

    # 3
    for i in range(2,len(nvar)):
        for j in range(2,min(i,k)+1):
            add_clause([-x[i],-r[i-1][j-1],r[i][j]])

    # 8
    for i in range(k+1,len(nvar)+1):
        add_clause([-x[i],r[i-1][k]])

    
    """
    vCounts = []
    vLits2 = []
    # First level
    vCounts.append(nvar[0])
    nvars += 1
    vLits2 = [-nvars]
    add_clause(vLits2)

    for i in range(1, k):
        vCounts.append(nvars)

    # Subsequent levels
    for j in range(1, len(nvar)):
        x = nvar[j]
        # Prohibit overflow (sum > k)
        vLits2 = [-vCounts[k - 1], -x]
        add_clause(vLits2)
        if j == len(nvar) - 1:
            break

        for i in range(k - 1, 0, -1):
            # Compute AND of x and i-1 of previous level
            nvars += 1
            a = nvars
            add_clause([-a, x])
            add_clause([-a, vCounts[i - 1]])
            add_clause([a, -x, -vCounts[i - 1]])

            # Compute OR of it and i of previous level
            nvars += 1
            new_or = nvars
            add_clause([-a, new_or])
            add_clause([-vCounts[i], new_or])
            add_clause([-new_or, a, vCounts[i]])
            # Keep it at i of this level
            vCounts[i] = new_or

        # Compute OR of x and 0 of previous level
        nvars+= 1
        new_or0 = nvars
        add_clause([-x, new_or0])
        add_clause([-vCounts[0], new_or0])
        add_clause([-new_or0, x, vCounts[0]])
        vCounts[0] = new_or0
    """

def create_array(dim1: int, dim2: int, dim3: int):
    global nvars
    array = [[[0 for _ in range(dim3)] for _ in range(dim2)] for _ in range(dim1)]
    for k in range(dim3):
        for j in range(dim2):
            for i in range(dim1):
                nvars += 1
                array[i][j][k] = nvars
    return array

def create_control_signal_array(DFG):
    global nvars
    array = [0 for q in range(DFG.nsels)]
    for q in range(DFG.nsels):
            nvars += 1
            array[q] = nvars
    return array

def add_initial_conditions(X: List[List[List[int]]]):
    # C1,2
    for j in range(CGRA.get_nnodes()):
        for i in range(DFG.get_ndata()):
            if j in CGRA.get_nodes()['mem'] and i < DFG.get_ninputs(): 
                add_clause([X[i][j][0]])
                # print(X[i][j][0])
            else:
                add_clause([-X[i][j][0]])
                # print(-X[i][j][0])


def add_final_conditions(X: List[List[List[int]]]):
    # C3
    for i in DFG.output_ids():
        add_clause([X[i][0][-1]])
        # print(X[i][0][ncycles-1])

def add_existence_constraints(X: List[List[List[int]]], Y: List[List[List[int]]], 
                            Z: List[List[List[int]]]):
    # C4
    for k in range(1, ncycles):
        for j in range(CGRA.get_nnodes()):
            for i in range(DFG.get_ndata()):
                clause = [-X[i][j][k], X[i][j][k-1]]
                # print(CGRA.incoms[j])
                for h in CGRA.incoms[j]:
                    # print(h)
                    clause.append(Y[i][h][k])
                clause.append(Z[i][j][k])
                add_clause(clause)
                # print(clause)

def add_communication_constraints(X: List[List[List[int]]], Y: List[List[List[int]]]):
    # C5
    for k in range(1, ncycles):
        for h in range(len(CGRA.get_edges())):
            for i in range(DFG.get_ndata()):
                clause = [-Y[i][h][k]]
                for j in CGRA.get_edges()[h][0]:   
                    clause.append(X[i][j][k-1])
                add_clause(clause)
                # print(clause)

def add_calculation_constraints(X: List[List[List[int]]], Y: List[List[List[int]]], 
                              Z: List[List[List[int]]], P: List[List[List[int]]], Q: List[int]):
    # C6,7
    global nvars
    for k in range(1, ncycles):
        for j in CGRA.get_pes():
            for i in range(DFG.get_ndata()):
                tmp_var = nvars
                for z in range(len(DFG.operands[i])):
                    s = DFG.operands[i][z]
                    nvars += 1
                    for operand in s:
                        if operand < 0:
                            clause = [-nvars, P[-operand-1][j][k]]
                        else:
                            clause = [-nvars, X[operand][j][k-1]]
                            for h in CGRA.incoms[j]:
                                clause.append(Y[operand][h][k])
                        add_clause(clause)
                        # print(clause)
                    if(DFG.nexs > 0):
                        # print(i,operand)
                        excond = DFG.exconds[i][z]
                        # print(excond)
                        for e in excond.items():
                            if e[1]:
                                # print(e[0])
                                clause = [-nvars, Q[e[0]]] 
                                add_clause(clause)
                                # print(clause)
                        
                clause = [-Z[i][j][k]]
                for z in range(tmp_var+1,nvars+1):
                    clause.append(z)
                add_clause(clause)
                # print(clause)
        for j in CGRA.get_mems():
            for i in range(DFG.get_ndata()):
                clause = [-Z[i][j][k]]
                add_clause(clause)
                # print(clause)

def add_capacity_constraints(X: List[List[List[int]]], Y: List[List[List[int]]], 
                           Z: List[List[List[int]]], nregs: int, nopes: int):
    # C8,9,10
    if nregs > 0:
        for k in range(1,ncycles):
            for j in CGRA.get_pes():
                nvar = []
                for i in range(DFG.get_ndata()):
                    nvar.append(X[i][j][k])
                atmost(nregs, nvar) 
                # print(k, j, comp_type, nvar)

    # print(len(all_clauses))
    for k in range(1, ncycles):
        for h in range(len(CGRA.get_edges())):
            if CGRA.get_edges()[h][2] > 0:
                nvar = []
                for i in range(DFG.get_ndata()):
                    nvar.append(Y[i][h][k])
                atmost(CGRA.get_edges()[h][2], nvar) 
    # print(len(all_clauses))

    if nopes > 0:
        for k in range(1, ncycles):
            for j in CGRA.get_pes():
                nvar = []
                for i in range(DFG.get_ndata()):
                    nvar.append(Z[i][j][k])
                atmost(1, nvar)  

def add_block_constraints(P: List[List[List[int]]],Q: List[int],
                          X: List[List[List[int]]], Y: List[List[List[int]]]):

    global all_clauses, nvars
    for k in range(ncycles):
        for j in CGRA.get_pes():
            for ex in DFG.get_exs():
                nc = len(ex[3])
                for i in range(nc):
                    nvars_ = nvars
                    for ii in range(nc):
                        d = ex[3][ii]
                        nvars+=1
                        clause = [-nvars]
                        # 13
                        if(d < 0):
                            clause.append(P[-d-1][j][k])
                        else:
                            if k > 0:
                                clause.append(X[d][j][k-1])
                            for h in CGRA.incoms[j]:
                                clause.append(Y[d][h][k])
                        add_clause(clause)
                        # print(clause)
                        clause = [-nvars]
                        s = ex[2]
                        if ex[0]:
                            # 12
                            s+= i * nc +ii
                        else:
                            # 11
                            s+= (ii-i+nc) %nc   
                        clause.append(Q[s])
                        add_clause(clause)
                        # print(clause)
                
                    clause = [-P[i+ex[1]][j][k]]   
                    for l in range( nvars_ + 1, nvars+1):
                        clause.append(l)     
                    add_clause(clause)
                    # print(clause)

                """
                for i in range(nc):
                    clause = []
                    for ii in range(nc):
                        s = ex[2] + i * nc + ii
                        clause.append(Q[s])
                    exact1(clause)
                    if not ex[0]:
                        break
                if ex[0]:
                    for ii in range(nc):
                        clause = []
                        for i in range(nc):
                            s = ex[2] + i * nc + ii
                            clause.append(Q[s])
                        # print(clause)
                        atmost(1,clause)
                """

def add_control_signal_capacity(Q: List[int]):
    for ex in DFG.exs:
        nc = len(ex[3])
        for i in range(nc):
            clause = []
            for ii in range(nc):
                s = ex[2] + i * nc + ii
                clause.append(Q[s])
            exact1(clause)
            if not ex[0]:
                break
        if ex[0]:
            for ii in range(nc):
                clause = []
                for i in range(nc):
                    s = ex[2] + i * nc + ii
                    clause.append(Q[s])
                # print(clause)
                atmost(1,clause)
                
            
def increment(X: List[List[List[int]]], Y: List[List[List[int]]], edges_cgra : List):
    for j in range(len(X[0])):
        assumptions = []
        for i in range(len(X)):
            for k in range(len(X[0][0])):
                assumptions.append(-X[i][j][k])
        check = solver.solve(assumptions=assumptions)
        if check:
            print('remove component', j)
            # If SAT, add these constraints permanently
            for assumption in assumptions:
                solver.add_clause([assumption])
        

    for j in range(len(Y[0])):    
        assumptions = []
        for i in range(len(Y)):
            for k in range(len(Y[0][0])):
                assumptions.append(-Y[i][j][k])
        check = solver.solve(assumptions=assumptions)
        if check:
            print('remove edge', j, 'from', edges_cgra[j].u, 'to', edges_cgra[j].v)
            # If SAT, add these constraints permanently
            for assumption in assumptions:
                solver.add_clause([assumption])
        solver.solve()  

def solve_mapping(num_cycles: int, nregs :int , nopes: int):
    global all_clauses, nvars, ncycles, var
    all_clauses = []
    nvars = 0 
    ncycles = num_cycles
    start_time = time.time()  # Start timing
    X = create_array(DFG.ndata, CGRA.nnodes, num_cycles)
    Y = create_array(DFG.ndata, len(CGRA.get_edges()), num_cycles)
    Z = create_array(DFG.ndata, CGRA.nnodes, num_cycles)
    P = create_array(DFG.nexs, CGRA.nnodes, num_cycles)
    Q = create_control_signal_array(DFG)

    add_initial_conditions(X)
    add_final_conditions(X)
    add_existence_constraints(X, Y, Z)
    add_communication_constraints(X, Y)
    add_calculation_constraints(X, Y, Z, P, Q)
    add_block_constraints(P, Q, X, Y)
    add_control_signal_capacity(Q)
    add_capacity_constraints(X, Y, Z, nregs, nopes)
    setup_time = time.time() - start_time
    print(f"Setup time: {setup_time:.3f} seconds")
    def print_clauses():
        global solver
        for clause in all_clauses:
            print(clause)
    
    if fmaxsat:
        hard_clauses = all_clauses
        soft_clauses = []
        for j in CGRA.get_pes():
            nvars += 1
            for k in range(ncycles):
                for i in range(DFG.get_ndata()):
                    hard_clauses.append([-nvars, -X[i][j][k]])   
            soft_clauses.append((1,[nvars]))
            var[j] = nvars
        
        for h in range(len(CGRA.get_edges())):
            nvars += 1
            for k in range(ncycles):
                for i in range(DFG.get_ndata()):
                    hard_clauses.append([-nvars, -Y[i][h][k]])   
            soft_clauses.append((1,[nvars]))
            var[CGRA.nnodes + h] = nvars    
        
        solve_start = time.time()
        output = maxsat.MaxSAT(hard_clauses, soft_clauses)
        solve_time = time.time() - solve_start
        
        for line in output.split('\n'):
            if line.startswith('v '):
                print(f"Found solution line: {line[:50]}...")
                binary_string = line[2:].strip()
                true_vars = set()
                if " " in binary_string:
                    try:
                        for val in binary_string.split():
                            val_int = int(val)
                            if val_int > 0:  # Positive literals represent true variables
                                true_vars.add(val_int)
                    except ValueError:
                        # Not integers, try as space-separated binary values
                        for i, val in enumerate(binary_string.split()):
                            if val == '1':
                                true_vars.add(i + 1)  # 1-indexed
                else:
                    # No spaces - treat as continuous binary string
                    for i, val in enumerate(binary_string):
                        if val == '1':
                            true_vars.add(i + 1)  # 1-indexed
                for j in CGRA.get_pes():
                    if j in var and var[j] in true_vars:
                        print(f"Component {j} is removed")     
                for h in range(len(CGRA.get_edges())):
                    if CGRA.nnodes + h in var and var[CGRA.nnodes + h] in true_vars:
                        print(f"Remove edge {h} from {CGRA.get_edges()[h][0]} to {CGRA.get_edges()[h][1]} is removed")
        print(f"MaxSAT solving time: {solve_time:.3f} seconds")
        return
    
    print(f"Number of variables: {nvars}")
    print(f"Number of clauses: {len(all_clauses)}")
    
    solve_start = time.time()
    sat = solver.solve()
    solve_time = time.time() - solve_start
    
    
    if sat:
        print("SAT")
        model = solver.get_model()        
        if(fincrement):
            increment_time_start = time.time()
            print("Incremental solving")
            increment(X, Y)
            increment_time = time.time() - increment_time_start
            print(f"Incremental solving time: {increment_time:.3f} seconds")
        if finterpret:
            model = solver.get_model()
            interpret_solution(model, X, Y, Z)
    else:
        print("UNSAT")
        
    print(f"Solving time: {solve_time:.3f} seconds")

def reduce_image(image: List[List[List[int]]]):
    fimage = [[[0 for _ in range(len(image[k][j]))] for j in range(CGRA.nnodes + len(CGRA.get_edges()))] for k in range(ncycles)]

    # mark output data
    # print(DFG.output_ids)
    for idx in DFG.output_ids():
        f = False
        for i in range(len(image[ncycles - 1][0])):
            if idx == image[ncycles - 1][0][i]:
                fimage[ncycles - 1][0][i] = 1
                f = True
                break
        if not f:
            print('incomplete sol')

    # propagate backwards
    for k in range(ncycles - 1, 0, -1):
        # propagate from mems
        for j in CGRA.get_mems():
            for i in range(len(image[k][j])):
                if not fimage[k][j][i]:
                    continue
                idi = image[k][j][i]
                f = False
                
                if not tempnodes or not tempnodes[j]:
                    # mark data in mem at previous cycle if possible
                    for ii in range(len(image[k - 1][j])):
                        idii = image[k - 1][j][ii]
                        if idi == idii:
                            fimage[k - 1][j][ii] = 1
                            f = True
                            break
                    if f:
                        continue
                
                # mark data coming
                for h in CGRA.incoms[j]:
                    for ii in range(len(image[k][CGRA.nnodes + h])):
                        idii = image[k][CGRA.nnodes + h][ii]
                        if idi == idii:
                            fimage[k][CGRA.nnodes + h][ii] = 1
                            f = True
                            break
                    if f:
                        break
                if not f:
                    print('incomplete sol')
        
        # propagate from pes
        for j in CGRA.get_pes():
            for i in range(len(image[k][j])):
                if not fimage[k][j][i]:
                    continue
                idi = image[k][j][i]
                f = False
                
                if not tempnodes or not tempnodes[j]:
                    # mark data in the pe at previous cycle if possible
                    for ii in range(len(image[k - 1][j])):
                        idii = image[k - 1][j][ii]
                        if idi == idii:
                            fimage[k - 1][j][ii] = 1
                            f = True
                            break
                    if f:
                        continue
                
                # mark data coming if possible
                for h in CGRA.incoms[j]:
                    for ii in range(len(image[k][CGRA.nnodes + h])):
                        idii = image[k][CGRA.nnodes + h][ii]
                        if idi == idii:
                            fimage[k][CGRA.nnodes + h][ii] = 1
                            f = True
                            break
                    if f:
                        break
                if f:
                    continue
                # mark data required for the operation
                for s in DFG.operands[idi]:
                    ff = True
                    # check if all operands are ready
                    for d in s:
                        ff = False
                        # in the node
                        for ii in range(len(image[k - 1][j])):
                            idii = image[k - 1][j][ii]
                            if d == idii:
                                ff = True
                                break
                        if ff:
                            continue
                        # coming
                        for h in CGRA.incoms[j]:
                            for ii in range(len(image[k][CGRA.nnodes + h])):
                                idii = image[k][CGRA.nnodes + h][ii]
                                if d == idii:
                                    ff = True
                                    break
                            if ff:
                                break
                        if ff:
                            continue
                        break
                    if not ff:
                        continue
                    # mark operands
                    for d in s:
                        ff = False
                        # in the node
                        for ii in range(len(image[k - 1][j])):
                            idii = image[k - 1][j][ii]
                            if d == idii:
                                fimage[k - 1][j][ii] = 1
                                ff = True
                                break
                        if ff:
                            continue
                        # coming
                        for h in CGRA.incoms[j]:
                            for ii in range(len(image[k][CGRA.nnodes + h])):
                                idii = image[k][CGRA.nnodes + h][ii]
                                if d == idii:
                                    fimage[k][CGRA.nnodes + h][ii] = 1
                                    ff = True
                                    break
                            if ff:
                                break
                        if ff:
                            continue
                        print('incomplete sol')
                    f = True
                    break
                if not f:
                    print('incomplete sol')
        """
        # propagate bypass
        for g in range(2):
            for h in bypass:
                for i in range(len(image[k][nnodes + h[0]])):
                    if not fimage[k][nnodes + h[0]][i]:
                        continue
                    idi = image[k][nnodes + h[0]][i]
                    f = False
                    for h2 in h[1]:
                        for ii in range(len(image[k][nnodes + h2])):
                            idii = image[k][nnodes + h2][ii]
                            if idi == idii:
                                fimage[k][nnodes + h2][ii] = 1
                                f = True
                                break
                        if f:
                            break
        """
        # propagate from communication paths
        for h in range(len(CGRA.get_edges())):
            for i in range(len(image[k][CGRA.nnodes + h])):
                if not fimage[k][CGRA.nnodes + h][i]:
                    continue
                idi = image[k][CGRA.nnodes + h][i]
                f = False
                """
                if h in bypass:
                    for h2 in bypass[h]:
                        for ii in range(len(image[k][CGRA.nnodes + h2])):
                            idii = image[k][CGRA.nnodes + h2][ii]
                            if idi == idii:
                                f = True
                                break
                        if f:
                            break
                """
                if f:
                    continue
                for j in CGRA.get_edges()[h][0]:
                    for ii in range(len(image[k - 1][j])):
                        idii = image[k - 1][j][ii]
                        if idi == idii:
                            fimage[k - 1][j][ii] = 1
                            f = True
                            break
                if not f:
                    print('incomplete sol')

    # update image
    
    image_old = copy.deepcopy(image)
    for k in range(ncycles):
        for j in range(CGRA.nnodes + len(CGRA.get_edges())):
            image[k][j].clear()
            for i in range(len(image_old[k][j])):
                if fimage[k][j][i]:
                    image[k][j].append(image_old[k][j][i])
    return image
    

def interpret_solution(model: List[int], X: List[List[List[int]]], Y: List[List[List[int]]], 
                      Z: List[List[List[int]]]):

    # read model
    model_set = set(abs(x) for x in model if x > 0)
    image = [[[] for _ in range(len(CGRA.get_edges()) + CGRA.get_nnodes()) ] for _ in range(ncycles)]
    
    for k in range(ncycles):
        for j in range(CGRA.get_nnodes()):
            for i in range(DFG.get_ndata()):
                if X[i][j][k] in model_set:
                    image[k][j].append(i)
        for h in range(len(CGRA.get_edges())):
            for i in range(DFG.get_ndata()):
                if Y[i][h][k] in model_set:
                    image[k][CGRA.nnodes+h].append(i)
    if freduce:
        image = reduce_image(image)
    for k in range(ncycles):
        print(f"\nCycle {k}:")
        for j in range(len(CGRA.get_edges())):
            for i in image[k][CGRA.nnodes+j]:
                print(f"node {i} is communicated in path {j} from component {CGRA.get_edges()[j][0]} to component {CGRA.get_edges()[j][1]}")
        for j in range(CGRA.get_nnodes()):
            # print(k,j,image[k][j])
            for i in image[k][j]:
                print(f"node {i} exists in component {j}")
                
def increment(X, Y):
    for j in CGRA.get_pes():
        assumptions = []
        for k in range(ncycles):
            for i in range(DFG.get_ndata()):
                assumptions.append(-X[i][j][k])
        check = solver.solve(assumptions=assumptions)
        if check:
            print('remove component', j)
            for assumption in assumptions:
                solver.add_clause([assumption])     
        solver.solve()   
    
    for h in range(len(CGRA.get_edges())):
        assumptions = []
        for k in range(ncycles):
            for i in range(DFG.get_ndata()):
                assumptions.append(-Y[i][h][k])
        check = solver.solve(assumptions=assumptions)
        if check:
            print('remove edge', h, 'from', CGRA.get_edges()[h][0], 'to', CGRA.get_edges()[h][1])
            for assumption in assumptions:
                solver.add_clause([assumption])     
        solver.solve()

def main():
    global fincrement, solver, CGRA, DFG, fincrement, finterpret, freduce, fmaxsat
    fincrement = config.FINCREMENT
    finterpret = config.FINTERPRET
    freduce = config.FREDUCE
    fmaxsat = config.FMAXSAT
    use_hierarchical = config.USE_HIERARCHICAL
    num_partitions = config.NUM_PARTITIONS
    
    solver = Glucose4(incr=fincrement)
    
    # Read CGRA architecture
    CGRA = graph.Graph()
    CGRA.create_node("mem", "_extmem")
    CGRA.read(config.CGRA_FILE)
    
    # Read DFG
    DFG = dfg.Dfg()
    DFG.read(config.DFG_FILE)
    
    MAC = True
    DFG.gen_operands(MAC, True)
    
    num_cycles = config.DEFAULT_NUM_CYCLES
    nregs = config.DEFAULT_NREGS
    nopes = config.DEFAULT_NOPES
    
    if use_hierarchical and hierarchical_available:
        print("=== Using Hierarchical Mapping Approach ===")
        solve_mapping_hierarchical(num_cycles, nregs, nopes, num_partitions)
    else:
        print("=== Using Regular Mapping Approach ===")
        solve_mapping(num_cycles, nregs, nopes)

if __name__ == "__main__":
    main()

