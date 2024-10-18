# Incremental Search!!! SAT
import time
import random
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

# Following is an example of a wff with 3 variables, 3 literals/clause, and 4 clauses
Num_Vars=3
Num_Clauses=4
wff=[[1,-2,-2],[2,3,3],[-1,-3,-3],[-1,-2,3],[1,2,-3]]


# Following is an example of a wff with 3 variables, 3 literals/clause, and 8 clauses
Num_Clauses=8
wff=[[-1,-2,-3],[-1,-2,3],[-1,2,-3],[-1,2,3],[1,-2,-3],[1,-2,3],[1,2,-3],[1,2,3]]

# NEW CODE BELOW!!!
def check(wff, assignment, var_index):
    # check if the current assignment satisfies the wff (at least partially)
    for clause in wff:
        clause_sat = False
        unchecked_lits = False
        for literal in clause:
            var = abs(literal)
            if var < var_index:
                if (literal > 0 and assignment[var] == 1) or (literal < 0 and assignment[var] == 0):
                    clause_sat = True
                    break
            else:
                unchecked_lits = True
        if not clause_sat and not unchecked_lits:
            return False  # unsatisfiable, prune branch
    return True

def solve_recursive(wff, assignment, var_index, Nvars, Nclauses):
    if var_index > Nvars:
        return check(wff, assignment, var_index)
    
    if not check(wff, assignment, var_index):
        return False  # prune this branch

    # try true on current variable
    assignment[var_index] = 1
    if solve_recursive(wff, assignment, var_index + 1, Nvars, Nclauses):
        return True
    
    # try false on current variable
    assignment[var_index] = 0
    if solve_recursive(wff, assignment, var_index + 1, Nvars, Nclauses):
        return True
    
    return False  # backtrack if neither works
# NEW CODE ABOVE!!!

def build_wff(Nvars,Nclauses,LitsPerClause):
    wff=[]
    for i in range(1,Nclauses+1):
        clause=[]
        for j in range(1,LitsPerClause+1):
            var=random.randint(1,Nvars)
            if random.randint(0,1)==0: var=-var
            clause.append(var)
        wff.append(clause)
    return wff

def test_wff(wff,Nvars,Nclauses):
    Assignment=list((0 for x in range(Nvars+2)))
    start = time.time() # Start timer
    SatFlag=solve_recursive(wff, Assignment, 1, Nvars, Nclauses) # new function
    end = time.time() # End timer
    exec_time=int((end-start)*1e6)
    return [wff,Assignment,SatFlag,exec_time]

def run_cases(TestCases,ProbNum,resultsfile,tracefile,cnffile):
    # TestCases: list of 4-tuples describing problem
    #   0: Nvars = number of variables
    #   1: NClauses = number of clauses
    #   2: LitsPerClause = Literals per clause
    #   3: Ntrials = number of trials
    # ProbNum: Starting number to be given to 1st output run
    # resultsfile: path to file to hold output
    # tracefile: path to file to hold output
    # cnffile: path to file to hold output
    # For each randomly built wff, print out the following list
    #   Problem Number
    #   Number of variables
    #   Number of clauses
    #   Literals per clause
    #   Result: S or U for satisfiable or unsatisfiable
    #   A "1"
    #   Execution time
    #   If satisfiable, a binary string of assignments
    if not(ShowAnswer):
        print("S/U will NOT be shown on cnf file")
    f1=open(resultsfile+".csv",'w')
    f2=open(tracefile+".csv",'w')
    f3=open(cnffile+".cnf","w")
    #initialize counters for final line of output
    Nwffs=0
    Nsat=0
    Nunsat=0

    for i in range(0,len(TestCases)):
        TestCase=TestCases[i]
        Nvars=TestCase[0]
        NClauses=TestCase[1]
        LitsPerClause=TestCase[2]
        Ntrials=TestCase[3]
        #Now run the number of trials for this wff configuration
        Scount=Ucount=0
        AveStime=AveUtime=0
        MaxStime=MaxUtime=0
        for j in range(0,Ntrials):
            #generate next trial case for this configuration
            Nwffs=Nwffs+1
            random.seed(ProbNum)
            wff = build_wff(Nvars,NClauses,LitsPerClause)
            results=test_wff(wff,Nvars,NClauses)
            wff=results[0]
            Assignment=results[1]
            Exec_Time=results[3]
            if results[2]:
                sat_vars.append(Nvars)
                sat_times.append(Exec_Time)
                y='S'
                Scount=Scount+1
                AveStime=AveStime+Exec_Time
                MaxStime=max(MaxStime,Exec_Time)
                Nsat=Nsat+1
            else:
                unsat_vars.append(Nvars)
                unsat_times.append(Exec_Time)
                y='U'
                Ucount=Ucount+1
                AveUtime=AveUtime+Exec_Time
                MaxUtime=max(MaxUtime,Exec_Time)
                Nunsat=Nunsat+1
            x=str(ProbNum)+','+str(Nvars)+','+str(NClauses)+','+str(LitsPerClause)
            x=x+str(NClauses*LitsPerClause)+','+y+',1,'+str(Exec_Time)
            if results[2]:
                for k in range(1,Nvars+1):
                    x=x+','+str(Assignment[k])
            print(x)
            f1.write(x+'\n')
            f2.write(x+'\n')
            #Add wff to cnf file
            if not(ShowAnswer):
                y='?'
            x="c "+str(ProbNum)+" "+str(LitsPerClause)+" "+y+"\n"
            f3.write(x)
            x="p cnf "+str(Nvars)+" "+str(NClauses)+"\n"
            f3.write(x)
            for i in range(0,len(wff)):
                clause=wff[i]
                x=""
                for j in range(0,len(clause)):
                    x=x+str(clause[j])+","
                x=x+"0\n"
                f3.write(x)
            #Increment problem number for next iteration
            ProbNum=ProbNum+1
        counts='# Satisfied = '+str(Scount)+'. # Unsatisfied = '+str(Ucount)
        maxs='Max Sat Time = '+str(MaxStime)+'. Max Unsat Time = '+str(MaxUtime)
        aves='Ave Sat Time = '+str(AveStime/Ntrials)+'. Ave UnSat Time = '+str(AveUtime/Ntrials)
        print(counts)
        print(maxs)
        print(aves)
        f2.write(counts+'\n')
        f2.write(maxs+'\n')
        f2.write(aves+'\n')
    x=cnffile+",TheBoss,"+str(Nwffs)+","+str(Nsat)+","+str(Nunsat)+","+str(Nwffs)+","+str(Nwffs)+"\n"
    f1.write(x)
    f1.close()
    f2.close()
    f3.close()

# Following generates several hundred test cases of 10 different wffs at each size
# and from 4 to 22 variables, 10 to 240 clauses, and 2 to 10 literals per clause 
TestCases=[
    [4,10,2,10],
    [8,16,2,10],
    [12,24,2,10],
    [16,32,2,10],
    [18,36,2,10],
    [20,40,2,10],
    [22,44,2,10],
    [24,48,2,10],
    [4,20,3,10],
    [8,40,3,10],
    [12,60,3,10],
    [16,80,3,10],
    [18,90,3,10],
    [20,100,3,10],
    [22,110,3,10],
    [24,120,3,10],
    [4,40,4,10],
    [8,80,4,10],
    [12,120,4,10],
    [16,160,4,10],
    [18,180,4,10],
    [20,200,4,10],
    [22,220,4,10],
    [24,240,4,10],
    [4,40,5,10],
    [8,80,5,10],
    [12,120,5,10],
    [16,160,5,10],
    [18,180,5,10],
    [20,200,5,10],
    [22,220,5,10],
    [24,240,5,10],
    [4,40,6,10],
    [8,80,6,10],
    [12,120,6,10],
    [16,160,6,10],
    [18,180,6,10],
    [20,200,6,10],
    [22,220,6,10],
    [24,240,6,10]]

trace=True
ShowAnswer=True # If true, record evaluation result in header of each wff in cnffile
ProbNum = 3
resultsfile = r'resultsfile'
tracefile = r'tracefile'
cnffile = r'cnffile' # Each of these list entries describes a series of random wffs to generate

#run_cases(TC2,ProbNum,resultsfile,tracefile,cnffile)
#run_cases(SAT2,ProbNum,resultsfile,tracefile,cnffile)

# NEW CODE BELOW!!!
sat_vars = []
sat_times = []
unsat_vars = []
unsat_times = []
run_cases(TestCases,ProbNum,resultsfile,tracefile,cnffile)

# finding curves
def exp_func(x, a, b):
    return a * np.exp(b * x)
sat_popt, sat_pcov = curve_fit(exp_func, sat_vars, sat_times)
unsat_popt, unsat_pcov = curve_fit(exp_func, unsat_vars, unsat_times)

sat_curve_text = "SAT curve: " + str(round(sat_popt[0], 3)) + " * e^(" + str(round(sat_popt[1], 3)) + "n)"
unsat_curve_text = "UNSAT curve: " + str(round(unsat_popt[0], 3)) + " * e^(" + str(round(unsat_popt[1], 3)) + "n)"

plt.text(0.05, 0.95, sat_curve_text, transform=plt.gca().transAxes, fontsize=10, verticalalignment='top', color='green')
plt.text(0.05, 0.90, unsat_curve_text, transform=plt.gca().transAxes, fontsize=10, verticalalignment='top', color='red')
plt.scatter(sat_vars, sat_times, color='green', label='SAT', marker='o')
plt.scatter(unsat_vars, unsat_times, color='red', label='UNSAT', marker='x')
plt.xlabel('# of Variables')
plt.ylabel('Execution Time (µs)')
plt.title('SAT INCREMENTAL SOLVER DATA')
plt.legend()
plt.savefig("incr_SAT_figure.png")
plt.show()
