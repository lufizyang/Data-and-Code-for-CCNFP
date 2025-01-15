import numpy as np
import scipy
import math
import os
import time
import psutil

from scipy.optimize import Bounds, minimize
np.set_printoptions(suppress=True)


class SLSQPScipy():
    def __init__(self, tail, head, capa, supd, func, jac=None, tol=None) -> None:
        self.tail = tail
        self.head = head
        self.capa = capa
        self.supd = supd
        self.func = func
        self.jac = jac
        self.tol = tol
        self.NodeNum = len(self.supd)
        self.EdgeNum = len(self.tail)
        self.edge2index_dict = {
            (self.tail[i],self.head[i]) : i
            for i in range(self.EdgeNum)
        }

    def GenerateCons(self):
        self.__bounds = Bounds(
            np.zeros(self.EdgeNum), self.capa
        )
        self.arcind = { i : [] for i in range(self.NodeNum) }
        self.arcond = { i : [] for i in range(self.NodeNum) }
        for i in range(self.EdgeNum):
            self.arcind[self.head[i]].append(self.tail[i])
            self.arcond[self.tail[i]].append(self.head[i])
        
        def generate_cons(ai,bi):
            def eq_cons(x):
                return (ai * x).sum() - bi
            return eq_cons
        ais = []
        jacs = []
        for i in range(self.NodeNum):
            ai = np.zeros(self.EdgeNum)
            for j in self.arcind[i]:
                ai[self.edge2index_dict[(j,i)]] = -1
            for j in self.arcond[i]:
                ai[self.edge2index_dict[(i,j)]] = 1
            ais.append(generate_cons(ai,self.supd[i]))
            jacs.append(ai.tolist())
        self.eq_cons = {
            'type': 'eq',
            'fun': lambda x: np.array([ais[i](x) for i in range(self.NodeNum)]),
            'jac': lambda x: np.array(jacs),
        }
        return self.__bounds, self.eq_cons
    
    def solve(self):
        bounds, eq_cons = self.GenerateCons()
        x0 = np.ones(self.EdgeNum) * 0.5
        res = minimize(self.func, x0, method='SLSQP', constraints=[eq_cons],bounds=bounds,options={'ftol':1e-8, 'disp': True})
        # print(np.around(res.x,2))
        return res.x


if __name__ == "__main__":
    # Paramaters Setting
    obj_type = "log" # four types of objectives constructed in this file, they are: log, sqrt, exp, mix
    
    # data load
    path = r'../test_case/'
    with open(path + "edta.txt",'r') as fine:
        tail = list(map(int, fine.readline().split('\t')))
        head = list(map(int, fine.readline().split('\t')))
        capa = list(map(float, fine.readline().split('\t')))
        fixc = list(map(float, fine.readline().split('\t')))
        unic = list(map(float, fine.readline().split('\t')))
    with open(path + "ndta.txt",'r') as finn:
        supd = list(map(float, finn.readline().split('\t')))
    capa = list(map(int, capa))
    unic = list(map(float, unic))
    fixc = list(map(float, fixc))
    supd = list(map(int, supd))


    ### objective function construction
    ## Log-type Obj. construction
    def cons_log_a(a):
        def log_a(x):
            # print(x)
            return sum([math.log(x[i] + 1, a[i]) for i in range(len(a))])
        return log_a    
    def jac_log_a(a):
        def jac_log(x):
            return np.array([
                1 / (math.log(a[i]) * (x[i] + 1)) for i in range(len(a))
            ])
        return jac_log
    
    ## Sqrt-type Obj. construction
    def cons_sqrt(a):
        def sqrt_x(x):
            return sum([pow(x[i] + 1, 1.0 / a[i]) - 1 for i in range(len(x))])
        return sqrt_x   
    def jac_sqrt(a):
        def jac_sqrt_x(x):
            return np.array([
                pow(x[i] + 1, 1.0 / a[i] - 1) / a[i] for i in range(len(a))
            ])
        return jac_sqrt_x

    ## Exp-type Obj. construction
    def cons_aexp(a):
        def aexpx(x):
            return sum([
                1.0 / (1.0 + pow(a[i], -x[i])) - 0.5 for i in range(len(a))
            ])
        return aexpx   
    def jac_exp_a(a):
        def jac_exp_x(x):
            return np.array([
                1.0 / pow(1.0 + pow(a[i], -x[i]), 2) * math.log(a[i]) * pow(a[i], -x[i]) for i in range(len(a))
            ])
        return jac_exp_x

    ## mix-type Obj. construction
    lu = [] # Auxiliary variable: Use to define the mix-type objective function.
    n1,n2,n3 = 0,0,0 # Counting the numbers of above three types obj. that used in the mix-type
    for i in range(len(tail)):
        if fixc[i] < 4:
            lu.append(0)
            n1 += 1
        elif fixc[i] > 7.5:
            lu.append(2)
            n3 += 1
        else:
            lu.append(1)
            n2 += 1

    def cons_mix_(a,lu):
        def mix_x(x):
            res = np.zeros(len(x))
            for i in range(len(x)):
                if lu[i] == 0:
                    res[i] = pow(x[i] + 1, 1 / a[i]) - 1
                elif lu[i] == 1:
                    res[i] = 1.0 / (1.0 + pow(a[i], -x[i])) - 1.0/2
                else:
                    res[i] = math.log(x[i] + 1) / math.log(a[i])
            return res.sum()
        return mix_x   
    def cons_mix_jac(a,lu):
        def jac_x(x):
            jac = np.zeros(len(x))
            for i in range(len(x)):
                if lu[i] == 0:
                    jac[i] = pow(x[i] + 1, 1 / a[i] - 1) / a[i]
                elif lu[i] == 1:
                    jac[i] = 1.0 / pow(1.0 + pow(a[i], -x[i]),2) * math.log(a[i]) * pow(a[i],-x[i])
                else:
                    jac[i] = 1 / (math.log(a[i]) * (x[i] + 1))
            return jac
        return jac_x
    
    ### Obj. setting by paras. obj_type
    if obj_type == "log":
        fucl = cons_log_a([i + 1 for i in unic])
        jacl = jac_log_a([i + 1 for i in unic])
    elif obj_type == "sqrt":
        fucl = cons_sqrt([i + 1 for i in unic])
        jacl = jac_sqrt([i + 1 for i in unic])
    elif obj_type == "exp":
        fucl = cons_aexp([i + 1 for i in unic])
        jacl = jac_exp_a([i + 1 for i in unic])
    elif obj_type == "mix":
        fucl = cons_mix_([i + 1 for i in unic],lu)
        jacl = cons_mix_jac([i + 1 for i in unic],lu)
    else:
        print("Do not set the dis-existed Obj. type!")
        exit()

    ### Output the basic info
    print("The number of node and edge are ", len(supd), ", ", len(unic))
    print("The type of objective is: ", obj_type)
    print()
    tm = psutil.Process(os.getpid()).memory_info().rss / 1024
    # SLSQP running part
    sls = SLSQPScipy(tail, head, capa, supd, fucl)
    start = time.time()
    sol_sls = sls.solve()
    end = time.time()

    # Testing all the constraints 
    print()
    if sls.eq_cons['fun'](sol_sls).all() == 0.0:
        print("All constraints have been satisfied.")
    else:
        print("There is at least one constraint that not be satisfied!")
    print("The objective of SLSQP is: ", fucl(sol_sls))
    print("The solving time of SLSQP is: ", np.around(end - start, 3), "s")
    print("The Memory Usage of SLSQP is: ", (psutil.Process(os.getpid()).memory_info().rss) / 1024 - tm, "kb")
    # with open("res_temp.txt",'w') as fo:
    #     for i in range(len(sol_sls)):
    #         fo.write(str(sol_sls[i]))
    #         fo.write(',')






