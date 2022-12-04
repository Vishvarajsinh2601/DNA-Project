import sys
import csv
import math
from math import factorial

import matplotlib.pyplot as plt
import gmpy2 as g
from gmpy2 import mpfr

g.get_context().precision = 1000

def printXY(name, XY, labels):
    """ Write a csv file called name+".csv" using
        data in XY organized as tuples (X0,Y0), (X1,Y1),
        ... (Xi, Yi), where each Xi and Yi are vectors of a
        common length, n.
    """
    with open(name+".csv","w") as file:
        wr = csv.writer(file)
        row = []
        for l in labels:
            row.append("x")
            row.append(l)
        wr.writerow(row)
        row = []
        l = 0
        for x,y in XY:
            l = max(l, len(x))
        for i in range(l):
            row = []
            for x,y in XY:
                if i < len(x):
                    row.append(x[i])
                    row.append(y[i])
                else:
                    row.append(' ')
                    row.append(' ')
            wr.writerow(row)

    
def makeCSV(name,X,X_label,Ys,Ys_labels):
    """ Write a csv file called name+".csv" using
        X as independent variable and Ys as a set of dependent variables.
        xlabel is the label for X, and ylabels for Y.
        CSV file contains lines like this, where m is length of X and Y and n is length
        of each member of Y

    """
    with open(name+'.csv', 'w') as file:
        wr = csv.writer(file)
        wr.writerow([X_label]+Ys_labels)
        for i,x in enumerate(X):
            row = [x]
            for y in Ys:
                row.append("{:1.5e}".format(y[i]))
            wr.writerow(row) 


def plot(X,Ys,labels,Xlabel="",Ylabel="",title=""):
    """ function plotting Y vs X with corresponding labels. """
    for Y,label in zip(Ys,labels):
        plt.plot(X,Y,label=label)
    plt.Xlabel(Xlabel)
    plt.Ylabel(Ylabel)
    plt.suptitle(title)
    plt.legend()
    plt.show()

def exp_func(base,e):
    # exponential function
    if base==mpfr("0.0"):
        if e==mpfr(0):
            return mpfr("1")
        else:
            return mpfr("0")
    return g.exp(g.mul(g.log(mpfr(base)),mpfr(e)))


def PMF(q, N, k):
    """ Binomial PMF. """

    tmp = g.mul(exp_func(q,k),g.mul(g.bincoef(N,k),exp_func(1-q,N-k)))
    return tmp

def PDF(q, N, k):
    """ Binomial CDF """

    tmp_list = [mpfr("0")]
    for i in range(0,k+1):
        tt1 = g.mul(exp_func(q,i),g.mul(g.bincoef(N,i),exp_func(1-q,N-i)))
        tmp_list.append( tt1 )   
    tmp1 = g.fsum(tmp_list)
    return tmp1

def plotPMF(N):
    X = [ _ for _ in range(0,N+1,1) ]
    Y = [
        [ PMF(0.01,N,_) for _ in X ],
        [ PMF(0.25,N,_) for _ in X ],
        [ PMF(0.5,N,_) for _ in X ],
        [ PMF(0.7,N,_) for _ in X ],
        [ PMF(0.9,N,_) for _ in X ]
         ]
    plot(X,Y,[ "q={} N={}".format(qi,N) for qi in [0.01,0.25,0.5,0.7,0.9] ],
         "Random Variable", "Probability", "Binomial PMF with N={} and varying q".format(N))


def plotPDF(N):
    X = [ _ for _ in range(0,N+1,1) ]
    Y = [
        [ PDF(0.01,N,_) for _ in X ],
        [ PDF(0.25,N,_) for _ in X ],
        [ PDF(0.5,N,_) for _ in X ],
        [ PDF(0.7,N,_) for _ in X ],
        [ PDF(0.9,N,_) for _ in X ]
         ]

    plot(X,Y,[ "q={} N={}".format(qi,N) for qi in [0.01,0.25,0.5,0.7,0.9] ], "Random Variable", "Probability",
         "Binomial CDF with N={} and varying q".format(N))

def plot_binomials(N):

    plotPDF(N)
    plotPMF(N)

class codeward:
    # modelling a codeward
    def __init__(self,length,num,es=1e-3):
        self.length = length
        self.num = num
        self.es = es

    def P_error(self,k=1):
        """ Compute the likelihood of a codeword being in error based on es and length """
        return g.sub(mpfr("1"),PDF(self.es,self.length,k-1))

class RSCodec:

    rs_code_objects = {}

    @staticmethod
    def create(n,k,d,es=1e-3,ee=1e-3):
        """  RS[n,k,d] instance
        """
        if (n,k,d,es,ee) in RSCodec.rs_code_objects:
            return RSCodec.rs_code_objects[(n,k,d,es,ee)]
        else:
            rs = RSCodec(n,k,d,es,ee)
            RSCodec.rs_code_objects[(n,k,d,es,ee)] = rs
            return rs

    def __init__(self,n,k,d,es=1e-3,ee=1e-3):
        """ init of RS[n,k,d] having symbol error rate(es) and error rate(ee)
        """
        self.q = 4
        self.n = n
        self.k = k
        self.d = d        
        self.t = int((d-1)/2)
        self.symbol_err_rate = es
        self.erasure_err_rate = ee
        self.result = mpfr("0")
        self.has_result = False
    def label(self):
        return "RS[{},{},{}] psym={:1.2e} perasure={:1.2e}".format(self.n,self.k,self.d,self.symbol_err_rate, self.erasure_err_rate)
    
    def R(self):
        return self.k / self.n

    def R_index(self, M):
        return (self.k - math.log(M,256))/float(self.n)
    
    def P_symbol_error(self):
        return g.sub(mpfr("1"),PDF(self.symbol_err_rate,self.n,self.t))
        
    def P_random_errors(self,i):
        return PMF(self.symbol_err_rate,self.n,i)
        
    def P_random_erasures(self,x):
        return g.sub(mpfr("1"),PDF(self.erasure_err_rate,self.n,x))
    
    def P_erasure_error(self):
        tmp = mpfr('0')
        for i in range(self.t+1):
            tmp = g.add(tmp, g.mul(self.P_random_errors(i),self.P_random_erasures(2*(self.t-i))))
        return tmp
    
    def P_result(self):
        if not self.has_result:
            self.result = g.add(self.P_symbol_error(),self.P_erasure_error())
            self.has_result = False
        return self.result

class RSCodeInnerOuter:

    def __init__(self,n_inner,k_inner,d_inner,n_outer,k_outer,d_outer,p_sub,p_strand_loss):
        self.rs_inner = RSCodec.create(n_inner,k_inner,d_inner,p_sub,mpfr("1e-15"))
        self.rs_outer = RSCodec.create(n_outer,k_outer,d_outer,self.rs_inner.P_result(),mpfr(p_strand_loss))

    def P_result(self):
        return self.rs_outer.P_result()

    def R_index(self, M):
        return (self.rs_outer.k)*(self.rs_inner.k - math.ceil(math.log(M, 256))) / (self.rs_inner.n*self.rs_outer.n)
    
    def R_raw(self):
        return (self.rs_inner.k * self.rs_outer.k)/(self.rs_inner.n*self.rs_inner.n)

    def getLabel(self):
        return "[{}*{},{}*{}]".format(self.rs_outer.n,self.rs_inner.n,self.rs_outer.k,self.rs_inner.k)

def result(L,cw_size=4,cw_er=1e-3,cw_es=1e-3,dropout=1e-3):
    n = int(L/cw_size)
    errs = []
    D = []
    for d in range(3,n,2):
        k = n - d
        D.append(d)
        rs_inner = RSCodec(n*cw_size,k*cw_size,d,cw_er,cw_es)
        rs_outer = RSCodec(255,255-31,31,rs_inner.P_result(),dropout)
        errs.append(g.log(rs_outer.P_result()))
    plot (D,[errs],["p"])    
    return


def RSC_sweep_length(n=255,k=[235],Length=[_ for _ in range(50,201,25)],p_se=[exp_func(10,_/10.0) for _ in range(-80,-10,1)],copies=[1],PL=lambda L, c: exp_func(g.mul(mpfr(1e-3),mpfr(L)),mpfr(c)), filename="sweep_length" ):

    X = Length
    Y = []
    Label = []
    for ki in k:
        for c in copies:
            for p in p_se:
                y = []
                Label.append( "RS(N={},k={}) perr={:1.1e} c={}".format(n,ki,float(p),c) )
                for l in Length:
                    rs_code = RSCodec(n,ki,n-ki,mpfr(p),PL(mpfr(l),c))
                    y.append(g.log10(rs_code.P_result()))
                Y.append(y)

    fig, ax = plt.subplots()
    plt.ylabel("Probabilty Decoding Error - Log Scale (10^y)")
    plt.xlabel("Strand Length (nt)")
    for y,label in zip(Y,Label):
        plt.plot(X,y,label=label)
    plt.legend()
    plt.show()
    makeCSV(filename,X,"Length",Y,Label)


def innerTarget(L,index,P_target,p_sym):
    """ Look for an inner RSCodec that is at least as excellent as P target, given L and
        p_sym's symbol error rate.
    """
    cw = codeward(4,256,p_sym)
    n = int(L/4)
    inner = None
    guess = int(n*p_sym)
    if guess % 2 == 0:
        guess += 1
    for d in range(guess,n-3,2):
        if n-d-index < 1:
            continue
        inner = RSCodec.create(n,n-d,d,es=mpfr(cw.P_error()),ee=mpfr(1e-15))
        if inner.P_result() <= P_target:
            return inner
    return inner

def getRS(M=1e9,L=200,P_res=mpfr(1e-11),p_sym=mpfr(1e-3),p_loss=mpfr(1e-3)):
    cw = codeward(4,256,p_sym)
    index = math.log(M,256)
    inner = innerTarget(L,index,P_res,p_sym)
    if inner==None:
        return None
    frac = g.factorial(int((inner.d-1)/2))
    outer = outerTarget(P_res,g.div(inner.P_result(),frac),
                                  g.mul(inner.P_result(),g.sub(mpfr(1),g.div(mpfr(1),mpfr(frac))))+p_loss)
    if outer==None:
        return None
    return RSCodeInnerOuter(inner.n,inner.k,inner.d,outer.n,outer.k,outer.d,cw.P_error(),p_loss)

def outerTarget(P_target,p_sym,p_loss):
    """ By adjusting the distance, look for an outer RSCodec that is at least as good as P
        target, assuming L=255, a symbol error rate of p sym, and a strand loss rate of p loss.
    """
    n=255
    outer = None
    guess = int(n*p_loss)
    if guess % 2 == 0:
        guess += 1
    for d in range(guess,n-3,2):
        outer = RSCodec.create(n,n-d,d,p_sym,p_loss)
        if outer.P_result() < P_target:
            return outer
    return outer
    

def RS_rate(L=[l for l in range(50,1000,20)],
                    pse=[1e-3], slope=[1e-3],
                    PL=lambda L, c, slope: exp_func(g.mul(mpfr(slope),mpfr(L)),mpfr(c))):

    XY = []
    XRE = []
    Label = []
    for s in slope:
        for p in pse:
            y = []
            re = []
            x = []
            Label.append( "RS p_break/nt={:1.1e} p_sym={:1.1e}".format(float(s),float(p)) )
            for l in L:
                rs = getRS(1e9,l,mpfr("1e-14"),p,PL(l,1,s))
                if rs == None or rs.P_result() > mpfr("1e-14"):
                    continue
                x.append(l)
                y.append(rs.R_index(1e9))
                re.append(g.log10(rs.P_result()))
                
            XY.append( [x,y] )
            XRE.append( [x,re] )
    
    fig, ax = plt.subplots()
    plt.ylabel("Information Density")
    plt.xlabel("Strand Length (nt)")
    for [x,y],label in zip(XY,Label):
        plt.plot(x,y,label=label)
    plt.legend()
    plt.show()
    printXY("RS_rate",XY,Label)

def infoDensityVSloss(L=[l for l in range(50,1000,20)],
                    pse=[1e-3], slope=[1e-3],
                    PL=lambda L, c, slope: exp_func(g.mul(mpfr(slope),mpfr(L)),mpfr(c))):

    XY = []
    XRE = []
    Label = []
    for l in L:
        for p in pse:
            y = []
            re = []
            x = []
            Label.append( "RS L={} p_sym={:1.1e}".format(l,float(p)) )
            for s in slope:
                rs = getRS(1e9,l,mpfr("1e-14"),p,PL(l,1,s))
                if rs == None or rs.P_result() > mpfr("1e-14"):
                    continue
                x.append(g.log10(s))
                y.append(rs.R_index(1e9))
                re.append(g.log10(rs.P_result()))
                
            XY.append( [x,y] )
            XRE.append( [x,re] )
    
    fig, ax = plt.subplots()
    plt.ylabel("Information Density")
    plt.xlabel("Log Probability of Strand Loss")
    ax.set_xlim(-2.5,-8.5)
    for [x,y],label in zip(XY,Label):
        plt.plot(x,y,label=label)
    plt.legend()
    plt.show()
    printXY("compare_density_loss",XY,Label)

def RSCode_copy(n,d,copy,p_sym,p_erasure):
    loss = exp_func(p_erasure,copy)
    rs = RSCodec(n,n-d+1,d,p_sym,loss)
    return rs
    
    
def comparedCopyandCode(p_err):
    label = []
    error = []
    rate = []
    rs = RSCode_copy(255,3,4,mpfr("0"),mpfr(p_err))
    
    label.append( rs.label() + " copies=4 " )
    error.append( g.log10(rs.P_result()) )
    rate.append( rs.R() / 4 )

    rs =outerTarget(rs.P_result(),mpfr("0"),mpfr(p_err))
    label.append( rs.label() + " copies=1" )
    error.append( g.log10(rs.P_result()) )
    rate.append( rs.R() )

    for l,e,r in zip(label,error,rate):
        print ("{}: {} {}".format(l,e,r))
    
def plotting_figures():
    RSC_sweep_length(n=255,
                k=[223],Length=[_ for _ in range(50,1025,25)],
                p_se=[mpfr(1e-2),mpfr(1e-3)],
                copies = [1,3,10],
                PL=lambda l,c: exp_func(g.mul(mpfr("0.0005"),mpfr(l)),mpfr(c)))

    RS_rate(slope=[1e-3,5e-4,1e-4],L=[l for l in range(50,1020,20)])
    infoDensityVSloss(slope=[ 10**(_/10.0) for _ in range(-80,-9, 10) ],L=[l for l in range(200,1100,200)])

if __name__ == "__main__":
    plotting_figures()
