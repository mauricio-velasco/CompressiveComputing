'''
Created on Aug 28, 2015

@author: velasco
'''
import pdb
import sympy as sp
import numpy as np
import matplotlib.pyplot as plt
from scipy.special._ufuncs import binom
import random as rnd
import math as math

def multinomial(Tot, vec):
    #Computes the special multinomial function \binom{d}{\alpha} defined in the paper
    if np.sum(vec)!=Tot:
        raise(AssertionError("The multinomial is not well defined"))
    else:
        if len(vec)==1: return 1
        else:
            return(binom(Tot,vec[0])*multinomial(Tot-vec[0],np.delete(vec,[0])))

def multinomial_test():
    #TODO: Convert in test in some testing framework
    for _ in range(10):
        V=np.zeros(5)
        Den = 1
        for m in range(5):
            V[m]=rnd.randint(1,5)
            Den = Den * math.factorial(V[m])
        Tot=np.sum(V)
        Res = math.factorial(Tot)/Den
        assert (multinomial(Tot,V) == Res)
    print("MULTINOMIAL PASO")


def degree(mon):
    #Computes the degree of a multi-index
    res =0
    for k in mon:
        res+=k
    return res

def oneMoreMon(monomialList, d):
    Result=[]
    for mon in monomialList:
        t=degree(mon)
        for k in range(d+1-t):
            #pdb.set_trace()
            monNew=mon.copy()
            monNew.append(k)
            Result.append(monNew)
    return Result    

def multiIndices(n,d):
    #Computes all exponent vectors of monomials of degree AT MOST d in n variables
    Monoms=[[]]
    for k in range(n):
        Monoms=oneMoreMon(Monoms,d)
    return Monoms

def oneMoreElem(sizeSub, sizeTot, elementList):
    Result = []
    for mon in elementList:
        curSum = np.sum(mon)
        AllowOne = True
        AllowZero = True
        if curSum + sizeTot-len(mon) == sizeSub:
            AllowZero = False
        if curSum == sizeSub:
            AllowOne = False
        if AllowOne:
            monNew = mon.copy()
            monNew.append(1)
            Result.append(monNew)
        if AllowZero:
            monNew = mon.copy()
            monNew.append(0)
            Result.append(monNew)        
    return Result

def subsets(sizeSub, sizeTot, indices =False):
    Sets = [[]]
    for _ in range(sizeTot):
        Sets= oneMoreElem(sizeSub,sizeTot, Sets)
    Results = Sets
    if indices:
        Results = []
        for mySet in Sets:
            myList = [i for i, x in enumerate(mySet) if x == 1]
            Results.append(myList)    
    return Results

def sampleSumOfSquaresMonomials(n,d,sampleSize):
    Mons = multiIndices(n,d)
    NewMons=[]
    square = lambda x: 2*x
    for mon in Mons:
        mon=list(map(square,mon))
        NewMons.append(mon)
    pointSamples = uniformPointSample(n,sampleSize)
    ValuesMatrix = evaluateMonomialsAtPoints(NewMons, pointSamples )
    return (np.sum(ValuesMatrix,axis=0))


def uniformPointSample(n,sampleSize):
    #Create a sample of sampleSize points, uniformly in the unit sphere in R^n. The rows of the return matrix are the samples.
    mu,sigma = 0, 1
    SampleMatrix = np.empty((n,sampleSize))
    for sample in range(n):
        s = np.random.normal(mu,sigma,sampleSize)
        SampleMatrix[sample]=np.array(s)
    FinalMatrix = np.empty((sampleSize,n))
    for index in range(sampleSize):
        Col = SampleMatrix[:,index]
        FinalMatrix[index] = Col/np.linalg.norm(Col) 
    return FinalMatrix

def uniformPointSample_Test():
    pointSample = uniformPointSample(3, 100)
    for point in pointSample:
        assert(abs(np.linalg.norm(point)-1.0)<0.0001)

def regularPointSampleCircle(sampleSize):
    #Produces regularly spaced points in a circle
    Results = []
    for j in range(sampleSize):
        Results.append(np.array([math.cos(2*math.pi*j/sampleSize), math.sin(2*math.pi*j/sampleSize)]))
    return np.array(Results)


def evaluateMonomialsAtPoints(monomials, pointSamples):
    #Returns the measurement matrix of a given set of monomials (multi-indices) at a set of points.#The rows of the return matrix are indexed by monomials and the columns are indexed by the points
    allValues=[]
    for point in pointSamples:
        pointValues = []
        for mon in monomials:
            value = 1
            for k in range(len(mon)):
                value = value * np.power(point[k],mon[k])
            pointValues.append(value)
        allValues.append(pointValues)
    return np.transpose(np.array(allValues))


def MonomialMeasurementMatrices(m,d,pointSet,sampleSize):
    #Returns a set of sampleSize \Phi measurement matrices evaluated on the points of the pointSet
    n = len(pointSet[0])
    Mons = multiIndices(n, d)
    monValues = np.matrix(evaluateMonomialsAtPoints(Mons, pointSet)) #Rows are monomials, columns are the points    
    Results = []
    for _ in range(sampleSize):
        SampleIndexSet = np.random.choice(len(Mons), m, replace=False)
        Results.append(monValues[SampleIndexSet])
    return Results

    
def PhiMeasurementMatricesOld( m, d , pointSet, sampleSize):
    #Returns an array of randomized measurement matrices of size m x pointset length
    n = len(pointSet[0])
    Mons = multiIndices(n, d)
    monValues = np.matrix(evaluateMonomialsAtPoints(Mons, pointSet)) #Rows are monomials, columns are the points
    #Next, the coefficient matrix Rows are the m measurements, columns are the coefficients of the monomials
    Results = []
    for _ in range(sampleSize):
        index = 0
        CoeffsMatrix = np.empty((len(Mons),m))
        for mon in Mons:
            monL=mon.copy()
            pdb.set_trace()
            monL.append(d-np.sum(monL))
            mu = 0
            sigma = np.math.sqrt(multinomial(d,monL)/pow(2,d))#MUY IMPORTANTE: el input del simulador es la desviacion estandard!!  
            s = np.random.normal(mu,sigma, m)
            CoeffsMatrix[index]=np.array(s)
            index+=1
        CoeffsMatrix = np.matrix(np.transpose(CoeffsMatrix))*monValues
        Results.append((1/math.sqrt(m)) * CoeffsMatrix)
    return Results

def PhiMeasurementMatrices( m, d , pointSet, sampleSize, rescaled=False):
    #Returns an array of randomized measurement matrices of size m x pointset length
    n = len(pointSet[0])
    Mons = multiIndices(n, d)
    monValues = np.matrix(evaluateMonomialsAtPoints(Mons, pointSet)) #Rows are monomials, columns are the points
    #Next we produce a random sample for the coefficients
    cov = []
    for mon1 in Mons:
        row=[]
        for mon2 in Mons:
            if mon1==mon2:
                monL=mon1.copy()
                monL.append(d-np.sum(monL))
                row.append(multinomial(d,monL)/pow(2,d))
            else:
                row.append(0)
        cov.append(row)
    mean = np.zeros(len(Mons))
    Results = []
    for _ in range(sampleSize):
        coeffs = np.random.multivariate_normal(mean, cov,m)
        M=(1/math.sqrt(m))* np.matrix(coeffs)*monValues
        if rescaled:
            sF = computeScalingFactor(np.transpose(M)*M)
            M = math.sqrt(sF)*M
        Results.append(M)
    
    return Results
 
 
def estimateDelta(k,measurementMatrix, numSets):
    shape = np.shape(measurementMatrix)
    Results = []
    for _ in range(numSets):
        sampleIndexSet = np.random.choice(shape[1], k, replace=False)
        subMatrix = measurementMatrix[:,sampleIndexSet]
        eigv = np.linalg.eigvalsh(np.transpose(subMatrix)*subMatrix)
        Results.append(max(1-np.min(eigv), np.max(eigv)-1))
    return(max(Results))

def estimateDelta_test():
    I = np.identity(20)
    I = np.matrix(I)
    for k in range(1,10):
        print(estimateDelta(k,I,5))
    


def estimateDeltakonRangeForGivenMatrices(ListOfSizesk, measurementMatrices, numSets ):
    Xs = []
    Ys = []
    for k in ListOfSizesk:
        for curMatrix in measurementMatrices:
            Xs.append(k)
            Ys.append(estimateDelta(k,curMatrix,numSets))
    Results=[Xs,Ys]
    return Results

def firstKDeltakonRangeForGivenMatrices(ListOfSizesk, measurementMatrices, first = False):
    Xs = []
    Ys = []
    for k in ListOfSizesk:
        for curMatrix in measurementMatrices:
            Xs.append(k)
            if first:
                for set in subsets(k,3):#TODO
                    subMatrix = curMatrix[:,range(k)]
            else:
                subMatrix = curMatrix[:,range(k)]
                
            Ys.append(computeDelta(np.transpose(subMatrix)*subMatrix))
    Results=[Xs,Ys]
    return Results

def firstKDeltakonRangeForGivenMatricesRescaled(ListOfSizesk, measurementMatrices):
    Xs = []
    Ys = []
    for k in ListOfSizesk:        
        for curMatrix in measurementMatrices:
            Xs.append(k)
            subMatrix = curMatrix[:,range(k)]
            Ys.append(computeDelta(np.transpose(subMatrix)*subMatrix))
    Results=[Xs,Ys]
    return Results

def DeltakonRangeForGivenMatrices(ListOfSizesk, measurementMatrices, first=False):
    Xs = []
    Ys = []
    for k in ListOfSizesk:        
        for curMatrix in measurementMatrices:
            Xs.append(k)
            if first:
                subMatrix = curMatrix[:,range(k)]
                Ys.append(computeDelta(np.transpose(subMatrix)*subMatrix))
            else:
                deltas = []
                for mySet in subsets(k,curMatrix.shape[1],indices = True):
                    subMatrix = curMatrix[:,mySet]
                    deltas.append(computeDelta(np.transpose(subMatrix)*subMatrix))
                Ys.append(np.max(deltas))
                
    Results=[Xs,Ys]
    return Results




def PhiMeasurementMatrices_Mean(n,numPoints, m, d, k , sampleSize):
    pointSet = uniformPointSample(n, numPoints)
    Results = []
    measurementMatrices = PhiMeasurementMatrices(m, d , pointSet, sampleSize)
    for A in measurementMatrices:
        subMatrix = A[:,range(k)]
        Results.append(np.transpose(subMatrix)*subMatrix)
    
    mean = np.average(Results, axis = 0)
    V = Vmatrix(d,pointSet[range(k)])
    return [mean,V]

def Vmatrix(d, pointSet):
    M= len(pointSet)
    Result = np.zeros((M,M))
    for s in range(M):
        for t in range(M):
            Result[s,t]=np.power((np.dot(pointSet[s],pointSet[t])+1)/2,d)
    return Result

def computeDelta(squareSymmMatrix):
    eigV=np.linalg.eigvalsh(squareSymmMatrix)
    return max(1-np.min(eigV),np.max(eigV)-1)

def computeScalingFactor(squareSymmMatrix):
    #WARNING: This confers unfair advantage because if a constant 
    eigV=np.linalg.eigvalsh(squareSymmMatrix)
    return 2/(np.max(eigV)+np.min(eigV))


    

def Example1DeltaforPhiOnCircle(Ds, Ms, FigureFileName):
    # We wish to sample on a regular grid in two dimensions
    n = 2
    pointSet = regularPointSampleCircle(10)
    fig = plt.figure()
    figCounter = 1
    #We choose degree 50 for the mean to be in a somewhat interesting region
    for d in Ds:#degree of the polynomials
        for m in Ms:#number of measurements
            # We choose the k`s we wish to study
            SizesOfK = range(2,10)
            # Next, we plot the delta of the theoretical mean V
            VDeltas = []
            for k in SizesOfK:
                V=Vmatrix(d,pointSet[np.arange(k)])
                #pdb.set_trace()
                VDeltas.append(computeDelta(V))
            #Creation of picture
            sizes = np.full(len(SizesOfK), 100, dtype=int)
        
            #Next, we plot the deltas from our random sample, 50 measurement matrices at each point
            PhimeasurementMatrices = PhiMeasurementMatrices(m, d, pointSet, 20)    
            PhiDeltasData = DeltakonRangeForGivenMatrices(SizesOfK, PhimeasurementMatrices)
            #Next, we renormalize the matrices
            RescaledPhiMeasurementMatrices=[]
            for A in PhimeasurementMatrices:
                sF = computeScalingFactor(np.transpose(A)*A)
                RescaledPhiMeasurementMatrices.append(math.sqrt(sF)*A)
        
            RescaledPhiDeltasData =  DeltakonRangeForGivenMatrices(SizesOfK, PhimeasurementMatrices) 
        
            #figure we want to return
            fig.add_subplot(len(Ms),len(Ds),figCounter)
            figCounter+=1
            plt.title("m="+str(m)+" d="+str(d))
            plt.scatter(PhiDeltasData[0],PhiDeltasData[1])
            plt.scatter(SizesOfK,VDeltas, color = "red", s = sizes, alpha= 0.9)#Adds red theoretical mean of V to picture
    plt.savefig(FigureFileName)
    plt.show()


def Example2DeltavsMonforPhiOnCircle():
    # We wish to sample on a regular gr
    n = 2
    pointSet = regularPointSampleCircle(10)
    #We choose degree 50 for the mean to be in a somewhat interesting region
    d = 25
    #number of measurements
    m = 300 
    # We choose the k`s we wish to study
    SizesOfK = range(2,10)

    # Next, we plot the delta of the theoretical mean V
    VDeltas = []
    for k in SizesOfK:
        V=Vmatrix(d,pointSet[np.arange(k)])
        #pdb.set_trace()
        VDeltas.append(computeDelta(V))
    #Creation of picture
    sizes = np.full(len(SizesOfK), 100, dtype=int)

    #Next, we plot the deltas from our random sample, 
    PhimeasurementMatrices = PhiMeasurementMatrices(m, d, pointSet, 20)    
    PhiDeltasData = DeltakonRangeForGivenMatrices(SizesOfK, PhimeasurementMatrices)
    #Next, we generate the monomial matrices
    MonMeasurementMatrices = MonomialMeasurementMatrices(m,d,pointSet,20)
    MonDeltasData = DeltakonRangeForGivenMatrices(SizesOfK, MonMeasurementMatrices)
    RescaledMonMeasurementMatrices=[]
    for A in MonMeasurementMatrices:
        sF = computeScalingFactor(np.transpose(A)*A)
        RescaledMonMeasurementMatrices.append(math.sqrt(sF)*A)
    RescaledMonMeasurementsData = DeltakonRangeForGivenMatrices(SizesOfK, RescaledMonMeasurementMatrices)
    #figure we want to return
    fig = plt.figure()
    fig.add_subplot(111)
    plt.scatter(PhiDeltasData[0],PhiDeltasData[1], color ="blue")
    plt.scatter(RescaledMonMeasurementsData[0],RescaledMonMeasurementsData[1], color = "green")
    plt.scatter(SizesOfK,VDeltas, color = "red", s = sizes, alpha= 0.9)#Adds red theoretical mean of V to picture
    
    plt.show()




#TODO: Implement rescaling 2/(m+M) in the quadratic form.

"PhiMeasurementMatrices_TestMean(n, numPoints, m, d, k, sampleSize)"
#mean = PhiMeasurementMatrices_Mean(2, 6, 60, 20, 3, 10000)


#First example: Delta of randomized matrices on a uniformly spaced grid on the circle
Ms=[50,200]
Ds=[10,20]
Example1DeltaforPhiOnCircle(Ds,Ms,"RandomizedDelta.pdf")



"""
fig = plt.figure()
ax = fig.add_subplot(111)
x_points = range(0,9)
y_points = range(0,9)
p = ax.plot(x_points, y_points, 'b')
ax.set_xlabel('x-points')
ax.set_ylabel('y-points')
ax.set_title('Simple XY point plot')
fig.show()
"""
"""
x=np.arange(4)
y=np.arange(4)
xerr= x+0.1
plt.figure()
plt.errorbar(x, y, xerr=[xerr,2*xerr], yerr=0.4)
plt.title("Simplest errorbars, 0.2 in x, 0.4 in y")
plt.show()
pdb.set_trace()
"""


"""
N = 50
x = np.random.rand(N)
y = np.random.rand(N)
colors = np.random.rand(N)
area = np.pi * (15 * np.random.rand(N))**2 # 0 to 15 point radiuses

plt.scatter(x, y, s=area, c=colors, alpha=0.5)
plt.show()


multinomial_test()
pointSet = uniformSample(3, 5)
#Res = AleatMeasurementMatrices(2, 3, pointSet, 2)
Res = randomMonomialMeasurementMatrices(2, 3, pointSet, 2)
pdb.set_trace()
B = np.mean(Res,axis=0)
C=B/B[0,0]
print(C)
print(np.linalg.eigvalsh(C))
"""
"""
A = sampleSumOfSquaresMonomials(3, 10, 10000)
B = sampleSumOfSquaresMonomials(3, 7, 10000)
C = sampleSumOfSquaresMonomials(3, 5, 10000)

print(np.mean(A))
B = sampleSumOfSquaresMonomials(3, 7, 1000)
C = sampleSumOfSquaresMonomials(3, 11, 1000)
D = sampleSumOfSquaresMonomials(3, 16, 1000)
pdb.set_trace()
"""
#Hay concentracion de la medida alrededore de la media pero solo en dimension mayor o igual a tres
"""
print(A)
plt.hist(A,bins=200, histtype="stepfilled",normed=True)
plt.hist(B,bins=200, histtype="stepfilled",normed=True)
plt.hist(C,bins=200, histtype="stepfilled",normed=True)
plt.show()
"""
