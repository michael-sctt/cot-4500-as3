# Michael Scott
# COT4500 Assignment 3

# Imports
import numpy as np

# Problem 1 Function
def EulerMethod(funct, rangeStart, rangeStop, iters, initPoint):

    # Calculate step size
    stepSize = (rangeStop - rangeStart) / iters

    # Make array for calculated values
    sols = []
    sols.append(initPoint)
    
    # Retrieve initial t and y
    t = 0
    y = sols[0]

    # Calculate values
    for i in range(0, iters + 1):
        sols.append(sols[i] + (stepSize * eval(funct)))
        t += stepSize
        y = sols[i + 1]

    # Last value calculated is the one we're looking for, print that out (to 5 decimal places)
    print("%.5f" % sols[-1])

# Problem 2 Function
def RungeKutta(funct, rangeStart, rangeStop, iters, initPoint):

    # Calculate step size
    stepSize = (rangeStop - rangeStart) / iters
    timeSteps = np.arange(rangeStart, rangeStop + stepSize, stepSize)

    # Make array for calculated values
    sols = []
    sols.append(initPoint)
    
    # Retrieve initial t and y
    t = timeSteps[0]
    y = sols[0]

    # Calculate values
    for i in range(0, iters):
        k1 = stepSize * eval(funct)
        
        t = t + (stepSize / 2)
        y = y + (k1 / 2)
        k2 = stepSize * eval(funct)

        y = y - (k1 / 2) + (k2 / 2)
        k3 = stepSize * eval(funct)

        t = timeSteps[i + 1]
        y = y - (k2 / 2) + k3
        k4 = stepSize * eval(funct)

        y = y - k3
        sols.append(y + (1 / 6) * (k1 + 2*k2 + 2*k3 + k4))
        y = sols[i + 1]

    # Last value calculated is the one we're looking for, print that out (to 5 decimal places)
    print("%.5f" % sols[-1])

# Problem 3 Function
def GaussWithBackSub(a):
    
    # Step 1
    for k in range(0, len(a)):
        for j in range(k + 1, len(a)):
            lam = a[j][k] / a[k][k]
            for i in range(0, len(a[k])):
                a[j][i] -= lam * (a[k][i])

    length = len(a)
    x = np.zeros(length)
    x[length - 1] = np.trunc(a[length-1][length] / a[length-1][length-1])
    for i in range(length - 2, -1, -1):
        x[i] = a[i][length]
        for j in range(i + 1, length):
            x[i] = x[i] - a[i][j] * x[j]
    
        x[i] = np.trunc(x[i] / a[i][i])

    # Output coming out backwards, unknown why
    x = np.flip(x)
    print(x)

# Problem 4 Function
def LUDecomp(mat):
    length = len(mat)
    L = np.zeros((length, length))
    U = np.zeros((length, length))

    # Make LU Matrices
    for i in range(0, length):
        L[i][i] = 1

        # Make Upper Triangle Matrix
        for j in range(i, length):
            U[i][j] = mat[i][j]
            for k in range(0, i):
                U[i][j] -= L[i][k] * U[k][j]

        # Make Lower Triangle Matrix
        for j in range(i + 1, length):
            L[j][i] = mat[j][i]
            for k in range(0, i):
                L[j][i] -= L[j][k] * U[k][i]
            L[j][i] /= U[i][i]

    print(L)
    print()
    print(U)
        
# Problem 5 Function
def detDiagDom(mat):
    for i in range(0, len(mat)):
        maxVal = mat[i][i]
        sum = 0
        for j in range(0, len(mat[i])):
            if i == j:
                continue
            sum += abs(mat[i][j])
            if sum > maxVal:
                return False
            
    # If code reaches this point, matrix is diagonally dominant
    return True

# Problem 6 Function
def detPosDef(mat):

    # Check if transpose is equal to input
    if (np.array_equal(mat, np.transpose(mat))):
    
        # Check if matrix is symmetric
        eigenvalues = np.linalg.eig(mat)[0]
        for val in eigenvalues:
            if val < 0:
                return False
        return True

    return False

# Driver Code
# Problem 1
funct = "t - y**2"
EulerMethod(funct, 0, 2, 10, 1)
print()

# Problem 2 (utilizes same function as problem 1)
RungeKutta(funct, 0, 2, 10, 1)
print()

# Problem 3
a = [[2, -1, 1, 6], [1, 3, 1, 0], [-1, 5, 4, 3]]
GaussWithBackSub(a)
print()

# Problem 4
mat4 = [[1, 1, 0, 3], [2, 1, -1, 1], [3, -1, -1, 2], [-1, 2, 3, -1]]
det = np.linalg.det(mat4)  # Very useful numpy function!
print("%.5f" % det)
print()
LUDecomp(mat4)
print()

# Problem 5
mat5 = [[9, 0, 5, 2, 1], [3, 9, 1, 2, 1], [0, 1, 7, 2, 3], [4, 2, 3, 12, 2], [3, 2, 4, 0, 8]]
print(detDiagDom(mat5))
print()

# Problem 6
mat6 = [[2, 2, 1], [2, 3, 0], [1, 0, 2]]
print(detPosDef(mat6))
