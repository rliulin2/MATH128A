import math
import cmath
import sympy
import numpy as np
import matplotlib.pyplot as plt

def bisectionMethod(f, tupleRange, tolerance, maxIter):
    leftSide = tupleRange[0]; rightSide = tupleRange[1]
    f_left = f(leftSide)
    f_right = f(rightSide)
    if abs(f_left) <= tolerance:
        return leftSide
    elif abs(f_right) <= tolerance:
        return rightSide
    else:
        midpt = (rightSide + leftSide) / 2
        diff = (rightSide - leftSide) / 2
        f_mid = f(midpt)
        iter = 1
        while f_mid != 0 and diff > tolerance and iter < maxIter:
            iter += 1
            signLeftSide = int(math.copysign(1, f_left))
            signMidpt = int(math.copysign(1, f_mid))
            if signLeftSide * signMidpt < 0:
                rightSide = midpt
                f_right = f_mid
            else:
                leftSide = midpt
                f_left = f_mid
            midpt = (rightSide + leftSide) / 2
            diff = (rightSide - leftSide) / 2
            f_mid = f(midpt)
        print("Iterations exceeded maximum")
        return midpt

def fixedPtIter(f, initialPt, tolerance, maxIter):
    iter = 0
    while iter < maxIter:
        p_1 = f(initialPt)
        if abs(p_1 - initialPt) < tolerance:
            return p_1
        else:
            iter += 1
            initialPt = p_1
    return "Iterations exceeded maximum"

def newtonMethod(f, f_prime, initialPt, tolerance, maxIter):
    iter = 0
    while iter < maxIter:
        p_1 = initialPt - f(initialPt) / f_prime(initialPt)
        if abs(p_1 - initialPt) < tolerance:
            return p_1
        else:
            initialPt = p_1
            iter += 1
    return "Iterations exceeded maximum"

def secantMethod(f, p_0, p_1, tolerance, maxIter):
    iter = 1
    q_0 = f(p_0)
    q_1 = f(p_1)
    while iter < maxIter:
        p_2 = p_1 - q_1 * (p_1 - p_0) / (q_1 - q_0)
        if abs(p_2 - p_1) < tolerance:
            return p_2
        else:
            p_0 = p_1
            q_0 = q_1
            p_1 = p_2
            q_1 = f(p_2)
            iter += 1
    return "Iterations exceeded maximum"

def steffensenMethod(f, initialPt, tolerance, maxIter):
    iter = 0
    while iter < maxIter:
        p_1 = f(initialPt)
        p_2 = f(p_1)
        p = initialPt - ((p_1 - initialPt) ** 2) / (p_2 - 2 * p_1 + initialPt)
        if abs(p - initialPt) < tolerance:
            return p
        else:
            iter += 1
            initialPt = p
    return p, "Iterations exceeded maximum"

def mullerMethod(f, p_0, p_1, p_2, tolerance, maxIter):
    iter = 0
    while iter < maxIter:
        h_1 = p_1 - p_0
        h_2 = p_2 - p_1
        δ_1 = (f(p_1) - f(p_2)) / h_1
        δ_2 = (f(p_2) - f(p_1)) / h_2
        d = (δ_2 - δ_1) / (h_2 + h_1)
        b = δ_2 + h_2 * d
        D = cmath.sqrt(b ** 2 - 4 * f(p_2) * d)
        if abs(b - D) < abs(b + D):
            E = b + D
        else:
            E = b - D
        h = -2 * f(p_2) / E
        p = p_2 + h
        if abs(h) < tolerance:
            return p
        else:
            p_0 = p_1
            p_1 = p_2
            p_2 = p
            h_1 = p_1 - p_0
            h_2 = p_2 - p_1
            δ_1 = (f(p_1) - f(p_0)) / h_1
            δ_2 = (f(p_2) - f(p_1)) / h_2
            d = (δ_2 - δ_1) / (h_2 + h_1)
            iter += 1
    return "Iterations exceeded maximum"

def hornerMethod(coeffArray, initialPt):
    y = coeffArray[0]
    z = coeffArray[0]
    for j in range(len(coeffArray) - 1, 0, -1):
        y = initialPt * y + coeffArray[j]
        z = initialPt * z + y
    y = initialPt * y + coeffArray[-1]
    return y, z

def lagrangeInterpolation(f, maxDeg, arrayOfPts):
    x = sympy.symbols("x")
    allLagrangeTerms = []
    for pt in arrayOfPts:
        allLagrangeTerms.append(x - pt)

    coeffPolys = []
    for i in range(maxDeg + 1):
        coeffPoly = 1
        for j in range(maxDeg + 1):
            if j == i:
                pass
            else:
                coeffPoly *= allLagrangeTerms[j] / (arrayOfPts[i] - arrayOfPts[j])
        coeffPolys.append(coeffPoly)

    finalPoly = 0
    for i in range(maxDeg + 1):
        finalPoly += coeffPolys[i] * f(arrayOfPts[i])

    return sympy.lambdify(x, finalPoly)

def nevilleInterpolation(arrayX, arrayY, f = None):
    x = sympy.symbols("x")
    n = len(arrayX)
    Qvals = np.empty((n, n), dtype = sympy.core.add.Add)
    if f == None:
        for i in range(n):
            Qvals[i, 0] = arrayY[i]
    else:
        for i in range(n):
            Qvals[i, 0] = f(arrayX[i])
    
    for i in range(1, n):
        for j in range(1, i + 1):
            Qvals[i, j] = ((x - arrayX[i - j]) * Qvals[i, j - 1]) - ((x - arrayX[i]) * Qvals[i - 1, j - 1]) / (arrayX[i] - arrayX[i - j])

    return Qvals

def newtonDividedDiff(arrayX, arrayY, f = None):
    n = len(arrayX)
    Fvals = np.empty((n, n))

    if f == None:
        for i in range(n):
            Fvals[i, 0] = arrayY[i]
    else:
        for i in range(n):
            Fvals[i, 0] = f(arrayX[i])

    for i in range(1, n):
        for j in range(1, i + 1):
            Fvals[i, j] = (Fvals[i, j - 1] - Fvals[i - 1, j - 1]) / (arrayX[i] - arrayX[i - j])

    output = []
    for i in range(n):
        output.append(Fvals[i, i])

    return output

def dividedDiffIntoPoly(arrayOfCoeffs, arrayX):
    x = sympy.symbols("x")
    P = 0

    for i in range(len(arrayOfCoeffs)):
        xterm = 1
        for j in range(i):
            xterm *= x - arrayX[j]
        P += arrayOfCoeffs[i] * xterm

    return sympy.lambdify(x, P)

def hermiteInterpolation(arrayX, arrayY, arrayYPrime, f = None, f_prime = None):
    n = len(arrayX) - 1

    arrayZ = np.empty(2 * n + 2)
    Qvals = np.empty((2 * n + 2, 2 * n + 2))
    for i in range(n + 1):
        arrayZ[2 * i] = arrayX[i]
        arrayZ[2 * i + 1] = arrayX[i]

        if f == None and f_prime == None:
            Qvals[2 * i, 0] = arrayY[i]
            Qvals[2 * i + 1, 0] = arrayY[i]
            Qvals[2 * i + 1, 1] = arrayYPrime[i]
        else:
            Qvals[2 * i, 0] = f(arrayX[i])
            Qvals[2 * i + 1, 0] = f(arrayX[i])
            Qvals[2 * i + 1, 1] = f_prime(arrayX[i])

        if i != 0:
            Qvals[2 * i, 1] = (Qvals[2 * i, 0] - Qvals[2 * i - 1, 0]) / (arrayZ[2 * i] - arrayZ[2 * i - 1])

    for i in range(2, 2 * n + 2):
        for j in range(2, i + 1):
            Qvals[i, j] = (Qvals[i, j - 1] - Qvals[i - 1, j - 1]) / (arrayZ[i] - arrayZ[i - j])

    coeffs = []
    for i in range(2 * n + 2):
        coeffs.append(Qvals[i, i])

    return coeffs

def hermiteIntoPoly(arrayOfCoeffs, arrayX):
    x = sympy.symbols("x")
    n = (len(arrayOfCoeffs) - 2) // 2

    H = 0
    xterm = 1
    for i in range(2 * n + 2):
        H += arrayOfCoeffs[i] * xterm
        xterm *= x - arrayX[i // 2]
    
    return sympy.lambdify(x, H)

def naturalCubicSpline(arrayX, arrayY):
    n = len(arrayX) - 1

    arrayh = []
    for i in range(n):
        arrayh.append(arrayX[i + 1] - arrayX[i])

    arrayα = []
    for i in range(1, n):
        arrayα.append(3 / arrayh[i] * (arrayY[i + 1] - arrayY[i]) - 3 / arrayh[i - 1] * (arrayY[i] - arrayY[i - 1]))

    l_0 = 1
    μarray = [0] * n; zarray = [0] * (n + 1)
    for i in range(1, n):
        l_0 = 2 * (arrayX[i + 1] - arrayX[i - 1]) - arrayh[i - 1] * μarray[i - 1]
        μarray[i] = arrayh[i] / l_0
        zarray[i] = (arrayα[i - 1] - arrayh[i - 1] * zarray[i - 1]) / l_0
    
    barray = [0] * n; carray = [0] * (n + 1); darray = [0] * n
    for j in range(n - 1, -1, -1):
        carray[j] = zarray[j] - μarray[j] * carray[j + 1]
        barray[j] = (arrayY[j + 1] - arrayY[j]) / arrayh[j] - arrayh[j] * (carray[j + 1] + 2 * carray[j]) / 3
        darray[j] = (carray[j + 1] - carray[j]) / (3 * arrayh[j])

    output = []
    for j in range(n):
        output.append([arrayY[j], barray[j], carray[j], darray[j]])
    
    return output

def clampedCubicSpline(arrayX, arrayY, firstLastDerivs):
    n = len(arrayX) - 1

    arrayh = []
    for i in range(n):
        arrayh.append(arrayX[i + 1] - arrayX[i])

    arrayα = [0] * (n + 1)
    arrayα[0] = 3 * (arrayY[1] - arrayY[0]) / arrayh[0] - 3 * firstLastDerivs[0]
    arrayα[n] = 3 * firstLastDerivs[1] - 3 * (arrayY[n] - arrayY[n - 1]) / arrayh[n - 1]
    for i in range(1, n):
        arrayα[i] = 3 / arrayh[i] * (arrayY[i + 1] - arrayY[i]) - 3 / arrayh[i - 1] * (arrayY[i] - arrayY[i - 1])

    l_0 = 2 * arrayh[0]; μarray = [0] * n; zarray = [0] * (n + 1)
    μarray[0] = 0.5; zarray[0] = arrayα[0] / l_0
    for i in range(1, n):
        l_0 = 2 * (arrayX[i + 1] - arrayX[i - 1]) - arrayh[i - 1] * μarray[i - 1]
        μarray[i] = arrayh[i] / l_0
        zarray[i] = (arrayα[i] - arrayh[i - 1] * zarray[i - 1]) / l_0
    l_0 = arrayh[n - 1] * (2 - μarray[n - 1])
    zarray[n] = (arrayα[n] - arrayh[n - 1] * zarray[n - 1]) / l_0

    barray = [0] * n; carray = [zarray[n]] * (n + 1); darray = [0] * n
    for j in range(n - 1, -1 , -1):
        carray[j] = zarray[j] - μarray[j] * carray[j + 1]
        barray[j] = (arrayY[j + 1] - arrayY[j]) / arrayh[j] - arrayh[j] * (carray[j + 1] + 2 * carray[j]) / 3
        darray[j] = (carray[j + 1] - carray[j]) / (3 * arrayh[j])

    output = []
    for j in range(n):
        output.append([arrayY[j], barray[j], carray[j], darray[j]])
    
    return output

def splineEval(abcdArray, index, x, x_j):
    coeffs = abcdArray[index]
    return coeffs[0] + coeffs[1] * (x - x_j) + coeffs[2] * (x - x_j) ** 2 + coeffs[3] * (x - x_j) ** 3

def bézierCurve(arrayPts, arrayLeftGuide, arrayRightGuide):
    n = len(arrayPts) - 1
    output = []

    for i in range(0, n):
        coeffs = [0] * 8
        coeffs[0] = arrayPts[i][0]
        coeffs[4] = arrayPts[i][1]
        coeffs[1] = 3 * (arrayLeftGuide[i][0] - arrayPts[i][0])
        coeffs[5] = 3 * (arrayLeftGuide[i][1] - arrayPts[i][1])
        coeffs[2] = 3 * (arrayPts[i][0] + arrayRightGuide[i][0] - 2 * arrayLeftGuide[i][0])
        coeffs[6] = 3 * (arrayPts[i][1] + arrayRightGuide[i][1] - 2 * arrayLeftGuide[i][1])
        coeffs[3] = arrayPts[i + 1][0] - arrayPts[i][0] + 3 * arrayLeftGuide[i][0] - 3 * arrayRightGuide[i][0]
        coeffs[7] = arrayPts[i + 1][1] - arrayPts[i][1] + 3 * arrayLeftGuide[i][1] - 3 * arrayRightGuide[i][1]
        output.append(coeffs)

    return output

def bézierCoeffsIntoPoly(arrayOfCoeffs, index):
    t = sympy.symbols("t")
    coeffs = arrayOfCoeffs[index]

    xt = coeffs[0] + coeffs[1] * t + coeffs[2] * t ** 2 + coeffs[3] * t ** 3
    yt = coeffs[4] + coeffs[5] * t + coeffs[6] * t ** 2 + coeffs[7] * t ** 3

    return xt, yt

def compositeTrapezoidal(f, tupleRange, n):
    h = (tupleRange[1] - tupleRange[0]) / n

    _output = f(tupleRange[0]) + f(tupleRange[1])

    x_j = tupleRange[0]
    for _ in range(1, n):
        x_j += h
        _output += 2 * f(x_j)

    return _output * h / 2

def compositeSimpson(f, tupleRange, n):
    if n % 2 == 1:
        raise ValueError("n must be even")

    h = (tupleRange[1] - tupleRange[0]) / n

    XI0 = f(tupleRange[0]) + f(tupleRange[1])
    XI1 = 0; XI2 = 0

    x_j = tupleRange[0]
    for i in range(1, n):
        x_j += h
        if i % 2 == 0:
            XI2 += f(x_j)
        else:
            XI1 += f(x_j)
    
    XI = h * (XI0 + 2 * XI2 + 4 * XI1) / 3
    return XI

def rombergIntegration(f, tupleRange, n):
    output = []
    arrayR = np.empty((2, n))
    
    h = tupleRange[1] - tupleRange[0]

    arrayR[0, 0] = h / 2 * (f(tupleRange[0]) + f(tupleRange[1]))
    output.append(arrayR[0, 0])

    for i in range(2, n + 1):
        trapezoidSum = 0
        for k in range(1, 2 ** (i - 2) + 1):
            trapezoidSum += f(tupleRange[0] + (k - 0.5) * h)
        arrayR[1, 0] = 1 / 2 * (arrayR[0, 0] + h * trapezoidSum)
        
        for j in range(2, i + 1):
            arrayR[1, j - 1] = arrayR[1, j - 2] + (arrayR[1, j - 2] - arrayR[0, j - 2]) / (4 ** (j - 1) - 1)

        for j in range(1, i + 1):
            output.append(arrayR[1, j - 1])

        h = h / 2

        for j in range(1, i + 1):
            arrayR[0, j - 1] = arrayR[1, j - 1]
    
    return output

def simpsonDoubleIntegral(f, tupleRange, c, d, m, n):
    b = tupleRange[1]; a = tupleRange[0]
    h = (b - a) / n

    arrayJ = [0, 0, 0]; arrayK = [0, 0, 0]

    for i in range(0, n + 1):
        x = a + i * h
        HX = (d(x) - c(x)) / m
        arrayK[0] = f(x, c(x)) + f(x, d(x))
        arrayK[1] = 0
        arrayK[2] = 0

        for j in range(1, m):
            y = c(x) + j * HX
            Q = f(x, y)

            if j % 2 == 0:
                arrayK[1] += Q
            else:
                arrayK[2] += Q

        L = (arrayK[0] + 2 * arrayK[1] + 4 * arrayK[2]) * HX / 3

        if i == 0 or i == n:
            arrayJ[0] += L
        elif i % 2 == 0:
            arrayJ[1] += L
        else:
            arrayJ[2] += L

    J = h * (arrayJ[0] + 2 * arrayJ[1] + 4 * arrayJ[2]) / 3

    return J

def eulerMethod(f, tupleRange, N, initialCond):
    a = tupleRange[0]; b = tupleRange[1]
    h = (b - a) / N
    t = a
    w = initialCond

    output = [[t, w]]

    for _ in range(1, N + 1):
        w += h * f(t, w)
        t += h
        output.append([t, w])

    return output

def modEulerMethod(f, tupleRange, N, initialCond):
    output = [[tupleRange[0], initialCond]]
    h = (tupleRange[1] - tupleRange[0]) / N

    t_i = tupleRange[0]
    for _ in range(N):
        w_i = output[-1][1]
        val = w_i + h / 2 * (f(t_i, w_i) + f(t_i + h, w_i + h * f(t_i, w_i)))
        t_i += h

        output.append([t_i, val])

    return output

def rungeKuttaFehlberg(f, tupleRange, initialCond, tol, hmax, hmin):
    t = tupleRange[0]; b = tupleRange[1]
    w = initialCond
    h = hmax
    flag = 1

    output = [[t, w, h]]

    while flag == 1:
        K_1 = h * f(t, w)
        K_2 = h * f(t + h / 4, w + K_1 / 4)
        K_3 = h * f(t + 3 * h / 8, w + 3 * K_1 / 32 + 9 * K_2 / 32)
        K_4 = h * f(t + 12 * h / 13, w + 1932 * K_1 / 2197 - 7200 * K_2 / 2197 + 7296 * K_3 / 2197)
        K_5 = h * f(t + h, w + 439 * K_1 / 216 - 8 * K_2 + 3680 * K_3 / 513 - 845 * K_4 / 4104)
        K_6 = h * f(t + h / 2, w - 8 * K_1 / 27 + 2 * K_2 - 3544 * K_3 / 2565 + 1859 * K_4 / 4104 - 11 * K_5 / 40)

        R = abs(K_1 / 360 - 128 * K_3 / 4275 - 2197 * K_4 / 75240 + K_5 / 50 + 2 * K_6 / 55) / h

        if R <= tol:
            t += h
            w += 25 * K_1 / 216 + 1408 * K_3 / 2565 + 2197 * K_4 / 4104 - K_5 / 5
            output.append([t, w, h])
        
        δ = 0.84 * (tol / R) ** 0.25

        if δ <= 0.1:
            h *= 0.1
        elif δ >= 4:
            h *= 4
        else:
            h *= δ

        if h > hmax:
            h = hmax

        if t >= b:
            flag = 0
        elif t + h > b:
            h = b - t
        elif h < hmin:
            flag = 0
            ValueError("min h exceeded")
    
    return output

def rungeKutta4(f, tupleRange, N, initialCond):
    a = tupleRange[0]; b = tupleRange[1]
    h = (b - a) / N
    t = a
    w = initialCond

    output = [[t, w]]

    for _ in range(1, N + 1):
        K_1 = h * f(t, w)
        K_2 = h * f(t + h / 2, w + K_1 / 2)
        K_3 = h * f(t + h / 2, w + K_2 / 2)
        K_4 = h * f(t + h, w + K_3)

        w += (K_1 + 2 * K_2 + 2 * K_3 + K_4) / 6
        t += h

        output.append([t, w])

    return output

def adamsBashforth2(f, tupleRange, N, initialConds):
    w_0 = initialConds[0]
    w_1 = initialConds[1]

    t_0 = tupleRange[0]
    h = (tupleRange[1] - tupleRange[0]) / N
    t_1 = t_0 + h

    outputs = initialConds.copy()

    for _ in range(1, N):
        w_0 = outputs[-2]
        w_1 = outputs[-1]
        approx = w_1 + h * (3 * f(t_1, w_1) - f(t_0, w_0)) / 2
        outputs.append(approx)
        t_0 = t_1
        t_1 += h

    return outputs

def adamsBashforth3(f, tupleRange, N, initialConds):
    w_0 = initialConds[0]
    w_1 = initialConds[1]
    w_2 = initialConds[2]

    t_0 = tupleRange[0]
    h = h = (tupleRange[1] - tupleRange[0]) / N
    t_1 = t_0 + h
    t_2 = t_1 + h

    outputs = initialConds.copy()

    for _ in range(2, N):
        w_0 = outputs[-3]
        w_1 = outputs[-2]
        w_2 = outputs[-1]
        approx = w_2 + h * (23 * f(t_2, w_2) - 16 * f(t_1, w_1) + 5 * f(t_0, w_0)) / 12
        outputs.append(approx)
        t_0 = t_1
        t_1 = t_2
        t_2 += h

    return outputs

def adamsBashforth4(f, tupleRange, N, initialConds):
    w_0 = initialConds[0]; w_1 = initialConds[1]; w_2 = initialConds[2]; w_3 = initialConds[3]

    t_0 = tupleRange[0]
    h = (tupleRange[1] - tupleRange[0]) / N
    t_1 = t_0 + h; t_2 = t_1 + h; t_3 = t_2 + h

    outputs = initialConds.copy()

    for _ in range(3, N):
        w_0 = outputs[-4]; w_1 = outputs[-3]; w_2 = outputs[-2]; w_3 = outputs[-1]
        approx = w_3 + h * (55 * f(t_3, w_3) - 59 * f(t_2, w_2) + 37 * f(t_1, w_1) - 9 * f(t_0, w_0)) / 24
        outputs.append(approx)
        t_0 = t_1; t_1 = t_2; t_2 = t_3; t_3 += h

    return outputs

def adamsBashforth5(f, tupleRange, N, initialConds):
    w_0 = initialConds[0]; w_1 = initialConds[1]; w_2 = initialConds[2]; w_3 = initialConds[3]; w_4 = initialConds[4]

    t_0 = tupleRange[0]
    h = (tupleRange[1] - tupleRange[0]) / N
    t_1 = t_0 + h; t_2 = t_1 + h; t_3 = t_2 + h; t_4 = t_3 + h

    outputs = initialConds.copy()

    for _ in range(4, N):
        w_0 = outputs[-5]; w_1 = outputs[-4]; w_2 = outputs[-3]; w_3 = outputs[-2]; w_4 = outputs[-1]
        approx = w_4 + h * (1901 * f(t_4, w_4) - 2774 * f(t_3, w_3) + 2616 * f(t_2, w_2) - 1274 * f(t_1, w_1) + 251 * f(t_0, w_0)) / 720
        outputs.append(approx)
        t_0 = t_1; t_1 = t_2; t_2 = t_3; t_3 = t_4; t_4 += h
    
    return outputs

def adamsVarStepSizeP_C(f, tupleRange, initialCond, tol, hmax, hmin):
    def RK4(h, v_0, x_0, v_1, x_1, v_2, x_2, v_3, x_3):
        xs = [x_0, x_1, x_2, x_3]
        vs = [v_0, v_1, v_2, v_3]

        output = []

        for j in range(1, 4):
            K_1 = h * f(xs[j - 1], vs[j - 1])
            K_2 = h * f(xs[j - 1] + h / 2, vs[j - 1] + K_1 / 2)
            K_3 = h * f(xs[j - 1] + h / 2, vs[j - 1] + K_2 / 2)
            K_4 = h * f(xs[j - 1] + h, vs[j - 1] + K_3)

            vs[j] = vs[j - 1] + (K_1 + 2 * K_2 + 2 * K_3 + K_4) / 6
            xs[j] = xs[0] + j * h

            output.append([xs[j], vs[j]])

        return output

    a = tupleRange[0]; b = tupleRange[1]
    t_0 = a; w_0 = initialCond; h = hmax; flag = 1; last = 0

    outputs = [[t_0, w_0]]

    rk4out = RK4(h, w_0, t_0, 0, 0, 0, 0, 0, 0)
    nflag = 1; i = 4; t = rk4out[2][0] + h
    t_1 = rk4out[0][0]; t_2 = rk4out[1][0]; t_3 = rk4out[2][0]
    w_1 = rk4out[0][1]; w_2 = rk4out[1][1]; w_3 = rk4out[2][1]
    tvals = [t_0, t_1, t_2, t_3]; wvals = [w_0, w_1, w_2, w_3]

    while flag == 1:
        WP = wvals[i - 1] + h * (55 * f(tvals[i - 1], wvals[i - 1]) - 59 * f(tvals[i - 2], wvals[i - 2]) + 37 * f(tvals[i - 3], wvals[i - 3]) - 9 * f(tvals[i - 4], wvals[i - 4])) / 24
        WC = wvals[i - 1] + h * (9 * f(t, WP) + 19 * f(tvals[i - 1], wvals[i - 1]) - 5 * f(tvals[i - 2], wvals[i - 2]) + f(tvals[i - 3], wvals[i - 3])) / 24

        σ = 19 * abs(WC - WP) / (270 * h)

        if σ <= tol:
            wvals.append(WC); tvals.append(t)
            if nflag == 1:
                for j in range(i - 3, i - 1, -1):
                    outputs.append([j, tvals[j], wvals[j], h])
            else:
                outputs.append([i, tvals[i], wvals[i], h])

            if last == 1:
                flag = 0
            else:
                i += 1
                nflag = 0

                if σ <= 0.1 * tol or tvals[i - 1] + h > b:
                    q = (tol / (2 * σ)) ** 0.25
                    if q > 4:
                        h *= 4
                    else:
                        h *= q
                    
                    if h > hmax:
                        h = hmax
                    
                    if tvals[i - 1] + 4 * h > b:
                        h = (b - tvals[i - 1]) / 4
                        last = 1
                    
                    rk4out = RK4(h, wvals[i - 1], tvals[i - 1], 0, 0, 0, 0, 0, 0)
                    
                    tvals.append(rk4out[0][0]); tvals.append(rk4out[1][0]); tvals.append(rk4out[2][0])
                    wvals.append(rk4out[0][1]); wvals.append(rk4out[1][1]); wvals.append(rk4out[2][1])

                    nflag = 1; i += 3
        else:
            q = (tol / (2 * σ)) ** 0.25
            
            if q < 0.1:
                h *= 0.1
            else:
                h *= q
            
            if h < hmin:
                flag = 0
                raise ValueError("minimum h exceeded")
            else:
                if nflag == 1:
                    i -= 3
                    rk4out = RK4(h, wvals[i - 1], tvals[i - 1], wvals[i], tvals[i], wvals[i + 1], tvals[i + 1], wvals[i + 2], tvals[i + 2])
                    
                    tvals[i] = rk4out[0][0]; tvals[i + 1] = rk4out[1][0]; tvals[i + 2] = rk4out[2][0]
                    wvals[i] = rk4out[0][1]; wvals[i + 1] = rk4out[1][1]; wvals[i + 2] = rk4out[2][1]

                    i += 3
                    nflag = 1

        t = tvals[i - 1] + h

    return outputs

def rungeKuttaSystem(fArr, tupleRange, N, initialCondArr):
    a = tupleRange[0]; b = tupleRange[1]
    h = (b - a) / N; m = len(fArr)

    t = a
    outputs = []

    kArr = np.empty((4, m))

    wArr = initialCondArr.copy()
    output = [t]; output.extend(wArr)
    outputs.append(output.copy())

    for i in range(1, N + 1):
        for j in range(0, m):
            fInputs = output
            kArr[0, j] = h * fArr[j](fInputs)

        fInputs = np.empty(m + 1)
        fInputs[0] = t + h / 2
        for k in range(1, m + 1):
            fInputs[k] = output[k] + kArr[0, k - 1] / 2
        for j in range(0, m):
            kArr[1, j] = h * fArr[j](fInputs)
        for k in range(1, m + 1):
                fInputs[k] = output[k] + kArr[1, k - 1] / 2
        for j in range(0, m):
            kArr[2, j] = h * fArr[j](fInputs)
        fInputs[0] = t + h
        for k in range(1, m + 1):
            fInputs[k] = output[k] + kArr[2, k - 1]
        for j in range(0, m):
            kArr[3, j] = h * fArr[j](fInputs)
        
        for j in range(1, m + 1):
            output[j] = output[j] + (kArr[0, j - 1] + 2 * kArr[1, j - 1] + 2 * kArr[2, j - 1] + kArr[3, j - 1]) / 6

        t += h
        output[0] = t

        outputs.append(output.copy())

    return outputs

def adams4PC(f, tupleRange, N, initialCond):
    a = tupleRange[0]; b = tupleRange[1]
    h = (b - a) / N

    t0 = a; w0 = initialCond

    outputs = [[t0, w0]]

    for _ in range(1, 4):
        K_1 = h * f(t0, w0)
        K_2 = h * f(t0 + h / 2, w0 + K_1 / 2)
        K_3 = h * f(t0 + h / 2, w0 + K_2 / 2)
        K_4 = h * f(t0 + h, w0 + K_3)

        w0 += (K_1 + 2 * K_2 + 2 * K_3 + K_4) / 6
        t0 += h

        outputs.append([t0, w0])

    w3 = outputs[3][1]; t3 = outputs[3][0]; w2 = outputs[2][1]; t2 = outputs[2][0]; w1 = outputs[1][1]; t1 = outputs[1][0]; w0 = outputs[0][1]; t0 = outputs[0][0]
    for i in range(4, N + 1):
        t = a + i * h

        w = w3 + h * (55 * f(t3, w3) - 59 * f(t2, w2) + 37 * f(t1, w1) - 9 * f(t0, w0)) / 24
        w = w3 + h * (9 * f(t, w) + 19 * f(t3, w3) - 5 * f(t2, w2) + f(t1, w1)) / 24

        outputs.append([t, w])

        t0 = t1; t1 = t2; t2 = t3; w0 = w1; w1 = w2; w2 = w3; t3 = t; w3 = w

    return outputs    

def trapezoidalWithNewton(f, f_y, tupleRange, N, initialCond, tol, maxIter):
    a = tupleRange[0]; b = tupleRange[1]
    h = (b - a) / N
    t = a; w = initialCond

    outputs = [[t, w]]

    for _ in range(1, N + 1):
        k1 = w + h * f(t, w) / 2
        w0 = k1; j = 1; flag = 0

        while flag == 0:
            w = w0 - (w0 - h / 2 * f(t + h, w0) - k1) / (1 - h / 2 * f_y(t + h, w0))
            if abs(w - w0) < tol:
                flag = 1
            else:
                j += 1; w0 = w
                if j > maxIter:
                    raise ValueError("Exceeded max iterations")

        t += h

        outputs.append([t, w])
    
    return outputs

def backwardEuler(f, f_y, tupleRange, N, initialCond, tol, maxIter):
    a = tupleRange[0]; b = tupleRange[1]
    h = (b - a) / N
    t = a; w = initialCond

    outputs = [[t, w]]

    for _ in range(1, N + 1):
        k1 = w + h * f(t, w)
        w0 = k1; j = 1; flag = 0

        while flag == 0:
            w = w0 - (w0 - w - h * f(t + h, w0)) / (1 - h * f_y(t + h, w0))
            if abs(w - w0) < tol:
                flag = 1
            else:
                j += 1; w0 = w
                if j > maxIter:
                    raise ValueError("Exceeded max iterations")

        t += h

        outputs.append([t, w])
    
    return outputs