import math

class params:
    def __init__(self, root_tol, func_tol, maxit):
        self.root_tol = root_tol
        self.func_tol = func_tol
        self.maxit = maxit

class info:
    def __init__(self, flag):
        self.flag = flag

def bisectionMethod(f, tupleRange, tol, maxIter): #copied from homework/textbook
    leftSide = tupleRange[0]; rightSide = tupleRange[1]
    f_left = f(leftSide)
    f_right = f(rightSide)
    if abs(f_left) <= tol: #unused
        return leftSide
    elif abs(f_right) <= tol: #unused
        return rightSide
    else:
        midpt = (rightSide + leftSide) / 2
        diff = (rightSide - leftSide) / 2
        f_mid = f(midpt)
        iter = 1
        while f_mid != 0 and diff > tol and iter < maxIter:
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
        return leftSide, rightSide, midpt

def IQI(f, x0, x1, x2):
    y0 = f(x0); y1 = f(x1); y2 = f(x2)
    approx = (y1 * y2 * x0) / ((y0 - y1) * (y0 - y2)) + (y0 * y2 * x1) / ((y1 - y0) * (y1 - y2)) + (y0 * y1 * x2) / ((y2 - y0) * (y2 - y1))
    return approx

def modifiedZeroin(f, initialInterval, params):
    func_tol = params.func_tol; root_tol = params.root_tol
    a = initialInterval[0]; b = initialInterval[1]
    iter = 0; iterpIter = 0

    x0 = a; x1 = b; x2 = (a + b) / 2

    while iter <= params.maxit:
        if x1 - x0 < root_tol:
            if abs(f(x2)) < func_tol:
                return [x2, info(0)]
            else:
                return [x2, info(1)]

        x3 = IQI(f, x0, x1, x2); iterpIter += 1
        try:
            if abs(f(x3)) < func_tol:
                return [x3, info(0)]
        except ValueError:
            pass
        
        if x3 < x0 or x3 > x1:
            x0, x1, x2 = bisectionMethod(f, (x0, x1), func_tol, 3)
            if iterpIter == 1:
                f3 = f(x2)

            #iterpIter = 0
        else:
            if x3 < x2:
                x0 = x0; x1 = x2; x2 = x3
            else:
                x0 = x2; x1 = x1; x2 = x3

            if iterpIter == 1:
                f3 = f(x3)

        if iterpIter > 5 and abs(f(x3)) > abs(f3 / 2):
            x0, x1, x2 = bisectionMethod(f, (x0, x1), func_tol, 3)
            #iterpIter = 0

        iter += 1

    if abs(f(x3)) < func_tol:
        return [x3, info(0)]
    else:
        return [x3, info(1)]

###########################
parameters = params(1e-9, 1e-9, 100)
def f1(x):
    return math.sqrt(x) - math.cos(x)
def f2(x):
    return 3 * (x + 1) * (x - 0.5) * (x - 1)
def f3(x):
    return math.log(x - 1)
def f4(x):
    return math.pi + 5 * math.sin(x / 2) - x
def f5(x):
    return (x - 4) ** 7 * (x - 2) ** 11

ans1 = modifiedZeroin(f1, [0, 1], parameters)
print(ans1[0], "flag:", ans1[1].flag)
print("close to 0?", f1(ans1[0]))

ans2 = modifiedZeroin(f2, [-2, 0.2], parameters)
print(ans2[0], "flag:", ans2[1].flag)
print("close to 0?", f2(ans2[0]))

ans3 = modifiedZeroin(f3, [1.5, 1e3], parameters)
print(ans3[0], "flag:", ans3[1].flag)
print("close to 0?", f3(ans3[0]))

ans4 = modifiedZeroin(f4, [0, 6.3], parameters)
print(ans4[0], "flag:", ans4[1].flag)
print("close to 0?", f4(ans4[0]))

ans5 = modifiedZeroin(f5, [1.8, 3.4], parameters)
print(ans5[0], "flag:", ans5[1].flag)
print("close to 0?", f5(ans5[0]))