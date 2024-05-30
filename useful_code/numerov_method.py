# method solves d^2X(t)/dx^2 = f(t)X(t) for boundary conditions X(t_1) = X_1, X(t_2) = X_2
# to h^6

def w(x,t,h,f):
    return (1 - f(t)*(h**2)/12)*x

def wInv(w,t,h,f):
    return w/(1 - f(t)*(h**2)/12)

def iterateNumerov(x_t,x_tMinush,t,h,f):
    w_t = w(x_t,t,h,f)
    w_tMinush = w(x_tMinush,t,h,f)
    w_rPlush = 2*w_t - w_tMinush + x_t*f(t)*(h**2)
    return wInv(w_rPlush,t,h,f)

def runNumerov(t_0,h,x_0,x_1,f,N):
    xList = [x_0,x_1]
    tList = [t_0,t_0 + h]
    for i in range(N):
        x_i = xList[-1]
        x_iMinus1 = xList[-2]
        t = tList[-1]
        xList.append(iterateNumerov(x_i,x_iMinus1,t,h,f))
        tList.append(t+h)
    return tList, xList