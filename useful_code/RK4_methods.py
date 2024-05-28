

#should update to adjust h automatically, see book, create new function
def iterateRK4(t_n,y_n,h,f):
    """
    Iterates by RK4 method.
    """
    k_1 = f(t_n, y_n)
    k_2 = f(t_n + h/2, y_n + h*k_1/2)
    k_3 = f(t_n + h/2, y_n + h*k_2/2)
    k_4 = f(t_n + h, y_n + h*k_3)

    return y_n + h*(k_1 + 2*k_2 + 2*k_3 + k_4)/6


def simulate(t_0,y_0,numberOfRuns,h,f):
    """
    Uses RK4 method to find a set of points on the solution to dY/dt=f(t,Y) given an intial condition.
    """
    yList = [y_0]
    tList = [t_0]

    for x in range(numberOfRuns):
        yList.append(iterateRK4(tList[-1],yList[-1],f))
        tList.append(tList[-1]+h)
    return tList, yList

def transformToXYZ(lst,N):
    """
    Transfroms [np.array([1,2]),np.array([3,4])] to [[1,3],[2,4]]
    """
    return [[s[c] for s in lst] for c in range(N)]

