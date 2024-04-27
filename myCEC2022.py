import numpy as np
import math
INF = 1.0e99
EPS = 1.0e-14
E  = 2.7182818284590452353602874713526625
PI = 3.1415926535897932384626433832795029
ini_flag = 0


def ellips_func(x, f, nx, Os, Mr, s_flag, r_flag):  # Ellipsoidal
    z = sr_func(x, nx, Os, Mr, 1.0, s_flag, r_flag)  # shift and rotate
    for i in range(nx):
        f[0] += (10.0 ** (6.0 * i / (nx - 1))) * z[i] ** 2

def bent_cigar_func(x, f, nx, Os, Mr, s_flag, r_flag):  # Bent_Cigar
    z = sr_func(x, nx, Os, Mr, 1.0, s_flag, r_flag)  # shift and rotate
    f[0] = z[0] ** 2
    for i in range(1, nx):
        f[0] += (10.0 ** 6.0) * z[i] ** 2

def discus_func(x, f, nx, Os, Mr, s_flag, r_flag):  # Discus
    z = sr_func(x, nx, Os, Mr, 1.0, s_flag, r_flag)  # shift and rotate
    f[0] = (10.0 ** 6.0) * z[0] ** 2
    for i in range(1, nx):
        f[0] += z[i] ** 2

def rosenbrock_func(x, f, nx, Os, Mr, s_flag, r_flag):  # Rosenbrock's
    z = sr_func(x, nx, Os, Mr, 2.048 / 100.0, s_flag, r_flag)  # shift and rotate
    z = [zi + 1.0 for zi in z]  # shift to origin
    for i in range(nx - 1):
        tmp1 = z[i] ** 2 - z[i + 1]
        tmp2 = z[i] - 1.0
        f[0] += 100.0 * tmp1 ** 2 + tmp2 ** 2

def ackley_func(x, f, nx, Os, Mr, s_flag, r_flag):  # Ackley's
    sum1, sum2 = 0.0, 0.0
    z = sr_func(x, nx, Os, Mr, 1.0, s_flag, r_flag)  # shift and rotate
    for i in range(nx):
        sum1 += z[i] ** 2
        sum2 += cos(2.0 * pi * z[i])
    sum1 = -0.2 * sqrt(sum1 / nx)
    sum2 /= nx
    f[0] = E - 20.0 * exp(sum1) - exp(sum2) + 20.0

def griewank_func(x, f, nx, Os, Mr, s_flag, r_flag):
    s = 0.0
    p = 1.0

    sr_func(x, z, nx, Os, Mr, 600.0/100.0, s_flag, r_flag)

    for i in range(nx):
        s += z[i]*z[i]
        p *= math.cos(z[i]/math.sqrt(1.0+i))
    f[0] = 1.0 + s/4000.0 - p


def levy_func(x, f, nx, Os, Mr, s_flag, r_flag):
    f[0] = 0.0
    sr_func(x, z, nx, Os, Mr,1.0, s_flag, r_flag)

    w = [1.0 + (z[i] - 0.0)/4.0 for i in range(nx)]

    term1 = pow((sin(pi*w[0])),2)
    term3 = pow((w[nx-1]-1),2) * (1+pow((sin(2*pi*w[nx-1])),2))

    sum = 0.0

    for i in range(nx-1):
        wi = w[i]
        newv = pow((wi-1),2) * (1+10*pow((sin(pi*wi+1)),2))
        sum = sum + newv

    f[0] = term1 + sum + term3

def schaffer_F7_func(x, f, nx, Os, Mr, s_flag, r_flag):
    sr_func(x, z, nx, Os, Mr, 1.0, s_flag, r_flag)  # shift and rotate

    f[0] = 0.0
    for i in range(nx-1):
        z[i] = math.sqrt(z[i]*z[i] + z[i+1]*z[i+1])
        tmp = math.sin(50.0 * math.pow(z[i], 0.2))
        f[0] += math.pow(z[i], 0.5) + math.pow(z[i], 0.5) * tmp * tmp

    f[0] = f[0] * f[0] / (nx-1) / (nx-1)

def step_rastrigin_func(x, f, nx, Os, Mr, s_flag, r_flag):
    f[0] = 0.0

    for i in range(nx):
        if abs(y[i] - Os[i]) > 0.5:
            y[i] = Os[i] + math.floor(2 * (y[i] - Os[i]) + 0.5) / 2

    sr_func(x, z, nx, Os, Mr, 5.12/100.0, s_flag, r_flag)  # shift and rotate

    for i in range(nx):
        f[0] += (z[i]*z[i] - 10.0 * math.cos(2.0 * math.pi * z[i]) + 10.0)


def katsuura_func(x, f, nx, Os, Mr, s_flag, r_flag):
    f[0] = 1.0
    tmp3 = pow(1.0*nx, 1.2)

    sr_func(x, z, nx, Os, Mr, 5.0/100.0, s_flag, r_flag)  # shift and rotate

    for i in range(nx):
        temp = 0.0
        for j in range(1, 33):
            tmp1 = pow(2.0, j)
            tmp2 = tmp1 * z[i]
            temp += abs(tmp2 - math.floor(tmp2 + 0.5)) / tmp1
        f[0] *= pow(1.0 + (i+1) * temp, 10.0 / tmp3)

    tmp1 = 10.0 / nx / nx
    f[0] = f[0] * tmp1 - tmp1

def zakharov_func(x, f, nx, Os, Mr, s_flag, r_flag):
    sr_func(x, z, nx, Os, Mr, 1.0, s_flag, r_flag)  # shift and rotate

    f[0] = 0.0
    sum1 = 0.0
    sum2 = 0.0
    for i in range(nx):
        xi = z[i]
        sum1 += math.pow(xi, 2)
        sum2 += 0.5 * (i + 1) * xi

    f[0] = sum1 + math.pow(sum2, 2) + math.pow(sum2, 4)


def hf02(x, f, nx, Os, Mr, S, s_flag, r_flag):  # Hybrid Function 2
    cf_num = 3
    fit = [0.0] * cf_num
    G = [0] * cf_num
    G_nx = [0] * cf_num
    Gp = [0.4, 0.4, 0.2]

    tmp = 0
    for i in range(cf_num - 1):
        G_nx[i] = ceil(Gp[i] * nx)
        tmp += G_nx[i]
    G_nx[cf_num - 1] = nx - tmp

    for i in range(1, cf_num):
        G[i] = G[i - 1] + G_nx[i - 1]

    sr_func(x, z, nx, Os, Mr, 1.0, s_flag, r_flag)  # shift and rotate

    y = [0.0] * nx
    for i in range(nx):
        y[i] = z[S[i] - 1]

    bent_cigar_func(y[G[0]:G[0]+G_nx[0]], fit[0], G_nx[0], Os, Mr, 0, 0)
    hgbat_func(y[G[1]:G[1]+G_nx[1]], fit[1], G_nx[1], Os, Mr, 0, 0)
    rastrigin_func(y[G[2]:G[2]+G_nx[2]], fit[2], G_nx[2], Os, Mr, 0, 0)

    f[0] = sum(fit)

def hf10(x, f, nx, Os, Mr, S, s_flag, r_flag):
    cf_num = 6
    fit = [0.0] * 6
    G = [0] * 6
    G_nx = [0] * 6
    Gp = [0.1, 0.2, 0.2, 0.2, 0.1, 0.2]

    tmp = 0
    for i in range(cf_num-1):
        G_nx[i] = math.ceil(Gp[i] * nx)
        tmp += G_nx[i]
    G_nx[cf_num-1] = nx - tmp

    G[0] = 0
    for i in range(1, cf_num):
        G[i] = G[i-1] + G_nx[i-1]

    sr_func(x, z, nx, Os, Mr, 1.0, s_flag, r_flag)  # shift and rotate

    for i in range(nx):
        y[i] = z[S[i] - 1]

    i = 0
    hgbat_func(y[G[i]:], fit[i], G_nx[i], Os, Mr, 0, 0)
    i = 1
    katsuura_func(y[G[i]:], fit[i], G_nx[i], Os, Mr, 0, 0)
    i = 2
    ackley_func(y[G[i]:], fit[i], G_nx[i], Os, Mr, 0, 0)
    i = 3
    rastrigin_func(y[G[i]:], fit[i], G_nx[i], Os, Mr, 0, 0)
    i = 4
    schwefel_func(y[G[i]:], fit[i], G_nx[i], Os, Mr, 0, 0)
    i = 5
    schaffer_F7_func(y[G[i]:], fit[i], G_nx[i], Os, Mr, 0, 0)

    f[0] = 0.0
    for i in range(cf_num):
        f[0] += fit[i]

def rastrigin_func(x, f, nx, Os, Mr, s_flag, r_flag):
    f[0] = 0.0

    sr_func(x, z, nx, Os, Mr, 5.12/100.0, s_flag, r_flag)  # shift and rotate

    for i in range(nx):
        f[0] += (z[i]*z[i] - 10.0*math.cos(2.0*math.pi*z[i]) + 10.0)



def happycat_func(x, f, nx, Os, Mr, s_flag, r_flag):
    alpha = 1.0 / 8.0

    sr_func(x, z, nx, Os, Mr, 5.0/100.0, s_flag, r_flag)  # shift and rotate

    r2 = 0.0
    sum_z = 0.0
    for i in range(nx):
        z[i] = z[i] - 1.0  # shift to origin
        r2 += z[i] * z[i]
        sum_z += z[i]

    f[0] = pow(abs(r2 - nx), 2 * alpha) + (0.5 * r2 + sum_z) / nx + 0.5
    

def schwefel_func(x, f, nx, Os, Mr, s_flag, r_flag):
    f[0] = 0.0

    sr_func(x, z, nx, Os, Mr, 1000.0/100.0, s_flag, r_flag)  # shift and rotate

    for i in range(nx):
        z[i] += 4.209687462275036e+002
        if z[i] > 500:
            f[0] -= (500.0 - (z[i] % 500)) * math.sin(math.sqrt(500.0 - (z[i] % 500)))
            tmp = (z[i] - 500.0) / 100
            f[0] += tmp * tmp / nx
        elif z[i] < -500:
            f[0] -= (-500.0 + (abs(z[i]) % 500)) * math.sin(math.sqrt(500.0 - (abs(z[i]) % 500)))
            tmp = (z[i] + 500.0) / 100
            f[0] += tmp * tmp / nx
        else:
            f[0] -= z[i] * math.sin(math.sqrt(abs(z[i])))
    f[0] += 4.189828872724338e+002 * nx

def grie_rosen_func(x, f, nx, Os, Mr, s_flag, r_flag):
    f[0] = 0.0

    sr_func(x, z, nx, Os, Mr, 5.0/100.0, s_flag, r_flag)  # shift and rotate

    z[0] += 1.0  # shift to origin
    for i in range(nx-1):
        z[i+1] += 1.0  # shift to origin
        tmp1 = z[i]*z[i] - z[i+1]
        tmp2 = z[i] - 1.0
        temp = 100.0 * tmp1 * tmp1 + tmp2 * tmp2
        f[0] += (temp * temp) / 4000.0 - math.cos(temp) + 1.0

    tmp1 = z[nx-1]*z[nx-1] - z[0]
    tmp2 = z[nx-1] - 1.0
    temp = 100.0 * tmp1 * tmp1 + tmp2 * tmp2
    f[0] += (temp * temp) / 4000.0 - math.cos(temp) + 1.0


def escaffer6_func(x, f, nx, Os, Mr, s_flag, r_flag):
    sr_func(x, z, nx, Os, Mr, 1.0, s_flag, r_flag)  # shift and rotate

    f[0] = 0.0
    for i in range(nx-1):
        temp1 = math.sin(math.sqrt(z[i]*z[i] + z[i+1]*z[i+1]))
        temp1 = temp1 * temp1
        temp2 = 1.0 + 0.001 * (z[i]*z[i] + z[i+1]*z[i+1])
        f[0] += 0.5 + (temp1 - 0.5) / (temp2 * temp2)

    temp1 = math.sin(math.sqrt(z[nx-1]*z[nx-1] + z[0]*z[0]))
    temp1 = temp1 * temp1
    temp2 = 1.0 + 0.001 * (z[nx-1]*z[nx-1] + z[0]*z[0])
    f[0] += 0.5 + (temp1 - 0.5) / (temp2 * temp2)

def hgbat_func(x, f, nx, Os, Mr, s_flag, r_flag):
    alpha = 1.0 / 4.0

    sr_func(x, z, nx, Os, Mr, 5.0/100.0, s_flag, r_flag)  # shift and rotate

    r2 = 0.0
    sum_z = 0.0
    for i in range(nx):
        z[i] = z[i] - 1.0  # shift to origin
        r2 += z[i] * z[i]
        sum_z += z[i]

    f[0] = math.pow(math.fabs(math.pow(r2, 2.0) - math.pow(sum_z, 2.0)), 2 * alpha) + (0.5 * r2 + sum_z) / nx + 0.5


def cf01(x, f, nx, Os, Mr, r_flag):
    cf_num = 5
    fit = [0.0] * cf_num
    delta = [10, 20, 30, 40, 50]
    bias = [0, 200, 300, 100, 400]

    rosenbrock_func(x, fit[0], nx, Os[0*nx], Mr[0*nx*nx], 1, r_flag)
    fit[0] = 10000 * fit[0] / 1e+4

    ellips_func(x, fit[1], nx, Os[1*nx], Mr[1*nx*nx], 1, r_flag)
    fit[1] = 10000 * fit[1] / 1e+10

    bent_cigar_func(x, fit[2], nx, Os[2*nx], Mr[2*nx*nx], 1, r_flag)
    fit[2] = 10000 * fit[2] / 1e+30

    discus_func(x, fit[3], nx, Os[3*nx], Mr[3*nx*nx], 1, r_flag)
    fit[3] = 10000 * fit[3] / 1e+10

    ellips_func(x, fit[4], nx, Os[4*nx], Mr[4*nx*nx], 1, 0)
    fit[4] = 10000 * fit[4] / 1e+10

    cf_cal(x, f, nx, Os, delta, bias, fit, cf_num)

def cf02(x, f, nx, Os, Mr, r_flag):
    cf_num = 3
    fit = [0.0] * cf_num
    delta = [20, 10, 10]
    bias = [0, 200, 100]

    schwefel_func(x, fit[0], nx, Os[0*nx], Mr[0*nx*nx], 1, 0)
    rastrigin_func(x, fit[1], nx, Os[1*nx], Mr[1*nx*nx], 1, r_flag)
    hgbat_func(x, fit[2], nx, Os[2*nx], Mr[2*nx*nx], 1, r_flag)

    cf_cal(x, f, nx, Os, delta, bias, fit, cf_num)

def cf06(x, f, nx, Os, Mr, r_flag):
    cf_num = 5
    fit = [0.0] * cf_num
    delta = [20, 20, 30, 30, 20]
    bias = [0, 200, 300, 400, 200]

    escaffer6_func(x, fit[0], nx, Os[0*nx], Mr[0*nx*nx], 1, r_flag)
    fit[0] = 10000 * fit[0] / 2e+7

    schwefel_func(x, fit[1], nx, Os[1*nx], Mr[1*nx*nx], 1, r_flag)

    griewank_func(x, fit[2], nx, Os[2*nx], Mr[2*nx*nx], 1, r_flag)
    fit[2] = 1000 * fit[2] / 100

    rosenbrock_func(x, fit[3], nx, Os[3*nx], Mr[3*nx*nx], 1, r_flag)

    rastrigin_func(x, fit[4], nx, Os[4*nx], Mr[4*nx*nx], 1, r_flag)
    fit[4] = 10000 * fit[4] / 1e+3

    cf_cal(x, f, nx, Os, delta, bias, fit, cf_num)

def hf06(x, f, nx, Os, Mr, S, s_flag, r_flag):  # Hybrid Function 6
    cf_num = 5
    fit = [0.0] * cf_num
    G = [0] * cf_num
    G_nx = [0] * cf_num
    Gp = [0.3, 0.2, 0.2, 0.1, 0.2]

    tmp = 0
    for i in range(cf_num - 1):
        G_nx[i] = ceil(Gp[i] * nx)
        tmp += G_nx[i]
    G_nx[cf_num - 1] = nx - tmp

    for i in range(1, cf_num):
        G[i] = G[i - 1] + G_nx[i - 1]

    sr_func(x, z, nx, Os, Mr, 1.0, s_flag, r_flag)  # shift and rotate

    y = [0.0] * nx
    for i in range(nx):
        y[i] = z[S[i] - 1]

    katsuura_func(y[G[0]:G[0]+G_nx[0]], fit[0], G_nx[0], Os, Mr, 0, 0)
    happycat_func(y[G[1]:G[1]+G_nx[1]], fit[1], G_nx[1], Os, Mr, 0, 0)
    grie_rosen_func(y[G[2]:G[2]+G_nx[2]], fit[2], G_nx[2], Os, Mr, 0, 0)
    schwefel_func(y[G[3]:G[3]+G_nx[3]], fit[3], G_nx[3], Os, Mr, 0, 0)
    ackley_func(y[G[4]:G[4]+G_nx[4]], fit[4], G_nx[4], Os, Mr, 0, 0)

    f[0] = sum(fit)




def cf07(x, f, nx, Os, Mr, r_flag):  # Composition Function 4
    cf_num = 6
    fit = [0.0] * cf_num
    delta = [10, 20, 30, 40, 50, 60]
    bias = [0, 300, 500, 100, 400, 200]

    hgbat_func(x, fit[0], nx, Os[0*nx:0*nx+nx], Mr[0*nx*nx:0*nx*nx+nx*nx], 1, r_flag)
    fit[0] = 10000 * fit[0] / 1000

    rastrigin_func(x, fit[1], nx, Os[1*nx:1*nx+nx], Mr[1*nx*nx:1*nx*nx+nx*nx], 1, r_flag)
    fit[1] = 10000 * fit[1] / 1e+3

    schwefel_func(x, fit[2], nx, Os[2*nx:2*nx+nx], Mr[2*nx*nx:2*nx*nx+nx*nx], 1, r_flag)
    fit[2] = 10000 * fit[2] / 4e+3

    bent_cigar_func(x, fit[3], nx, Os[3*nx:3*nx+nx], Mr[3*nx*nx:3*nx*nx+nx*nx], 1, r_flag)
    fit[3] = 10000 * fit[3] / 1e+30

    ellips_func(x, fit[4], nx, Os[4*nx:4*nx+nx], Mr[4*nx*nx:4*nx*nx+nx*nx], 1, r_flag)
    fit[4] = 10000 * fit[4] / 1e+10

    escaffer6_func(x, fit[5], nx, Os[5*nx:5*nx+nx], Mr[5*nx*nx:5*nx*nx+nx*nx], 1, r_flag)
    fit[5] = 10000 * fit[5] / 2e+7

    cf_cal(x, f, nx, Os, delta, bias, fit, cf_num)


def shiftfunc(x, y, nx, Os):
    for i in range(nx):
        y[i] = x[i] - Os[i]

def rotatefunc(x, y, nx, Mr):
    for i in range(nx):
        y[i] = 0
        for j in range(nx):
            y[i] += x[j] * Mr[i*nx+j]

def sr_func(x, sr_x, nx, Os, Mr, sh_rate, s_flag, r_flag):  # shift and rotate
    if s_flag == 1:
        if r_flag == 1:
            shiftfunc(x, y, nx, Os)
            for i in range(nx):  # shrink to the original search range
                y[i] = y[i] * sh_rate
            rotatefunc(y, sr_x, nx, Mr)
        else:
            shiftfunc(x, sr_x, nx, Os)
            for i in range(nx):  # shrink to the original search range
                sr_x[i] = sr_x[i] * sh_rate
    else:
        if r_flag == 1:
            for i in range(nx):  # shrink to the original search range
                y[i] = x[i] * sh_rate
            rotatefunc(y, sr_x, nx, Mr)
        else:
            for i in range(nx):  # shrink to the original search range
                sr_x[i] = x[i] * sh_rate


def asyfunc(x, y, nx, beta):
    for i in range(nx):
        if x[i] > 0:
            y[i] = pow(x[i], 1.0 + beta * i / (nx - 1) * pow(x[i], 0.5))


def oszfunc(x, xosz, nx):
    for i in range(nx):
        if i == 0 or i == nx - 1:
            if x[i] != 0:
                xx = math.log(abs(x[i]))
            if x[i] > 0:
                c1 = 10
                c2 = 7.9
            else:
                c1 = 5.5
                c2 = 3.1
            if x[i] > 0:
                sx = 1
            elif x[i] == 0:
                sx = 0
            else:
                sx = -1
            xosz[i] = sx * math.exp(xx + 0.049 * (math.sin(c1 * xx) + math.sin(c2 * xx)))
        else:
            xosz[i] = x[i]


def cf_cal(x, f, nx, Os, delta, bias, fit, cf_num):
    w = [0.0] * cf_num
    w_max = 0
    w_sum = 0

    for i in range(cf_num):
        fit[i] += bias[i]
        w[i] = 0
        for j in range(nx):
            w[i] += pow(x[j] - Os[i*nx+j], 2.0)
        if w[i] != 0:
            w[i] = pow(1.0 / w[i], 0.5) * exp(-w[i] / 2.0 / nx / pow(delta[i], 2.0))
        else:
            w[i] = INF
        if w[i] > w_max:
            w_max = w[i]

    for i in range(cf_num):
        w_sum += w[i]

    if w_max == 0:
        for i in range(cf_num):
            w[i] = 1
        w_sum = cf_num

    f[0] = 0.0
    for i in range(cf_num):
        f[0] = f[0] + w[i] / w_sum * fit[i]

def cec22_test_func(x, f, nx, mx, func_num):
    cf_num = 12
    global ini_flag, n_flag, func_flag, M, OShift, y, z, x_bound, SS

    if func_num < 1 or func_num > 12:
        print(f"\nError: Test function {func_num} is not defined.\n")
        return

    if ini_flag == 1:
        if n_flag != nx or func_flag != func_num:
            ini_flag = 0

    if ini_flag == 0:
        M = None
        OShift = None
        y = None
        z = None
        x_bound = None
        y = np.zeros(nx)
        z = np.zeros(nx)
        x_bound = np.full(nx, 100.0)

        if nx not in [2, 10, 20]:
            print("\nError: Test functions are only defined for D=2,10,20.\n")
            return

        if nx == 2 and func_num in [6, 7, 8]:
            print("\nError:  NOT defined for D=2.\n")
            return

        # Load Matrix M
        file_name = f"input_data/M_{func_num}_D{nx}.txt"
        try:
            with open(file_name, "r") as fpt:
                if func_num < 9:
                    M = np.loadtxt(fpt)
                else:
                    M = np.loadtxt(fpt).reshape(cf_num, nx, nx)
        except FileNotFoundError:
            print(f"\n Error: Cannot open {file_name} for reading \n")
            return

        # Load shift_data
        file_name = f"input_data/shift_data_{func_num}.txt"
        try:
            with open(file_name, "r") as fpt:
                if func_num < 9:
                    OShift = np.loadtxt(fpt)
                else:
                    OShift = np.loadtxt(fpt).reshape(cf_num, nx)
        except FileNotFoundError:
            print(f"\n Error: Cannot open {file_name} for reading \n")
            return

        # Load Shuffle_data
        if func_num >= 6 and func_num <= 8:
            file_name = f"input_data/shuffle_data_{func_num}_D{nx}.txt"
            try:
                with open(file_name, "r") as fpt:
                    SS = np.loadtxt(fpt, dtype=int)
            except FileNotFoundError:
                print(f"\n Error: Cannot open {file_name} for reading \n")
                return

        n_flag = nx
        func_flag = func_num
        ini_flag = 1

    for i in range(mx):
        if func_num == 1:
            zakharov_func(x[i*nx:(i+1)*nx], f[i], nx, OShift, M, 1, 1)
            f[i] += 300.0
        elif func_num == 2:
            rosenbrock_func(x[i*nx:(i+1)*nx], f[i], nx, OShift, M, 1, 1)
            f[i] += 400.0
        # ... and so on for the other cases
        else:
            print("\nError: There are only 10 test functions in this test suite!\n")
            f[i] = 0.0