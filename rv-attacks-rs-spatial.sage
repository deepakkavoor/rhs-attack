import random
# random.seed(0)
import time
from quadtree import get_quadtree, xmin, xmax, ymin, ymax

# Overview:
# RS knows EN_ij values and p; RV knows p, alpha
# RS sends EN to RV
# Goal: RV can eliminate s, alpha, a_jh from N1, N2 to find quadrant vertices
# Refer to paper for exact equations 


# Default parameters from TRACE paper:
# number of quad tree nodes m: 50 or 100
# k1 = 512, k2 = 160, k3 = 75, k4 = 75 or
# k1=2048, k2 = 1000, k3 = 400, k4 = 400
# prime |p| = k1 = 512
# prime |alpha| = k2 = 160
# s is in Zp*
# |a_jh| = k3 = 75 (j=1..4, h=1..6)
# |r_ij| = k4 = 75
(k1,k2,k3,k4) = (512,160,75,75)
#(k1,k2,k3,k4) = (2048,1000,400,400)
m_max = 53
#m_max = 101


def modinv(x,p):
    # when p is prime. return inverse of x modulo p
    return pow(x, p-2, p)


# to generate large primes, use this website: https://asecuritysite.com/encryption/random3
# example (for k1=512, k2=160):
# p = 2610012493187751553850229270260260188412785728466825936469473849183260172904846529218774317876204989248256766513760785264139183325575284012588127483034783
# alpha = 805054899696897334540366748238814746338336769857

def rv_attack_rs_helper(EN1, EN2, p, alpha):

    # RV is given EN1, EN2 and tries to find N1, N2
    # create matrices A,b; goal is to later combine and solve AX=b 
    # as of now A:28x16, X:16x1, b:28x1, but partition A,b so that they correspond to their Ni variables
    # here X will have the 16 unknown quadrant vertices (x_Nij,y_Nij) (i=1,2;j=1...4) of N1,N2
    # first 8 elem of X (i.e. X[0:8]) will be from N1; X[8:16] will be from N2
    # X[0:8] will be of the format x_N11,y_N11,x_N12,y_N12,...,x_N14,y_N14
    # (assume 0 indexing for j. i=1,2;j=0..3). this means
    # x_Nij will correspond to X[8*(i-1)+2*j]
    # y_Nij will correspond to X[8*(i-1)+2*j+1]

    A = [[0 for _ in range(16)] for _ in range(28)]
    b = [0 for _ in range(28)]
    P = GF(p)
    i=0

    for j in range(4):
        # choose which row of i to fill
        # i should be a counter from 0 to 27
        jnxt = (j+1)%4

        # eqn (8) from paper
        A[i][8*(1-1)+2*jnxt+1] =  P(EN1[j][0]-EN2[j][0])
        A[i][8*(2-1)+2*j] = P(EN1[j][3]-EN2[j][3])
        b[i] = P(EN1[j][4]-EN2[j][4])
        i+=1

        # eqn (9) from paper
        A[i][8*(2-1)+2*jnxt+1] =  P(EN1[j][0]-EN2[j][0])
        A[i][8*(1-1)+2*j] = P(EN1[j][3]-EN2[j][3])
        b[i] = P(EN1[j][4]-EN2[j][4])
        i+=1

        # eqn (10) from paper
        A[i][8*(1-1)+2*jnxt] =  P(EN1[j][1]-EN2[j][1])
        A[i][8*(2-1)+2*j+1] = P(EN1[j][2]-EN2[j][2])
        b[i] = P(EN1[j][5]-EN2[j][5])
        i+=1

        # eqn (11) from paper
        A[i][8*(2-1)+2*jnxt] =  P(EN1[j][1]-EN2[j][1])
        A[i][8*(1-1)+2*j+1] = P(EN1[j][2]-EN2[j][2])
        b[i] = P(EN1[j][5]-EN2[j][5])
        i+=1

        # eqn (12) from paper
        A[i][8*(1-1)+2*j] =  P(EN1[j][1]-EN2[j][1])
        A[i][8*(2-1)+2*j] =  P(-EN1[j][1]+EN2[j][1])
        A[i][8*(1-1)+2*j+1] = P(-EN1[j][0]+EN2[j][0])
        A[i][8*(2-1)+2*j+1] = P(EN1[j][0]-EN2[j][0])
        b[i] = P(0)
        i+=1

        # eqn (13) from paper
        A[i][8*(1-1)+2*j+1] = P(EN1[j][2]-EN2[j][2])
        A[i][8*(2-1)+2*j+1] = P(-EN1[j][2]+EN2[j][2])
        A[i][8*(1-1)+2*jnxt] =  P(-EN1[j][1]+EN2[j][1])
        A[i][8*(2-1)+2*jnxt] =  P(EN1[j][1]-EN2[j][1])
        b[i] = P(0)
        i+=1

        # eqn (14) from paper
        A[i][8*(1-1)+2*jnxt] =  P(EN1[j][3]-EN2[j][3])
        A[i][8*(2-1)+2*jnxt] =  P(-EN1[j][3]+EN2[j][3])
        A[i][8*(1-1)+2*jnxt+1] = P(-EN1[j][2]+EN2[j][2])
        A[i][8*(2-1)+2*jnxt+1] = P(EN1[j][2]-EN2[j][2])
        b[i] = P(0)
        i+=1

    # first 8 cols of A correspond to N1
    A_N1 = [Ai[0:8] for Ai in A]
    A_N2 = [Ai[8:16] for Ai in A]

    # size of A_Ni is 28x8, size of b is 28x1

    return A_N1,A_N2,b


def rv_attack_rs_three_quad(EN1,EN2,EN3,p,alpha):
    # (just for experimental stats; use the function taking four ENi's)
    # try to find quadrant vertices of three Ni's 
    # call rv_attack_rs_helper for every two pair of ENi's
    # later create a (28+28+28)x(24) matrix A, (28+28+28)x1 matrix b
    # and solve AX=b to find 24x1 matrix X containing vertices of N1,N2,N3

    dummy_Ai = [0 for _ in range(8)]
    A = []

    A_N1,A_N2,b_N1N2 = rv_attack_rs_helper(EN1,EN2,p,alpha)
    for i in range(28):
        A.append(A_N1[i]+A_N2[i]+dummy_Ai)

    A_N2,A_N3,b_N2N3 = rv_attack_rs_helper(EN2,EN3,p,alpha)
    for i in range(28):
        A.append(dummy_Ai+A_N2[i]+A_N3[i])

    A_N1,A_N3,b_N1N3 = rv_attack_rs_helper(EN1,EN3,p,alpha)
    for i in range(28):
        A.append(A_N1[i]+dummy_Ai+A_N3[i])


    b = b_N1N2 + b_N2N3 + b_N1N3

    # print('A',A)
    # print(len(A),len(A[0]))
    # print('b',b)
    # print('Solving')

    A = matrix(GF(p), A)
    b = matrix(GF(p), b).transpose()
    X = A.solve_right(b)

    # starting from topleft clockwise
    X = [Xi[0] for Xi in X]
    # print('X', X)
    print('rank(A):',rank(A))
    return A,X, b


def rv_attack_rs_four_quad(EN1,EN2,EN3,EN4,p,alpha):
    # try to find quadrant vertices of three Ni's
    # call rv_attack_rs_helper for every two pair of ENi's
    # later create a (28+28+28)x(24) matrix A, (28+28+28)x1 matrix b
    # and solve AX=b to find 24x1 matrix X containing vertices of N1,N2,N3

    dummy_Ai = [0 for _ in range(8)]
    A = []

    A_N1,A_N2,b_N1N2 = rv_attack_rs_helper(EN1,EN2,p,alpha)
    for i in range(28):
        A.append(A_N1[i]+A_N2[i]+dummy_Ai+dummy_Ai)

    A_N2,A_N3,b_N2N3 = rv_attack_rs_helper(EN2,EN3,p,alpha)
    for i in range(28):
        A.append(dummy_Ai+A_N2[i]+A_N3[i]+dummy_Ai)

    A_N3,A_N4,b_N3N4 = rv_attack_rs_helper(EN3,EN4,p,alpha)
    for i in range(28):
        A.append(dummy_Ai+dummy_Ai+A_N3[i]+A_N4[i])

    A_N1,A_N3,b_N1N3 = rv_attack_rs_helper(EN1,EN3,p,alpha)
    for i in range(28):
        A.append(A_N1[i]+dummy_Ai+A_N3[i]+dummy_Ai)

    A_N2,A_N4,b_N2N4 = rv_attack_rs_helper(EN2,EN4,p,alpha)
    for i in range(28):
        A.append(dummy_Ai+A_N2[i]+dummy_Ai+A_N4[i])

    A_N1,A_N4,b_N1N4 = rv_attack_rs_helper(EN1,EN4,p,alpha)
    for i in range(28):
        A.append(A_N1[i]+dummy_Ai+dummy_Ai+A_N4[i])


    b = b_N1N2 + b_N2N3 + b_N3N4 + b_N1N3 + b_N2N4 + b_N1N4

    # print('A',A)
    # print(len(A),len(A[0]))
    # print('b',b)
    # print('Solving')

    A = matrix(GF(p), A)
    b = matrix(GF(p), b).transpose()
    X = A.solve_right(b)

    # starting from topleft clockwise
    X = [Xi[0] for Xi in X]
    # print('X', X)
    print('rank(A):',rank(A))
    return A,X, b
    


def rv_attack_rs_spatial():

    # Generated by RS: p, alpha, m, s, a_ji, N, EN
    # Fix some parameters
    p = random_prime(2**k1-1, False, 2**(k1-1))
    alpha = random_prime(2**k2-1, False, 2**(k2-1))

    s = random.randint(0,p-1)
    sinv = modinv(s,p)

    print('RS generating a_jh')
    a = []
    for j in range(4):
        a_j = []
        for h in range(6):
            a_j.append(random.randint(1<<(k3-1),1<<k3-1))
        a.append(a_j)

    print('RS generating N_i')
    N = get_quadtree(m_max)
    print('len of N', len(N))
    m = len(N)

    print('RS generating EN_i')
    EN = []
    for i in range(m):
        EN_i = []
        for j in range(4):
            j_nxt = (j + 1 ) % 4
            [x_Nij, y_Nij] = N[i][j]
            [x_Nijnxt, y_Nijnxt] = N[i][j_nxt]
            EN_ij1 = s*(x_Nij * alpha + a[j][0]) % p
            EN_ij2 = s*(y_Nij * alpha + a[j][1]) % p
            EN_ij3 = s*(x_Nijnxt * alpha + a[j][2]) % p
            EN_ij4 = s*(y_Nijnxt * alpha + a[j][3]) % p
            EN_ij5 = s*(x_Nij * y_Nijnxt * alpha + a[j][4]) % p
            EN_ij6 = s*(x_Nijnxt * y_Nij * alpha + a[j][5]) % p
            EN_ij = [EN_ij1, EN_ij2, EN_ij3, EN_ij4, EN_ij5, EN_ij6]
            EN_i.append(EN_ij)
        EN.append(EN_i)

    left = [i for i in range(0,m)]  # indices remaining Nis that are to be recovered
    done = []  # Nis that are already recovered
    success = True
    j=1
    while len(left)>0:

        if len(left) < 4:
            while len(left) < 4:
                left.append(random.choice(done))

        print('Iteration: j=',j,'no. left:',len(left),'out of',m)
        start_j = time.time()

        # RV gets four random ENi and tries to find corresponding four Ni
        j1 = random.choice(left)
        left.remove(j1)
        j2 = random.choice(left)
        left.remove(j2)
        j3 = random.choice(left)
        left.remove(j3)
        j4 = random.choice(left)
        left.remove(j4)
        EN1 = EN[j1]
        EN2 = EN[j2]
        EN3 = EN[j3]
        EN4 = EN[j4]
        N1=N[j1]
        N2=N[j2]
        N3 = N[j3]
        N4 = N[j4]
        #A,X,b = rv_attack_rs_three_quad(EN1,EN2,EN3,p,alpha)
        A,X,b = rv_attack_rs_four_quad(EN1,EN2,EN3,EN4,p,alpha)

        # for verifying
        X_orig = []
        #N_curr = N1+N2+N3
        N_curr = N1+N2+N3+N4
        for Ni in N_curr:
            X_orig.append(Ni[0])
            X_orig.append(Ni[1])
        print('X_orig==X',X_orig==X)
        success = success & (X_orig == X)
        # print('X',X)
        # print('X_orig',X_orig)
        X_orig = matrix(GF(p), X_orig).transpose()
        b_orig = A*X_orig
        # print('b_orig:',b_orig)
        # print('b',b) 
        print('b_orig==b',b_orig==b)
        print('Time for curr (sec):',time.time()-start_j)

        done.append(j1)
        done.append(j2)
        done.append(j3)
        done.append(j4)

        j+=1

    return success


if __name__ == '__main__':

    num_trials = 10
    num_success = 0
    max_time = 0
    avg_time = 0
    print('Params m:',m_max)
    print('Params (k1,k2,k3,k4):',(k1,k2,k3,k4))
    for i in range(1, num_trials+1):
        start = time.time()
        success = rv_attack_rs_spatial()
        if success:
            print('====================Trial', i, 'successful====================')
            num_success += 1
            max_time = max(max_time, time.time()-start)
            avg_time += (time.time()-start)
        else:
            print('====================Trial', i, 'unsuccessful====================')
        print('Max time for successful attack (sec):', max_time)
    if num_success > 0:
        avg_time /= num_success
    print('Final success rate:', num_success, '/', num_trials, '(', num_success/num_trials*100, '%)')
    print('Avg time for successful attack (sec):', avg_time)


