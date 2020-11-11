load('sidon_cryptosystem.sage')
'''
method that produces system of homogeneous equations for major attack. 
returns a 3-d array, every entry of which is equal to 0.
'''
def minrankmajor(q,k):
    basefield = GF(q)
    r = 2
    rk = r*k
    y, F, F_r, d, c = ConstructSidon2k(q, k)
    ##Construct the attacker's F_r
    y2, F2, F_r2, d2, c2 = ConstructSidon2k(q, k)
    matrixList, sidonbasis, mult_table, F_r_basis , origbasis = publicKey(y,q,F,F_r)

    ##Construct the attacker's basis of F and F_r
    iterations = 0
    sidonbasis = Matrix(basefield, k, lambda i,j: basefield.random_element()*0)

    ##basis of F
    while sidonbasis.is_invertible() == False and iterations < 500:
        element = F.random_element()
        sidonbasis = Matrix([list(vector(element^(q^i))) for i in range(k)])
        iterations += 1
    v = [F2(list(sidonbasis[i])) for i in range(k)]
    sidonbasis = [j*y2*y2^-1 for  j in v]
    sidonbasismat = Matrix([convertToLong(i) for i in sidonbasis])

    ##Extend F to F_r
    F_r2basis = Matrix(basefield, k, lambda i,j: basefield.random_element()*0)
    while F_r2basis.is_invertible() == False and iterations < 500:
        F_r2basis = sidonbasismat
        stackmat = Matrix(basefield, rk - k, rk, lambda i, j: basefield.random_element())
        F_r2basis = F_r2basis.stack(stackmat)
        iterations += 1
        print(F_r2basis.rank(), F_r2basis.nrows(), F_r2basis.ncols())

    F_r2basisels =  [convertFromLong(F_r2basis[i], F2, F_r2) for i in range(rk)]
    g = ['g' + str(i) for i in range(rk)]
    ##create a polynomial ring over the variables in order to do symbolic calculations 
    b = [["b" + str(j) + str(i) for i in range(rk)] for j in range(rk)]
    u = [["u" + str(j) + str(i) for i in range(k) for j in range(k)]]
    bvec = []
    uvec = []
    for i in b: 
        bvec += i
    for i in u: 
        uvec += i
    R = PolynomialRing(F_r2, g + bvec + uvec)
    indeterminates = R.gens()
    R.inject_variables()

    for i in range(rk): 
        g[i] = R(g[i])
    for i in range(k*k): 
        uvec[i] = R(uvec[i])
    for i in range(rk*rk):
        bvec[i] = R(bvec[i])
    y2 = vector(g).row()*vector(F_r2basisels).column()
    F_r2_q = vector([R(i^q) for i in F_r2basisels])
    v = [] 
    for i in range(k):
        v += vector(uvec[i*k:(i+1)*k]).row()*vector(F_r2basisels[0:k]).column() + y2*vector(uvec[i*k:(i+1)*k]).row()*vector(F_r2_q[0:k]).column() 
    B = [0 for i in range(rk)]
    for i in range(rk): 

        for j in range(rk):
            B[i] += R(bvec[i*rk + j])*R(F_r2_q[j])
    v_mult_table =[[v[i]*v[j] for i in range(len(v))] for j in range(len(v))]
    vectorized_eq_table = [[0 for i in range(k)] for j in range(k)]
    for i in range(k): 
        for j in range(k): 
            if i > j: 
                continue
            else: 
                eqvec = list(eq_table[i][j])
                sumvec = vector([0 for l in range(rk)])
                for term in eqvec: 
                    coeff = term[0]
                    monomial = term[1]
                    coeff = convertToLong(F_r2(coeff))
                    rep = coeff*F_r2basis.inverse()
                    sumvec += rep*monomial
            vectorized_eq_table[i][j] = sumvec
    return vectorized_eq_table