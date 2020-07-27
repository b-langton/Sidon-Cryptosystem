load('sidon_cryptosystem.sage')
##Function that returns the coeffiecient xth variable in the yth equation in the system
def getLinearizedCoeff(y, x, k, r, matrixList, c_matrix):
    if y % 100 == 0 and x == 0: 
        print(y)
    t1,a1,t2,a2,d,lx,ly,jx,jy = getIndices(y,x,k,r)

    if t1*r*k + a1 > t2*r*k + a2: 
      
        return 0
    c = c_matrix[a1][a2][d]
    coeff = c*(matrixList[t1][jy][ly]*matrixList[t2][jx][lx] - matrixList[t1][jy][lx]*matrixList[t2][jx][ly])

    coeff2 = 0
    if t1*r*k + a1 != t2*r*k + a2:
        
        temp = t1
        t1 = t2
        t2 = temp

        temp = a1
        a1 = a2
        a2 = temp
        
        c2 = c_matrix[a1][a2][d]
        coeff2 = c2*(matrixList[t1][jy][ly]*matrixList[t2][jx][lx] - matrixList[t1][jy][lx]*matrixList[t2][jx][ly])
           
    return coeff + coeff2

##Helper function that breaks y and x into the 9 indexes that specify an equation and variable in the system
def getIndices(y,x, k, r): 
    y = int(y)
    x = int(x)
    
    ##Get the variable, treating x like a base rk number
    rk  = r*k
    t1 = x//rk^3
    rem = x% rk^3 
    a1 = rem//rk^2
    rem = rem%rk^2
    t2 = rem//rk
    rem = rem%rk
    a2= rem 
    
    
    d = y//((k^2-k)*(k^2 - k)/4)
    rem = y%((k^2-k)*(k^2 - k)/4)
    l = rem//int(((k^2 - k)/2)) + 1
    rem = rem%((k^2 - k)/2)
    j = rem + 1 
    
    ##Get the equation. Indices must follow that lx > ly and jx > jy   
    lx  = 1
    jx = 1
    ly = 0
    jy = 0
    continuel = True
    continuej = True
    startl = k*(k-1)/2
    startj = k*(k-1)/2
    for i in range(k-1):
        if continuel == True:
            startl -= k-1 - i
            if l> startl:
                lx = k-1 - i 
                ly = l - startl - 1
                continuel = False
        if continuej == True:
            startj -= k-1 - i
            if j> startj:
                jx = k-1 - i 
                jy = j - startj - 1
                continuej = False
    return int(t1), int(a1), int(t2), int(a2), int(d), int(lx), int(ly), int(jx), int(jy)

def minorAttack(q,k):
    r = 2
    basefield = GF(q)
    rk = r*k
    y, F, F_r, d, c = ConstructSidon2k(q, k)
    matrixList, sidonbasis, mult_table, F_r_basis , origbasis = publicKey(y,q,F,F_r)
    basefield = GF(q)

    ##Construct a new random basis for F_r 
    F_r_basis_new = Matrix(basefield, rk, lambda i,j: basefield.random_element())
    while F_r_basis_new.is_invertible() == False: 
         F_r_basis_new = Matrix(basefield, rk, lambda i,j: basefield.random_element())
    cob_matrix = F_r_basis_new.inverse()
    basiselements = []
    ##represent the multiplication table of the new basis over that basis
    for i in range(rk): 
        basiselements += [convertFromLong(F_r_basis_new[i], F, F_r)]
    c_matrix = [[0 for i in range(rk)] for j in range(rk)]
    new_mult_table = [[i*j for i in basiselements] for j in basiselements]
    for i in range(rk): 
        for j in range(rk): 
            element = convertToLong(new_mult_table[i][j])
            c_matrix[i][j] = element*cob_matrix
    ##Construct the coefficient matrix for the system. Not efficient but works for small k 
    eqn_mat = Matrix(basefield, (k^3*(k-1)^2)/2, rk^4, lambda i,j: getLinearizedCoeff(i,j, k, r, matrixList, c_matrix))
    R = PolynomialRing(GF(q),'x', 36, order = "degrevlex")
    indeterminates = R.gens()
    R.inject_variables()
    pivot_cols = eqn_mat.transpose().pivot_rows()
    eqn_mat = eqn_mat.matrix_from_rows(eqn_mat.pivot_rows())
    indices = []
    for i in range(eqn_mat.ncols()):
        if eqn_mat.transpose()[i] == zero_vector(eqn_mat.nrows()):
            continue
        indices += [i]
    eqn_mat = eqn_mat.transpose().matrix_from_rows(indices).transpose()
    monomial_order = []
    for i in indices:
        t1,a1,t2,a2,d,lx,ly,jx,jy =  getIndices(0,i,k,r)
        monomial_order += [indeterminates[t1*r*k + a1]*indeterminates[t2*r*k + a2]]
    ideallist = list(eqn_mat*vector(monomial_order).column())
    ideallist = [R(i) for i in ideallist]
    I = ideal(ideallist)
    return I
            