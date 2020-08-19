"""Constructs a Sidon Space in Gq(rk, k) for r > 2 
   Paramters are q, k, and r defining F_q, F_q^k, and F_q^(rk) as used in the sidon space paper
   returns y, F_qn element defining the sidon space, F which is F_q^k, and F_n which is F_q^(rk)
"""
def ConstructSidon(q, k, r): 
    
    
    ##Define the field F_p^k and the polynomial ring over that field
    assert(r!=2)
    F.<z> = GF(q^k)
    R.<a> = F[]
    
    ##Construct another extension field of degree r over this ring using an irreducible polynomial of degree r over F_q^k
    irred_poly1 = R.irreducible_element(r)
    irred_poly2 = R.irreducible_element(r)
    F_r.<x> = F.extension(irred_poly1)[]
    F_n = F.extension(irred_poly1)
    ##find a root of this in F_q^rk
    roots = irred_poly2(x).roots()

    y = roots[0][0]
    
    ##returns a tuple containing y, the subfield F_q^k and the field F_q^rk
    ##The Sidon space is defined as u + u^p*y for u in F 
    return y, F, F_n
"""
Constructs a Sidon Space in Gq(2k, k), q cannot equal 2. 
Parameters q and k are the same as in ConstructSidon. 
Returns the same things as constructSidon, but additionally returns b and c, the linear and constant coefficients of the 
irreducible quadratic polynomial that y is a root of. These are necessary for factoring products of elements from the sidon space. 
"""
def ConstructSidon2k(q, k): 
    assert(q!=2)
    
    ##Define the field F_p^k and the polynomial ring
    F.<z> = GF(q^k) 
    R.<a> = F[] 
    
    p = F.characteristic()
    ##Find and irreducible polynomial of degree 2 with the constant term not in W_q-1 
    irred_poly = R.irreducible_element(2)
    iterations = 0
   
    while irred_poly.list()[0]^((q^k - 1)/(q-1)) == 1 and iterations < 100: 
        irred_poly = R.irreducible_element(2)
        iterations += 1
    
    assert(iterations != 100)
    
  
    ##Construct the second extension field as a polynomial ring over the first modded out by this irreducible element
    F_n  = F.extension(irred_poly)
    F_2k.<x> = F.extension(irred_poly)[]
    roots = irred_poly(x).roots() 
    
    ##Return a root of the polynomial as well as info about the two fields 
    y = roots[0][0]
    
    return y, F, F_n, irred_poly.list()[1], irred_poly.list()[0]
    
"""Factoring algorithm for Sidon spaces with r > 2 

    Parameters: 
            product: the product of two elements to be factored 
            y: The root of the irreducible polynomial that defines the sidon space 
            q: size of the base field
            F: F_q^k as described in the paper 
            F_r: F_q^n, the extension field the sidon space resides in

    Returns the two elements that were originally multiplied to get the product 
"""
def factor(product, y, q, F, F_r): 
    r = len(list(F_r.modulus())) - 1 
    p = F.characteristic()
    assert(r > 2)
    
    ##construct a basis from y 
    basis = [vector(y^n) for n in range(r)]
    mat = matrix(basis)
    
    ##calculate the change of basis matrix to go from the standard basis to our new basis, and compute 
    ##the representation of the product in that basis 
    cob_mat = mat.inverse() 
    product_representation = vector(product)*cob_mat
    
    ##Find the roots of the polynomial we derive from product_represntation, and then compute the original u and v
    assert(i == 0 for i in product_representation[3:])
    product_representation = product_representation[0:3]
    F_.<x> = F[]
    poly = F_(list(product_representation))
    roots = poly.roots()
    
    ##returns u and v up to multiplication from F_p, not the original things multiplied 
    if len(roots) == 1: 
        ans1 = (-1*1/roots[0][0]).nth_root(q - 1)
        ans2 = (-1*1/roots[0][0]).nth_root(q - 1)
        
    else: 
        ans1 = (-1*1/roots[0][0]).nth_root(F.characteristic() - 1)
        ans2 = (-1*1/roots[1][0]).nth_root(F.characteristic() - 1)
    
    if product/((ans1 + ans1^q*y)*(ans2 + ans2^q*y)) != 1: 
        ans1 = ans1*product/((ans1 + ans1^q*y)*(ans2 + ans2^q*y))
    return ans1, ans2

"""Factoring algorithm for when r = 2. 

    Parameters: 
        The first 5 are identical to the factor method 
        b: the coefficient of the linear term of the irrducible quadratic polynomial used to define the sidon space 
        c: the constant term of the same irreducible quadratic polynomial
    
    Returns the same thing is factor. 
"""
def factor2(product, y,q,F, F_r, b, c): 
    ##factoring algorithm for Sidon space with r = 2 
    p = F.characteristic()
    k = len(F.modulus().list()) - 1
    
    ##construct a basis for F_r that we will use to extract usefull info
    basis = [vector(y^n) for n in (0, 1)]
    mat = matrix(basis)
    cob_mat = mat.inverse() 
    product_representation = vector(product)*cob_mat
    
    
    ##Figure out the linear transformation x - cx^q 
    basis_transf = []
    identity = matrix.identity(k)
    
    for i in range(k): 
        transformed = vector(F(list(identity[i])) - c*F(list(identity[i]))^q)
        basis_transf += [transformed]
        
    transformation = matrix(basis_transf).inverse() 
    
    ##invert this transformation and use it to calculate uv     
    uv = vector(product_representation[0])*transformation
    uv = F(uv)
 
    ##Calculate the last two terms and create a quadratic polynomial
    F_.<s> = F[] 
    second_term = product_representation[1] + b*(uv)^q
    third_term = uv^q
    poly = F_([vector(uv), vector(second_term), vector(third_term)])
    
    ##get the roots of the polynomial and extract u and v
    roots = poly.roots()
    
    ##returns u and v up to multiplication from F_q, not the original things multiplied 
    if len(roots) == 1: 
        ans1 = (-1*1/roots[0][0]).nth_root(q - 1)
        ans2 = (-1*1/roots[0][0]).nth_root(q - 1)
        
    else: 
        ans1 = (-1*1/roots[0][0]).nth_root(q - 1)
        ans2 = (-1*1/roots[1][0]).nth_root(q - 1)
    
    if product/((ans1 + ans1^q*y)*(ans2 + ans2^q*y)) != 1: 
        ans1 = ans1*product/((ans1 + ans1^q*y)*(ans2 + ans2^q*y))
    return ans1, ans2
    
"""Generates the public key (and associated private key) using info from the given sidon space 

    Parameters: 
            y: The F_q^n element that defines the sidon space 
            q: Size of the basefield 
            F: F_q^k as used in paper 
            F_r: F_q^n, the extension field the sidon space lives in. 
    
    Returns: 
            matrixList: The public key, coefficient matrices of the multiplication table of the sidon space 
            sidonbasis: The basis of the sidon space used, as a list of F_q^n elements 
            mult_table: The multiplication table of sidonbasis
            F_r_basis: The basis of F_q^n used to get the coefficient matrices 
            origbasis: Basis of F_q^k used to construct sidonbasis 
"""
def publicKey(y,q, F, F_r): 
    p = F.characteristic()
    k = len(F.modulus().list()) - 1
    rk = (len(F_r.modulus().list()) - 1)*k
    basefield = GF(q)
    iterations = 0
    iterations2 = 0
    
    ##Construct bases for the sidon space as well as F_r
    sidonbasis = Matrix(basefield, k, lambda i,j: basefield.random_element())
    while sidonbasis.is_invertible() == False and iterations < 100:
        sidonbasis = Matrix(basefield, k, lambda i,j: basefield.random_element())
        iterations += 1
    F_r_basis = Matrix(basefield, rk, lambda i,j: basefield.random_element())
    while F_r_basis.is_invertible() == False and iterations2 < 100: 
         F_r_basis = Matrix(basefield, rk, lambda i,j: basefield.random_element())
         iterations2 += 1
    
    assert(iterations<100 and iterations2<100)
    
    ##Get the multiplication table from the basis of the sidon space
    v = [F(list(sidonbasis[i])) for i in range(k)]
    origbasis = sidonbasis
    sidonbasis = [j + j^q*y for j in v]
    mult_table = vector(sidonbasis).column()*vector(sidonbasis).row()
    cob_matrix = F_r_basis.inverse()
    print(cob_matrix)
    vec_list = [[0 for i in range(k)] for j in range(k)]
    
    ##Generate the public key M(V,B)
    for i in range(k): 
        for j in range(k):
            element = mult_table[i][j]
            long_representation = convertToLong(element)
            vec_list[i][j] = list(long_representation*cob_matrix) 
    matrixlist = [0 for i in range(rk)]
    
    ##Return the public key 
    for i in range(rk): 
        matrixlist[i] = Matrix(basefield, k, lambda l,j: vec_list[l][j][i])
    return matrixlist, sidonbasis, mult_table, F_r_basis, origbasis

"""Converts an element in Fq_n to a vector over F_q
    Parameter: An element of F_q^n 
    Returns: A vector over F_q of length n representing that element 
"""
def convertToLong(element):
    long_representation = []
   
    for l in range(len(list(element))): 
        long_representation += list(vector(list(element)[l]))
    long_representation = vector(long_representation)
    long_representation
    return long_representation

"""Converts an element from F_r represented as a vector over F_q to an element of F_r

    Parameters: 
            element: a vector of length n over F_q 
            F: F_q^k
            F_r: The extension field of F_q^k, F_q^n 
    Returns: The corresponding element of F_r 
"""
def convertFromLong(element, F, F_r):

    
    k = len(F.modulus().list()) - 1
    r = len(F_r.modulus().list()) - 1
    final_list = []
    
    for i in range(r): 
        final_list += [F(list(element)[i*k:(i+1)*k])]
    return(F_r(final_list))
    
def getIndices(y,x, k, r): 
    y = int(y)
    x = int(x)
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
            
            
        
def getLinearizedCoeff(y, x, k, r, matrixList, c_matrix):
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
    
    