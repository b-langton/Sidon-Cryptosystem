# Sidon-Cryptosystem
Codebase for: 
N. Raviv, B. Langton, and I. Tamo, "Multivariate Public Key Cryptosystem from Sidon Spaces," Public Key Cryptography (PKC), 2021.
## Usage, Implementation
All of the base code is written in Sage 9.0 . You must have this installed to be able to use the code. 
The useful, documented code is all located in the "Sage Code" Folder. This includes the implementation of the Sidon cryptosystem, as well as the implementations of the 3 major attacks included in the paper. 
With the exception of the algebraic attack, these implementations do not solve the MP system generated, they only return the associated ideal. They can be solved by running your Groebner basis 
algorithm of choice on the ideal. 

The "Notebooks" folder contains the experiments and testing. Most of this code is not commented or documented. The important experiments and code will be put into 
a separate folder. 

### Code Examples-- Constructing and using the sidon cryptosystem 
Constructing the cryptosystem is straightforward. We make the sidon space and then pass the necessary info into the constructor for the cryptosystem: 
```
load('sidon_cryptosystem.sage')
q = 7
basefield = GF(q)
k = 3
y, F, F_r, d, c = ConstructSidon2k(q, k)
key, sidonbasis, mult_table, F_r_basis , origbasis = publicKey(y,q,F,F_r)
```

All the other outputs are helpful in decrypting an encrypted message,
but the key (public key) is the only thing you need to send messages: 

```
a = vector(GF(q), [2,2,0])
b = vector(GF(q), [1,1,3])
message = encrypt(a,b,key)
```
a and b have now been encrypted into message by the public key. 
Technically the message being sent is the equivalence class containing a^t * b and b^t * a. 
To decrypt the message, we simply call: 

```
mat = decrypt(message, y,q, F, F_r, F_r_basis, sidonbasis, origbasis, d,c)
```

### Attacks 
While actually solving the systems of equations that appear in these attacks is difficult, framework is provided to obtain the multivariate polynomial systems. 
There are 3 main attacks implemented here: Direct (Section 5.2), Minrank Minor (Section 4.2), and the "Major" attack (Section 5.1). Getting the ideals associated to the systems is simple, just call the method with 
the appropriate parameters: 
```
##Ideal for algebraic attack
I = algebraicAttack(q,k,a,b,matrixList)

##Ideal for minrank minor attack. Uses a random sidon cryptosystem. 
I2 = minrankMinor(q,k)
```



