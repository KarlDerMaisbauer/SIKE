# -*- coding: utf-8 -*-
"""
Created on Wed Aug 17 14:09:53 2022

@author: David1
"""
import sympy
import copy
#Parameters

#434
#e1=0xD8
#e2=0x89

#503
#e1 = 0xFA
#e2 = 0x9F

#610
e1 = 0x131
e2 = 0xC0

p = (2**e1)*(3**e2)-1

#p=(2**63)*(3**41)-1
from collections import deque
class fp2:
    """ Element in the arithmetic field FP2(X^2+1).
        
        Parameters
        ----------
        arg1 : integer
        First element from Fp
        
        arg2 : integer
        Second element from Fp
        """
    def __init__(self, a0=0, a1=0):
        self.a0 = a0
        self.a1 = a1
    def add(self, P, Q):
        self.a0 = (P.a0 + Q.a0)
        self.a1 = (P.a1 + Q.a1)
    def sous(self, P, Q):
        self.a0 = (P.a0 - Q.a0)
        self.a1 = (P.a1 - Q.a1)
    def mul(self, P, Q):
        self.a0 = (P.a0 * Q.a0 - P.a1 * Q.a1)%p # 
        self.a1 = (P.a0 * Q.a1 + P.a1 * Q.a0)%p
    def neg(self, P):
        self.a0 = -P.a0
        self.a1 = -P.a1
    def inv(self, P):
        part = (P.a0 * P.a0 + P.a1 * P.a1)%p
        self.a0 = (P.a0 * sympy.mod_inverse(part, p))%p
        self.a1 = (-P.a1 * sympy.mod_inverse(part, p))%p
    def pow2(self, P):
        self.a0 = (P.a0 * P.a0 - P.a1 * P.a1 + 4*p**2)%p
        self.a1 = (P.a0 * P.a1 + P.a1 * P.a0)%p
    def cpy(self, P):
        self.a0 = copy.deepcopy(P.a0)
        self.a1 = copy.deepcopy(P.a1)
    def __str__(self):
        return "%s + %si" %(hex(self.a0),hex(self.a1))


    # added functions for testing
    def __eq__(self, eq):
        if isinstance(eq, self.__class__):
            return self.a0 == eq.a0 and self.a1 == eq.a1
        else:
            return False
    def mod(self):
            self.a0 = self.a0 % p
            self.a1 = self.a1 % p
    def mul_test(self, P, Q):
        self.a0 = (P.a0 * Q.a0 - P.a1 * Q.a1) 
        self.a1 = (P.a0 * Q.a1 + P.a1 * Q.a0)

class point:
    """ Definition and operation of a point in a Montgomery curves.
        
        Parameters
        ----------
        arg1 : fp2
        X coordinate in projective

        arg2 : fp2
        Z coordinate in projective
        
        x = X / Z
        """
    #functions for testing
    def mod(self):
        self.X.mod()
        self.Z.mod()
    
    def __init__(self, x = fp2(), z = fp2()):
        self.X = copy.deepcopy(x)
        self.Z = copy.deepcopy(z)

    def double(self,P,para):
        """
        Doubling of a Montgomery point in projective coordinates (X:Z).

        Parameters
        ----------
        P : Point
        
        para : Point
        Parameters of the curve : (A+24 : C24) = (A + 2C : 4C)

        Returns
        -------
        Self.

        """
        #initiate
        t0 = fp2()
        t1 = fp2()
        tb = fp2()
        xf = fp2()
        za = fp2()
        zf = fp2()
        #Doubling
        t0.sous(P.X,P.Z)#t0 = X - Z
        #print(t0)
        t1.add(P.X,P.Z)#t1 = X + Z
        tb.pow2(t0)#t0 = t0² = (X - Z)²
        t0.cpy(tb)
        tb.pow2(t1)#t1 = t1² = (X + Z)²
        t1.cpy(tb)
        za.mul(para.Z,t0)#za = C24*t0 = C24 * (X - Z)²
        xf.mul(za,t1)#xf = za*t1 = C24 * (X - Z)² * (X + Z)²
        #print(t1,t0)
        tb.sous(t1,t0)#t1 = t0 - t1 = (X + Z)² - (X - Z)²
        t1.cpy(tb)
        #print(t1)
        t0.mul(para.X,t1)#t0 = A24*t1 = A24 * (X + Z)² * (X - Z)²
        tb.add(za,t0)#za = za + t0
        za.cpy(tb)
        zf.mul(za,t1)#zf = za * t1
        #tb.pow2(P.Z)
        #t0.mul(xf,tb)
        #print("double :",conv(point(t0,zf)),t0," ori :",conv(point(xf,zf)),xf,zf,tb)
        #xf.cpy(t0)
        self.Z.cpy(zf)
        self.X.cpy(xf)
        del t0,t1,tb,xf,za,zf
        
    def puissance2_k(self,k,P,para):
        """
        k times doubling of a Montgomery point in projective coordinates (X:Z):
        <- [2^k] * P
        
        Parameters
        ----------
        arg1 : self
        Self.
        
        arg2 : integer
        Power of 2 to multiply.
        
        arg3 : Point
        The point to double k times.
        
        arg4 : Point
        Parameters of the curve : (A+24 : C24) = (A + 2C : 4C)
        """
        #INIT
        R = point()
        R.cpy(P)
        #case k = 0 : id function 
        if (k == 0):
            self.Z.cpy(R.Z)
            self.X.cpy(R.X)
        #case k = 1 : simple point doubling
        elif(k == 1):
            R.double(P,para)
            self.Z.cpy(R.Z)
            self.X.cpy(R.X)
        else:
            Rc = point()
            for i in range(k):
                Rc.double(R,para)
                R.cpy(Rc)
                #print(" 2^",i+1,"S =",R)
            self.Z.cpy(R.Z)
            self.X.cpy(R.X)
        del R
    
    def TPL(self,P,para):
        """
        Tripling of a Montgomery point in projective coordinates (X:Z).

        Parameters
        ----------
        P : Point
        projective Montgomery x-coordinates P = (X:Z), where x=X/Z
        
        para : Point
        Montgomery curve constants (A24+ : A24-) = (A+2C : A-2C).

        """
        #INIT
        t0 = fp2()
        tb = fp2()
        t1 = fp2()
        t2 = fp2()
        t3 = fp2()
        t4 = fp2()
        t5 = fp2()
        t6 = fp2()
        xf = fp2()
        zf = fp2()
        #computation
        t0.sous(P.X,P.Z)#t0 = Xp - Zp
        t2.pow2(t0)#t2 = t0² = (Xp - Zp)²
        t1.add(P.X,P.Z)#t1 = Xp + Zp
        t3.pow2(t1)#t3 = t1² = Xp + Zp
        t4.add(t0,t1)#t4 = t1 + t0 = Xp - Zp + Xp + Zp = 2 * Xp
        tb.sous(t1,t0)#t0 = t1 - t0
        t0.cpy(tb)
        t1.pow2(t4)#t1 = t4² = (2 * Xp)²
        tb.sous(t1,t3)#t1 = t1 - t3 = (Xp + Zp)(1 - (Xp + Zp))
        t1.cpy(tb)
        tb.sous(t1,t2)#t1 = t1 - t2 = (Xp + Zp)(1 - (Xp + Zp)) - (Xp - Zp)²
        t1.cpy(tb)
        t5.mul(t3,para.X)#t5 = t3 * A24+
        tb.mul(t5,t3)#t3 = t3 * t5
        t3.cpy(tb)
        t6.mul(t2,para.Z)#t6 = t2 * A24-
        tb.mul(t2,t6)#t2 = t2 * t6
        t2.cpy(tb)
        tb.sous(t2,t3)#t3 = t2 - t3
        t3.cpy(tb)
        t2.sous(t5,t6)#t2 = t5 - t6
        tb.mul(t1,t2)#t1 = t1 * t2
        t1.cpy(tb)
        t2.add(t3,t1)#t2 = t3 + t1
        tb.pow2(t2)# t2 = t2²
        t2.cpy(tb)
        xf.mul(t2,t4)#Xf = t2 * t4
        tb.sous(t3,t1)#t1 = t3 - t1
        t1.cpy(tb)
        tb.pow2(t1)#t1 = t1²
        t1.cpy(tb)
        zf.mul(t1,t0)#Zf = t1 * t0
        self.X.cpy(xf)
        self.Z.cpy(zf)
        
    def TPLe(self,k,P,para):
        """
        Computes [3^e](X:Z) on Montgomery curve with projective.


        Parameters
        ----------
        k : integer
        Number of point tripling to do.
        P : Point
        projective Montgomery x-coordinates P = (X:Z), where x=X/Z
        para : Point
        Montgomery curve constants (A24+ : A24-) = (A+2C : A-2C).

        Returns
        -------
        None.

        """
        #INIT
        R = point()
        R.cpy(P)
        #case k = 0 : id function 
        if (k == 0):
            self.Z.cpy(R.Z)
            self.X.cpy(R.X)
        #case k = 1 : simple point doubling
        elif(k == 1):
            R.TPL(P,para)
            self.Z.cpy(R.Z)
            self.X.cpy(R.X)
        else:
            Rc = point()
            for i in range(k):
                Rc.TPL(R,para)
                R.cpy(Rc)
                #print(" 2^",i+1,"S =",R)
            self.Z.cpy(R.Z)
            self.X.cpy(R.X)
        del R
    
    def ladder3pt(self,k,P,Q,QP,para):
        """
        Computes P + [k] * Q on Montgomery curve with projective.

        Parameters
        ----------
        k : integer
            DESCRIPTION.
        P : Point
        projective Montgomery x0-coordinates P = (X0:Z0), where x0=X0/Z0
        Q : Point
        projective Montgomery x1-coordinates Q = (X1:Z1), where x1=X1/Z1
        QP : Point
        projective Montgomery coordinates of the point Q - P
        para : Point
        Montgomery curve constants (A+2C: 4C).

        Returns
        -------
        None.

        """
        #INIT
        R0 = point()
        R1 = point()
        R2 = point()
        #S = copy of R
        S0 = point()
        S1 = point()
        S2 = point()
        A24 = point()
        a0 = fp2()
        #value assignement
        R0.cpy(Q)
        R1.cpy(P)
        R2.cpy(QP)
        #inverse of 4 in Fp
        inv4 = sympy.mod_inverse(4,p)
        #A24
        a0.add(para.X,fp2(2,0))
        A24.X.mul(a0,fp2(inv4,0))
        A24.Z = fp2(1,0)
        Tab = []
        r = k
        while (r!= 1):
            Tab.append(r%2)
            r=r//2
        Tab.append(1)
        n = len(Tab)
        #print(f'k is {hex(k)}')
        for i in range(n):
            condition = k>>i
            if condition&1 == 1:
                S0,S1 = DBLADD(R0,R1,R2,A24)
                R0.cpy(S0)
                R1.cpy(S1)
            else:
                S0,S2 = DBLADD(R0,R2,R1,A24)
                R0.cpy(S0)
                R2.cpy(S2)
        self.X.cpy(R1.X)
        self.Z.cpy(R1.Z)
        
    def isogeny_2_point(self,P,P2):
        """ Evaluation of the 2-isogeny on a point P.
        
        Parameters
        ----------
        arg1 : self
        Self.
        
        arg2 : Point
        projective Montgomery x-coordinates P = (X:Z), where x=X/Z
        
        arg3 : Point
        projective Montgomery x-coordinates P2 that has an order of 2. 
        """
        #INIT
        t0 = fp2()
        t1 = fp2()
        t2 = fp2()
        t3 = fp2()
        tc = fp2()
        xf = fp2()
        zf = fp2()
        #Computation
        t0.add(P2.X,P2.Z)#t0 = X(P2) + Z(P2)
        t1.sous(P2.X,P2.Z)#t1 = X(P2) - Z(P2)
        t2.add(P.X,P.Z)#t2 = X(P) + Z(P)
        t3.sous(P.X,P.Z)#t3 = X(P) - Z(P)
        tc.mul(t0,t3)#t0 = t0 * t3
        t0.cpy(tc)
        tc.mul(t1,t2)#t1 = t1 * t2
        t1.cpy(tc)
        t2.add(t0,t1)
        t3.sous(t0,t1)
        xf.mul(P.X,t2)#X = X(P) * t2
        zf.mul(P.Z,t3)#Z = Z(P) * t3
        self.X.cpy(xf)
        self.Z.cpy(zf)
        
    def isogeny_3_point(self,P,P3,K1,K2):
        """ Evaluation of the 3-isogeny on a point P.
        
        Parameters
        ----------
        arg1 : self
        Self.
        
        arg2 : Point
        projective Montgomery x-coordinates P = (X:Z), where x=X/Z
        
        arg3 : Point
        projective Montgomery x-coordinates P3 that has an order of 3
        
        arg4 : fp2
        1st constant from the iso curve
        
        arg5 : fp2
        2st constant from the iso curve
        """
        #INIT
        t0 = fp2()
        t1 = fp2()
        t2 = fp2()
        tc = fp2()
        xf = fp2()
        zf = fp2()
        #Computation
        t0.add(P.X,P.Z)#t0 = X(P) + Z(P)
        t1.sous(P.X,P.Z)#t1 = X(P) + Z(P)
        tc.mul(K1,t0)#t0 = K1 * t0
        t0.cpy(tc)
        tc.mul(K2,t1)#t1 = K2 * t1
        t1.cpy(tc)
        t2.add(t0,t1)#t2 = t0 + t1
        tc.sous(t1,t0)#t0 = t1 - t0
        t0.cpy(tc)
        tc.pow2(t2)#t2 = t2²
        t2.cpy(tc)
        tc.pow2(t0)#t0 = t0²
        t0.cpy(tc)
        xf.mul(P.X,t2)#xf = X(P) * t2
        zf.mul(P.Z,t0)#zf = Z(P) * t0
        self.X.cpy(xf)
        self.Z.cpy(zf)
        
    def isogeny_4_point(self,P,P4,K1,K2,K3):
        """ Evaluation of the 4-isogeny on a point P.
        
        Parameters
        ----------
        arg1 : self
        Self.
        
        arg2 : Point
        projective Montgomery x-coordinates P = (X:Z), where x=X/Z
        
        arg3 : Point
        projective Montgomery x-coordinates P4 that has an order of 4
        
        arg4 : fp2
        1st constant from the iso curve
        
        arg5 : fp2
        2st constant from the iso curve
        
        arg6 : fp2
        3st constant from the iso curve
        
        """
        #INIT
        t0 = fp2()
        t1 = fp2()
        tc = fp2()
        xP = fp2()
        zP = fp2()
        xf = fp2()
        zf = fp2()
        #Computation
        xP.cpy(P.X)
        zP.cpy(P.Z)
        t0.add(xP,zP)#t0 = X(P) + Z(P)
        t1.sous(xP,zP)#t1 = X(P) - Z(P)
        xP.mul(t0,K2)#X(P) = t0 * K2
        zP.mul(t1,K3)#Z(P) = t1 * K3
        tc.mul(t0,t1)#t0 = t0 * t1
        t0.cpy(tc)
        tc.mul(t0,K1)#t0 = t0 * K1
        t0.cpy(tc)
        t1.add(xP,zP)#t1 = X(P) + Z(P)
        tc.sous(xP,zP)#Z(P) = X(P) - Z(P)
        zP.cpy(tc)
        tc.pow2(t1)#t1 = t1²
        t1.cpy(tc)
        tc.pow2(zP)#Z(P) = Z(P)²
        zP.cpy(tc)
        xP.add(t0,t1)#X(P) = t0 + t1
        tc.sous(zP,t0)#t0 = Z(P) - t0
        t0.cpy(tc)
        xf.mul(xP,t1)#xf = X(P) * t1
        zf.mul(zP,t0)#zf = Z(P) * t0
        self.X.cpy(xf)
        self.Z.cpy(zf)
        
    def __str__(self):
        return "(%s : %s)"%(self.X,self.Z)
    def cpy(self,P):
        self.X.cpy(P.X)
        self.Z.cpy(P.Z)



def DBLADD(P,Q,QP,para):
        """
        Simultaneous doubling and differential addition.

        Parameters
        ----------
        P : Point
        projective Montgomery points P=(XP:ZP) such that xP=XP/ZP
        Q : Point
        projective Montgomery points Q=(XQ:ZQ) such that xQ=XQ/ZQ
        QP : Point
        affine difference point : QP = Q - P
        para : Point
        Montgomery curve constant : (A+24 : C24) = (A + 2C : 4C)

        Returns
        -------
        PP : Point
        projective Montgomery points PP <- 2*P = (X2P:Z2P) such that x(2P)=X2P/Z2P.
        PQ : Point
        projective Montgomery points PQ <- P+Q = (XPQ:ZPQ) such that = x(P+Q)=XPQ/ZPQ

        """
        #INIT
        PP = point()
        PQ = point()
        t0 = fp2()
        tb = fp2()
        t1 = fp2()
        t2 =fp2()
        #double and add
        t0.add(P.X,P.Z)#t0 = Xp + Zp
        t1.sous(P.X,P.Z)#t1 = Xp - Zp
        PP.X.pow2(t0)#X2P = (Xp + Zp)²
        t2.sous(Q.X,Q.Z)#t2 = Xq - Zq
        PQ.X.add(Q.X,Q.Z)#XPQ = Xp + Xq
        tb.mul(t0,t2)#t0 = t0 * t2
        t0.cpy(tb)
        PP.Z.pow2(t1)#Z2P = t1² = (Xp - Zp)²
        tb.mul(t1,PQ.X)#t1 = t1 * XPQ = (Xp - Zp)² * (Xp + Xq)
        t1.cpy(tb)
        t2.sous(PP.X,PP.Z)#t2 = X2P - Z2P
        tb.mul(PP.X,PP.Z)#tb = X2P * Z2P
        PP.X.cpy(tb)#X2P = tb
        PQ.X.mul(para.X,t2)#XPQ = A24 * t2 = A24 * (X2P - Z2P)
        PQ.Z.sous(t0,t1)#ZPQ = t0 - t1
        tb.add(PQ.X,PP.Z)#tb = XPQ + Z2P
        PP.Z.cpy(tb)#ZPP = tb
        PQ.X.add(t0,t1)#XPQ = t0 + t1
        tb.mul(PP.Z,t2)#tb = PP.Z * t2
        PP.Z.cpy(tb)#ZPP = tb
        tb.pow2(PQ.Z)#tb = ZPQ²
        PQ.Z.cpy(tb)#ZPQ = tb
        tb.pow2(PQ.X)#tb = XPQ²
        PQ.X.cpy(tb)#XPQ = tb
        tb.mul(QP.X,PQ.Z)#tb = XQP * ZPQ
        PQ.Z.cpy(tb)#ZPQ = tb
        tb.mul(QP.Z,PQ.X)#tb = ZQP * XPQ
        PQ.X.cpy(tb)#XPQ = tb
        return PP,PQ

def j_invariance(AC):
    """
    Computes the j-invariant of a Montgomery curve with projective constant.

    Parameters
    ----------
    AC : point
    Curve parameters in projective constant.

    Returns j=256*(A^2-3*C^2)^3/(C^4*(A^2-4*C^2)),
    which is the j-invariant of the Montgomery curve B*y^2=x^3+(A/C)*x^2+x 
    or (equivalently) j-invariant of B'*y^2=C*x^3+A*x^2+C*x.
    -------
    j : fp2
    j-invariant of the Montgomery curve.
    """
    #INIT
    j = fp2()
    t0 = fp2()
    t1 = fp2()
    tcopy = fp2()
    #Computation
    j.pow2(AC.X)# j = A²
    t1.pow2(AC.Z)# t1 = C²
    t0.add(t1,t1)#t0 = t1 + t1
    tcopy.sous(j,t0)# t0 = j - t0
    t0.cpy(tcopy)
    tcopy.sous(t0,t1)# t0 = t0 - t1
    t0.cpy(tcopy)
    j.sous(t0,t1)#j = t0 - t1
    tcopy.pow2(t1)#t1 = t1²
    t1.cpy(tcopy)
    tcopy.mul(j,t1)#j = j * t1
    j.cpy(tcopy)
    tcopy.add(t0,t0)#t0 = t0 + t0
    t0.cpy(tcopy)
    tcopy.add(t0,t0)#t0 = t0 + t0
    t0.cpy(tcopy)
    t1.pow2(t0)#t1 = t0²
    tcopy.mul(t0,t1)#t0 = t0 * t1
    t0.cpy(tcopy)
    tcopy.add(t0,t0)#t0 = t0 + t0
    t0.cpy(tcopy)
    tcopy.add(t0,t0)#t0 = t0 + t0
    t0.cpy(tcopy)
    tcopy.inv(j)#j = 1 / j
    j.cpy(tcopy)
    tcopy.mul(t0,j)#j = t0 * j
    j.cpy(tcopy)
    return j

def iso_2_curve(P2):
    """
    Computes the corresponding 2-isogeny of a projective Montgomery
    point (X2:Z2) of order 2.

    Parameters
    ----------
    P2 : Point
    projective point of order two P2 = (X2:Z2).

    Returns
    -------
    AC : Point
    the 2-isogenous Montgomery curve with projective coefficients A/C

    """
    AC = point()
    A24 = fp2()
    C24 = fp2()
    cpy = fp2()
    A24.pow2(P2.X)#A24 = X(P2)²
    C24.pow2(P2.Z)#Z24 = Z(P2)²
    cpy.sous(C24,A24)#A24 = C24 - A24
    A24.cpy(cpy)
    AC.X.cpy(A24)
    AC.Z.cpy(C24)    
    return AC

def iso_3_curve(P3):
    """
    Computes the corresponding 2-isogeny of a projective Montgomery
    point (X3:Z3) of order 3.

    Parameters
    ----------
    P3 : Point
    projective point of order two P2 = (X3:Z3).

    Returns
    -------
    AC : Point
    the 3-isogenous Montgomery curve with projective coefficients A/C
    
    K1 : fp2
    1st constant from the iso curve
        
    K2 : fp2
    2st constant from the iso curve

    """
    #INIT
    AC = point()
    K1 = fp2()
    K2 = fp2()
    Aup = fp2()
    Adown = fp2()
    cpy = fp2()
    t0 = fp2()
    t1 = fp2()
    t2 = fp2()
    t3 = fp2()
    t4 = fp2()
    #Computation
    K1.sous(P3.X, P3.Z)#K1 = X(P3) - Z(P3)
    t0.pow2(K1)#t0 = K1²
    K2.add(P3.X, P3.Z)#K2 = X(P3) + Z(P3)
    t1.pow2(K2)#t1 = K2²
    t2.add(t0,t1)#t2 = t1 + t0
    t3.add(K1,K2)#t3 = K1 + K2
    cpy.pow2(t3)#t3 = t3²
    t3.cpy(cpy)
    cpy.sous(t3,t2)#t3 = t3 - t2
    t3.cpy(cpy)
    t2.add(t1,t3)#t2 = t1 + t3
    cpy.add(t3,t0)#t3 = t3 + t0
    t3.cpy(cpy)
    t4.add(t3,t0)#t4 = t3 + t0
    cpy.add(t4,t4)#t4 = t4 + t4
    t4.cpy(cpy)
    cpy.add(t1,t4)#t4 = t1 + t4
    t4.cpy(cpy)
    Adown.mul(t2,t4)#A24- = t2 * t4
    t4.add(t1,t2)#t4 = t1 + t2
    cpy.add(t4,t4)#t4 = t4 + t4
    t4.cpy(cpy)
    cpy.add(t0,t4)#t4 = t0 + t4
    t4.cpy(cpy)
    Aup.mul(t3,t4)#A24+ = t3 * t4
    AC.X.cpy(Aup)
    AC.Z.cpy(Adown)    
    return AC,K1,K2

def iso_4_curve(P4):
    """
    Computes the corresponding 2-isogeny of a projective Montgomery
    point (X4:Z4) of order 4.

    Parameters
    ----------
    P4 : Point
    projective point of order two P4 = (X4:Z4).

    Returns
    -------
    AC : Point
    the 4-isogenous Montgomery curve with projective coefficients A/C
    
    K1 : fp2
    1st constant from the iso curve
        
    K2 : fp2
    2st constant from the iso curve
    
    K3 : fp2
    2st constant from the iso curve

    """
    #INIT
    AC = point()
    K1 = fp2()
    K2 = fp2()
    K3 = fp2()
    A24 = fp2()
    C24 = fp2()
    cpy = fp2()
    #Computation
    K2.sous(P4.X, P4.Z)#K2 = X(P4) - Z(P4)
    K3.add(P4.X, P4.Z)#K1 = X(P4) + Z(P4)
    K1.pow2(P4.Z)#K1 = Z(P4)²
    cpy.add(K1,K1)#K1 = K1 + K1
    K1.cpy(cpy)
    C24.pow2(K1)#C24 = K1²
    cpy.add(K1,K1)#K1 = K1 + K1
    K1.cpy(cpy)
    A24.pow2(P4.X)#A24 = X(P4)²
    cpy.add(A24,A24)#A24 = A24 + A24
    A24.cpy(cpy)
    cpy.pow2(A24)#A24 = A24²
    A24.cpy(cpy)
    AC.X.cpy(A24)
    AC.Z.cpy(C24)  
    return AC,K1,K2,K3







def e_2_iso(S, para):
    Z1 = point()
    Z2 = point()
    Z3 = point()
    K1 = fp2()
    K2 = fp2()
    K3 = fp2()
    T = point()
    A24 = point()
    A24.cpy(para)
    e = e1 - 2
    while(e >= 0):
        T.puissance2_k(e, S, A24)
        (A24, K1, K2, K3) = iso_4_curve(T)

        if e != 0:
            pt = point()
            pt.cpy(S)
            S.isogeny_4_point(pt, pt, K1, K2, K3)
        Z1.isogeny_4_point(Z1, Z1, K1, K2, K3)
        Z2.isogeny_4_point(Z2, Z2, K1, K2, K3)
        Z3.isogeny_4_point(Z3, Z3, K1, K2, K3)
        e = e - 2
    return A24, Z1, Z2, Z3

def e_2_iso_strat(S, para):
    Z1 = point()
    Z2 = point()
    Z3 = point()
    K1 = fp2()
    K2 = fp2()
    K3 = fp2()
    A24 = point()
    A24.cpy(para)
    stat = []
    if e1 == 0xD8:
        stat = [48, 28, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 13,
                 7, 4, 2, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 4, 2, 1, 1, 2, 1, 1, 2, 1, 1, 1, 21, 12, 7, 4, 2, 1, 1, 2, 1,
                 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2,
                 1, 1]
    elif e1 == 0xFA:
        stat = [61, 32, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1,
                16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 29, 16,
                8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 13, 8, 4, 2,
                1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 5, 4, 2, 1, 1, 2, 1, 1, 2, 1, 1, 1]
    elif e1 == 0x131:
        stat = [67, 37, 21, 12, 7, 4, 2, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 9, 5, 3, 2, 1,
                1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 16, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1,
                1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 33, 16, 8, 5, 2, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 2, 1,
                1 , 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 16, 8, 4, 2, 1, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8,
                4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1]
    elif e1 == 0x174:
        stat = [80, 48, 27, 15, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 7, 4, 2, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1,
                12, 7, 4, 2, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 21, 12, 7, 4, 2, 1, 1,
                2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1,
                1, 1, 2, 1, 1, 33, 20, 12, 7, 4, 2, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1,
                8, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 2, 1, 1, 16, 8, 4, 2, 1, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1,
                1, 8 , 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1]
    else:
        print("false val in e2")
        exit()
    queue = deque()
    
    pt = point()
    pt.cpy(S)
    elem = []
    elem.append(e1 // 2)
    elem.append(pt)
    queue.append(elem)
    i = 0
    while len(queue) != 0:
        element = queue.pop()
        #print(f'A24 {A24.X.a0}')
        #print(element)
        if element[0] == 1:
            (A24, K1, K2, K3) = iso_4_curve(element[1])
            new_q = deque()
            while len(queue) != 0:
                update = queue.popleft()
                updated_elem = []
                updated_elem.append(update[0] - 1)
                temp = point()
                temp.isogeny_4_point(update[1], update[1], K1, K2, K3)
                updated_elem.append(temp)
                new_q.append(updated_elem)
            queue = new_q
            Z1.isogeny_4_point(Z1, Z1, K1, K2, K3)
            Z2.isogeny_4_point(Z2, Z2, K1, K2, K3)
            Z3.isogeny_4_point(Z3, Z3, K1, K2, K3)

        elif 0 < stat[i] and stat[i] < element[0]:
            queue.append(element)
            temp = point()
            temp.puissance2_k(2*stat[i], element[1], A24)
            new_elem = []
            new_elem.append(element[0] - stat[i])
            new_elem.append(temp)
            queue.append(new_elem)
            i = i + 1
        else:
            #print(element[0])
            #print(stat[i])
            #print(i)
            #print(0x89)
            assert False, "false strat"
    #print(f'A24 {A24.X.a0}')
    return A24, Z1, Z2, Z3

def e_2_iso_strat_x(S, para, x1, x2, x3):
    Z1 = point()
    Z2 = point()
    Z3 = point()
    Z1.cpy(x1)
    Z2.cpy(x2)
    Z3.cpy(x3)


    K1 = fp2()
    K2 = fp2()
    K3 = fp2()
    A24 = point()
    A24.cpy(para)
    stat = []
    if e1 == 0xD8:
        stat = [48, 28, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 13,
                 7, 4, 2, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 4, 2, 1, 1, 2, 1, 1, 2, 1, 1, 1, 21, 12, 7, 4, 2, 1, 1, 2, 1,
                 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2,
                 1, 1]
    elif e1 == 0xFA:
        stat = [61, 32, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1,
                16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 29, 16,
                8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 13, 8, 4, 2,
                1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 5, 4, 2, 1, 1, 2, 1, 1, 2, 1, 1, 1]
    elif e1 == 0x131:
        stat = [67, 37, 21, 12, 7, 4, 2, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 9, 5, 3, 2, 1,
                1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 16, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1,
                1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 33, 16, 8, 5, 2, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 2, 1,
                1 , 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 16, 8, 4, 2, 1, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8,
                4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1]
    elif e1 == 0x174:
        stat = [80, 48, 27, 15, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 7, 4, 2, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1,
                12, 7, 4, 2, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 21, 12, 7, 4, 2, 1, 1,
                2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1,
                1, 1, 2, 1, 1, 33, 20, 12, 7, 4, 2, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1,
                8, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 2, 1, 1, 16, 8, 4, 2, 1, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1,
                1, 8 , 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1]
    else:
        print("false val in e2")
        exit()
    #print(f'strat used:\n{stat}')
    queue = deque()
    
    pt = point()
    pt.cpy(S)
    elem = []
    elem.append((e1 // 2))
    elem.append(pt)
    queue.append(elem)
    i = 0
    while len(queue) != 0:
        element = queue.pop()
        #print(f'A24 {A24.X.a0}')
        #print(element)
        if element[0] == 1:
            (A24, K1, K2, K3) = iso_4_curve(element[1])
            new_q = deque()
            while len(queue) != 0:
                update = queue.popleft()
                updated_elem = []
                updated_elem.append(update[0] - 1)
                temp = point()
                temp.isogeny_4_point(update[1], update[1], K1, K2, K3)
                updated_elem.append(temp)
                new_q.append(updated_elem)
            queue = new_q
            pt = point()
            pt.cpy(Z1)
            Z1.isogeny_4_point(pt, pt, K1, K2, K3)
            pt.cpy(Z2)
            Z2.isogeny_4_point(pt, pt, K1, K2, K3)
            pt.cpy(Z3)
            Z3.isogeny_4_point(pt, pt, K1, K2, K3)

        elif 0 < stat[i] and stat[i] < element[0]:
            queue.append(element)
            temp = point()
            temp.puissance2_k(2*stat[i], element[1], A24)
            new_elem = []
            new_elem.append(element[0] - stat[i])
            new_elem.append(temp)
            queue.append(new_elem)
            i = i + 1
        else:
            #print(element[0])
            #print(stat[i])
            #print(i)
            #print(0x89)
            assert False, "false strat"
    #print(f'A24 {A24.X.a0}')
    return A24, Z1, Z2, Z3

def e_2_iso_x(S, para, x1, x2, x3):
    Z1 = point()
    Z2 = point()
    Z3 = point()
    Z1.cpy(x1)
    Z2.cpy(x2)
    Z3.cpy(x3)
    K1 = fp2()
    K2 = fp2()
    K3 = fp2()
    T = point()
    A24 = point()
    A24.cpy(para)
    e = e1 - 2
    while(e >= 0):
        T.puissance2_k(e, S, A24)
        (A24, K1, K2, K3) = iso_4_curve(T)

        if e != 0:
            pt = point()
            pt.cpy(S)
            S.isogeny_4_point(pt, pt, K1, K2, K3)
        pt = point()
        pt.cpy(Z1)
        Z1.isogeny_4_point(pt, pt, K1, K2, K3)
        pt.cpy(Z2)
        Z2.isogeny_4_point(pt, pt, K1, K2, K3)
        pt.cpy(Z3)
        Z3.isogeny_4_point(pt, pt, K1, K2, K3)
        e = e - 2
    return A24, Z1, Z2, Z3

def e_2_iso_my_strat_x(S, para, x1, x2, x3):
    #strat = [48, 28, 16, 8, 4, 2, 1, 0, 0, 1, 0, 0, 2, 1, 0, 0, 1, 0, 0, 4, 2, 1, 0, 0, 1, 0, 0, 2, 1, 0, 0, 1, 0, 0, 8, 4, 2, 1, 0, 0, 
    #         1, 0, 0, 2, 1, 0, 0, 1, 0, 0, 4, 2, 1, 0, 0, 1, 0, 0, 2, 1, 0, 0, 1, 0, 0, 13, 7, 4, 2, 1, 0, 0, 1, 0, 0, 2, 1, 0, 0, 1, 0, 
    #         0, 3, 2, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 5, 4, 2, 1, 0, 0, 1, 0, 0, 2, 1, 0, 0, 1, 0, 0, 2, 1, 1, 0, 0, 0, 1, 0, 0, 21, 
    #         12, 7, 4, 2, 1, 0, 0, 1, 0, 0, 2, 1, 0, 0, 1, 0, 0, 3, 2, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 5, 3, 2, 1, 0, 0, 1, 0, 0, 1, 1, 
    #         0, 0, 0, 2, 1, 1, 0, 0, 0, 1, 0, 0, 9, 5, 3, 2, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 2, 1, 1, 0, 0, 0, 1, 0, 0, 4, 2, 1, 1, 
    #         0, 0, 0, 1, 0, 0, 2, 1, 0, 0, 1, 0, 0]
    strat = [67, 37, 21, 12,  7, 4, 2, 1, 0, 0, 1, 0, 0, 2, 1, 0, 0, 1, 0, 0, 3, 2, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0,  5, 3, 2, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 2, 1, 1, 0, 0, 0, 1, 0, 0,  9,  5, 3, 2, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 2, 1, 1, 0, 0, 0, 1, 0, 0, 4, 2, 1, 1, 0, 0, 0, 1, 0, 0, 2, 1, 0, 0, 1, 0, 0, 16,  9,  5, 3, 2, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 2, 1, 1, 0, 0, 0, 1, 0, 0, 4, 2, 1, 1, 0, 0, 0, 1, 0, 0, 2, 1, 0, 0, 1, 0, 0,  8, 4, 2, 1, 0, 0, 1, 0, 0, 2, 1, 0, 0, 1, 0, 0, 4, 2, 1, 0, 0, 1, 0, 0, 2, 1, 0, 0, 1, 0, 0, 33, 16,  8,  5, 2, 1, 1, 0, 0, 0, 1, 0, 0, 2, 1, 1, 0, 0, 0, 1, 0, 0, 4, 2, 1, 0, 0, 1, 0, 0, 2, 1, 0, 0, 1, 0, 0,  8, 4, 2, 1, 0, 0, 1, 0, 0, 2, 1, 0, 0, 1, 0, 0, 4, 2, 1, 0, 0, 1, 0, 0, 2, 1, 0, 0, 1, 0, 0, 16,  8, 4, 2, 1, 1, 0, 0, 0, 1, 0, 0, 2, 1, 0, 0, 1, 0, 0, 4, 2, 1, 0, 0, 1, 0, 0, 2, 1, 0, 0, 1, 0, 0,  8, 4, 2, 1, 0, 0, 1, 0, 0, 2, 1, 0, 0, 1, 0, 0, 4, 2, 1, 0, 0, 1, 0, 0, 2, 1, 0, 0, 1, 0, 0]
    length = 9
    strat_len = 303
    last_pos = 0
    queue = []
    queue.append(S)
    for i in range(length):
        queue.append(point())
    i = 0
    A24 = point()
    A24.cpy(para)
    Opt1 = point()
    Opt2 = point()
    Opt3 = point()
    Opt1.cpy(x1)
    Opt2.cpy(x2)
    Opt3.cpy(x3)
    #while(i < strat_len):
    for i in range(strat_len):
        K1 = fp2()
        K2 = fp2()
        K3 = fp2()
        curr_val = strat[i]
        if curr_val == 0:    
            (A24, K1, K2, K3) = iso_4_curve(queue[last_pos])
            
            for x in range(last_pos):
                up_pt = point()
                up_pt.cpy(queue[x])
                queue[x].isogeny_4_point(up_pt,up_pt ,K1,K2,K3)
            last_pos = last_pos - 1
            temp = point()
            temp.cpy(Opt1)
            Opt1.isogeny_4_point(temp,temp ,K1,K2,K3)
            temp.cpy(Opt2)
            Opt2.isogeny_4_point(temp,temp ,K1,K2,K3)
            temp.cpy(Opt3)
            Opt3.isogeny_4_point(temp,temp ,K1,K2,K3)
            
        else:
            new_pt = point()
            #print(queue[last_pos])
            #print(A24)
            new_pt.puissance2_k(2*strat[i], queue[last_pos], A24)
            last_pos = last_pos+1
            queue[last_pos] = new_pt

            #print(new_pt)
            #print(strat[i])
            #print("\n\n")
        #i = i + 1
    return A24, Opt1, Opt2, Opt3

def e_3_iso_strat(S, para):
    Z1 = point()
    Z2 = point()
    Z3 = point()
    K1 = fp2()
    K2 = fp2()
    A24 = point()
    A24.cpy(para)
    stat = []
    if e2 == 0x89:
        stat = [66, 33, 17, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 1, 2, 1, 1, 4, 2, 1, 1,
                2, 1, 1, 16, 8, 4, 2, 1, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 32,
                16, 8, 4, 3, 1, 1, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 16, 8, 4,
                2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1]
    elif e2 == 0x9F:
        stat = [71, 38, 21, 13, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 5, 4, 2, 1, 1, 2, 1, 1, 2, 1, 1, 1,
                9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 17, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1,
                4, 2, 1, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 33, 17, 9, 5, 3, 2, 1, 1,
                1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 16, 8, 4,
                2, 1, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1]
    elif e2 == 0xC0:
        stat = [86, 48, 27, 15, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 7, 4, 2, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 12,
                7, 4, 2, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 21, 12, 7, 4, 2, 1, 1, 2, 1, 1,
                3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1,
                1, 38, 21, 12, 7, 4, 2, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 9, 5, 3, 2, 1,
                1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 17, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1,
                8, 4, 2, 1, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1]
    elif e2 == 0xEF:
        stat = [112, 63, 32, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1,
                1, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 31, 16,
                8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 15, 8, 4, 2, 1,
                1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 7, 4, 2, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 49, 31, 16, 8, 4, 2, 1, 1, 2,
                1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 15, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2,
                1, 1, 2, 1, 1, 7, 4, 2, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 21, 12, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1,
                1, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1]
    else:
        print("false val in e3")
        exit()


    queue = deque()
    
    pt = point()
    pt.cpy(S)
    elem = []
    elem.append(e2)
    elem.append(pt)
    queue.append(elem)
    i = 0
    while len(queue) != 0:
        element = queue.pop()
        #print(f'A24 {A24.X.a0}')
        #print(element)
        if element[0] == 1:
            (A24, K1, K2) = iso_3_curve(element[1])
            new_q = deque()
            while len(queue) != 0:
                update = queue.popleft()
                updated_elem = []
                updated_elem.append(update[0] - 1)
                temp = point()
                temp.isogeny_3_point(update[1], update[1], K1, K2)
                updated_elem.append(temp)
                new_q.append(updated_elem)
            queue = new_q
            Z1.isogeny_3_point(Z1, Z1, K1, K2)
            Z2.isogeny_3_point(Z2, Z2, K1, K2)
            Z3.isogeny_3_point(Z3, Z3, K1, K2)

        elif 0 < stat[i] and stat[i] < element[0]:
            queue.append(element)
            temp = point()
            temp.TPLe(stat[i], element[1], A24)
            new_elem = []
            new_elem.append(element[0] - stat[i])
            new_elem.append(temp)
            queue.append(new_elem)
            i = i + 1
        else:
            #print(element[0])
            #print(stat[i])
            #print(i)
            #print(0x89)
            assert False, "false strat"
    #print(f'A24 {A24.X.a0}')
    return A24, Z1, Z2, Z3

def e_3_iso_strat_x(S, para, x1, x2, x3):
    Z1 = point()
    Z2 = point()
    Z3 = point()
    Z1.cpy(x1)
    Z2.cpy(x2)
    Z3.cpy(x3)

    K1 = fp2()
    K2 = fp2()
    A24 = point()
    A24.cpy(para)
    stat = []
    if e2 == 0x89:
        stat = [66, 33, 17, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 1, 2, 1, 1, 4, 2, 1, 1,
                2, 1, 1, 16, 8, 4, 2, 1, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 32,
                16, 8, 4, 3, 1, 1, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 16, 8, 4,
                2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1]
    elif e2 == 0x9F:
        stat = [71, 38, 21, 13, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 5, 4, 2, 1, 1, 2, 1, 1, 2, 1, 1, 1,
                9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 17, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1,
                4, 2, 1, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 33, 17, 9, 5, 3, 2, 1, 1,
                1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 16, 8, 4,
                2, 1, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1]
    elif e2 == 0xC0:
        stat = [86, 48, 27, 15, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 7, 4, 2, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 12,
                7, 4, 2, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 21, 12, 7, 4, 2, 1, 1, 2, 1, 1,
                3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1,
                1, 38, 21, 12, 7, 4, 2, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 9, 5, 3, 2, 1,
                1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 17, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1,
                8, 4, 2, 1, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1]
    elif e2 == 0xEF:
        stat = [112, 63, 32, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1,
                1, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 31, 16,
                8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 15, 8, 4, 2, 1,
                1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 7, 4, 2, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 49, 31, 16, 8, 4, 2, 1, 1, 2,
                1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 15, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2,
                1, 1, 2, 1, 1, 7, 4, 2, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 21, 12, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1,
                1, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1]
    else:
        print("false val in e3")
        exit()


    queue = deque()
    
    pt = point()
    pt.cpy(S)
    elem = []
    elem.append(e2)
    elem.append(pt)
    queue.append(elem)
    i = 0
    while len(queue) != 0:
        element = queue.pop()
        #print(f'A24 {A24.X.a0}')
        #print(element)
        if element[0] == 1:
            (A24, K1, K2) = iso_3_curve(element[1])
            new_q = deque()
            while len(queue) != 0:
                update = queue.popleft()
                updated_elem = []
                updated_elem.append(update[0] - 1)
                temp = point()
                temp.isogeny_3_point(update[1], update[1], K1, K2)
                updated_elem.append(temp)
                new_q.append(updated_elem)
            queue = new_q
            pt = point()
            pt.cpy(Z1)
            Z1.isogeny_3_point(pt, pt, K1, K2)
            pt.cpy(Z2)
            Z2.isogeny_3_point(pt, pt, K1, K2)
            pt.cpy(Z3)
            Z3.isogeny_3_point(pt, pt, K1, K2)

        elif 0 < stat[i] and stat[i] < element[0]:
            queue.append(element)
            temp = point()
            temp.TPLe(stat[i], element[1], A24)
            new_elem = []
            new_elem.append(element[0] - stat[i])
            new_elem.append(temp)
            queue.append(new_elem)
            i = i + 1
        else:
            #print(element[0])
            #print(stat[i])
            #print(i)
            #print(0x89)
            assert False, "false strat"
    #print(f'A24 {A24.X.a0}')
    return A24, Z1, Z2, Z3

def e_3_iso_x(S, para, x1, x2, x3):
    Z1 = point()
    Z2 = point()
    Z3 = point()
    Z1.cpy(x1)
    Z2.cpy(x2)
    Z3.cpy(x3)
    K1 = fp2()
    K2 = fp2()
    T = point()
    A24 = point()
    A24.cpy(para)
    e = e2 - 1
    while(e >= 0):
        T.TPLe(e, S, A24)
        (A24, K1, K2) = iso_3_curve(T)

        if e != 0:
            pt = point()
            pt.cpy(S)
            S.isogeny_3_point(pt, pt, K1, K2)
        Z1.isogeny_3_point(Z1, Z1, K1, K2)
        Z2.isogeny_3_point(Z2, Z2, K1, K2)
        Z3.isogeny_3_point(Z3, Z3, K1, K2)
        e = e - 1
    return A24, Z1, Z2, Z3

def change_params(e1_new, e2_new):
    global e1
    global e2
    global p
    e1 = e1_new
    e2 = e2_new
    p = (2**e1_new)*(3**e2_new)-1

def get_A(xp, xq, xpq):
    t0 = fp2()
    t1 = fp2()
    A = fp2()
    c = fp2()
    t1.add(xp, xq)  #1
    t0.mul(xp, xq)  #2
    A.mul(xpq, t1)  #3 
    c.cpy(A)
    A.add(c, t0)    #4
    c.cpy(t0)
    t0.mul(c, xpq)  #5
    c.cpy(A)
    one = fp2(1, 0)
    A.sous(A, one)   #6
    c.cpy(t0)
    t0.add(c, c)    #7
    c.cpy(t1)
    t1.add(c, xpq) #8
    c.cpy(t0)
    t0.add(c, c)    #9
    c.cpy(A)
    A.pow2(c)       #10
    c.cpy(t0)
    t0.inv(c)       #11
    c.cpy(A)
    A.mul(c, t0)        #12
    c.cpy(A)
    A.sous(c, t1)
    return A

def isogen_2(secred_key, xP2 = fp2(), xQ2 = fp2(), xR2 = fp2(), xP3 = fp2(), xQ3 = fp2(), xR3 = fp2()):
    A = point(fp2(6, 0), fp2(1,0))
    A_plus = point(fp2(8, 0), fp2(4,0))

    X1 = point(xP3, fp2(1, 0))
    X2 = point(xQ3, fp2(1, 0))
    X3 = point(xR3, fp2(1, 0))
    Xs = point()

    P  = point(xP2, fp2(1, 0))
    Q  = point(xQ2, fp2(1, 0))
    PQ = point(xR2, fp2(1, 0))
    Xs.ladder3pt(secred_key, P, Q, PQ, A)
    print('isogen_2 ladder erg:')
    print(f'{hex(Xs.X.a0)}\n{hex(Xs.X.a1)}\n{hex(Xs.Z.a0)}\n{hex(Xs.Z.a1)}')
    A_p_new = point()
    X1_n = point()
    X2_n = point()
    X3_n = point()

    A_p_new2 = point()
    X1_n2 = point()
    X2_n2 = point()
    X3_n2 = point()
    

    
    (A_p_new, X1_n, X2_n, X3_n) = e_2_iso_x(Xs, A_plus, X1, X2, X3)
    #(A_p_new2, X1_n2, X2_n2, X3_n2) = e_2_iso_strat_x(Xs, A_plus, X1, X2, X3)

    print(f'{X1_n == X1_n2}')#, "X1 erg wrong"
    print(f'{X2_n == X2_n2}')#, "X1 erg wrong"
    print(f'{X3_n == X3_n2}')#, "X1 erg wrong"
    print("e_2_iso_erg:")
    temp = fp2()
    x1_res = fp2()
    temp.inv(X1_n.Z)
    x1_res.mul(X1_n.X, temp)

    x2_res = fp2()
    temp.inv(X2_n.Z)
    x2_res.mul(X2_n.X, temp)

    x3_res = fp2()
    temp.inv(X3_n.Z)
    x3_res.mul(X3_n.X, temp)
    print(f'x1: {hex(x1_res.a0)}\n{hex(x1_res.a1)}')
    print(f'x2: {hex(x2_res.a0)}\n{hex(x2_res.a1)}')
    print(f'x3: {hex(x3_res.a0)}\n{hex(x3_res.a1)}')
    return (x1_res, x2_res, x3_res)

def isogen_3(secred_key, xP2 = fp2(), xQ2 = fp2(), xR2 = fp2(), xP3 = fp2(), xQ3 = fp2(), xR3 = fp2()):
    A = point(fp2(6, 0), fp2(1,0))
    A_plus = point(fp2(8, 0), fp2(4,0))

    X1 = point(xP2, fp2(1, 0))
    X2 = point(xQ2, fp2(1, 0))
    X3 = point(xR2, fp2(1, 0))
    Xs = point()

    P  = point(xP3, fp2(1, 0))
    Q  = point(xQ3, fp2(1, 0))
    PQ = point(xR3, fp2(1, 0))
    Xs.ladder3pt(secred_key, P, Q, PQ, A)
    print('isogen_3 ladder erg:')
    print(f'{hex(Xs.X.a0)}\n{hex(Xs.X.a1)}\n{hex(Xs.Z.a0)}\n{hex(Xs.Z.a1)}')
    A_p_new = point()
    X1_n = point()
    X2_n = point()
    X3_n = point()

    

    
    (A_p_new, X1_n, X2_n, X3_n) = e_3_iso_x(Xs, A_plus, X1, X2, X3)
    
    
    temp = fp2()
    x1_res = fp2()
    temp.inv(X1_n.Z)
    x1_res.mul(X1_n.X, temp)

    x2_res = fp2()
    temp.inv(X2_n.Z)
    x2_res.mul(X2_n.X, temp)

    x3_res = fp2()
    temp.inv(X3_n.Z)
    x3_res.mul(X3_n.X, temp)
    
    return (x1_res, x2_res, x3_res)

def isoex_2(x1, x2, x3, secret_key):
    A = point(get_A(x1, x2, x3), fp2(1,0))
    P  = point(x1, fp2(1, 0))
    Q  = point(x2, fp2(1, 0))
    PQ = point(x3, fp2(1, 0))
    Xs = point()
    Xs.ladder3pt(secret_key, P, Q, PQ, A)
    A_plus = point(fp2(), fp2(4,0))
    A_plus.X.add(A.X, fp2(2,0))

    A_plus_n = point()
    X1 = point()
    X2 = point()
    X3 = point()
    X1_n = point()
    X2_n = point()
    X3_n = point()

    (A_plus_n, X1_n, X2_n, X3_n) = e_2_iso_strat_x(Xs, A_plus, X1, X2, X3)


    temp = fp2()
    temp2 = fp2()

    temp.mul(A_plus_n.X, fp2(4,0))
    temp2.mul(A_plus_n.Z, fp2(2,0))
    A.X.sous(temp, temp2)
    A.Z.cpy(A_plus_n.Z)


    jinv = fp2()

    jinv = j_invariance(A)
    return jinv
    
def isoex_3(x1, x2, x3, secret_key):
    A = point(get_A(x1, x2, x3), fp2(1,0))
    P  = point(x1, fp2(1, 0))
    Q  = point(x2, fp2(1, 0))
    PQ = point(x3, fp2(1, 0))
    Xs = point()
    Xs.ladder3pt(secret_key, P, Q, PQ, A)
    A_plus = point(fp2(), fp2())
    A_plus.X.add(A.X, fp2(2,0))
    A_plus.Z.sous(A.X, fp2(2,0))

    A_plus_n = point()
    X1 = point()
    X2 = point()
    X3 = point()
    X1_n = point()
    X2_n = point()
    X3_n = point()

    (A_plus_n, X1_n, X2_n, X3_n) = e_3_iso_x(Xs, A_plus, X1, X2, X3)


    temp = fp2()
    #temp2 = fp2()

    temp.add(A_plus_n.X, A_plus_n.Z)
    #temp2.mul(A_plus_n.Z, fp2(2,0))
    A.X.mul(temp, fp2(2,0))
    A.Z.sous(A_plus_n.X, A_plus_n.Z)


    jinv = fp2()

    jinv = j_invariance(A)
    return jinv
    