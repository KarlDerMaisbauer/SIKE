

from sidh_434 import*
import sympy


bits = 9*64

def bit_not(n, numbits):
    return (1 << numbits) - 1 - n

path = "erg.txt"

def tohex(val, nbits):
  return hex((val + (1 << nbits)) % (1 << nbits))

def get_val(x):
    x_str = x.replace("\n", "")
    x_int = int(x_str, 0)
    l = list(x_str)
    if int(l[2],16) > 7:
        return - (bit_not(x_int, x_int.bit_length()) + 1)
    else:
        return x_int

def fp_add(x, y, res):
    add1 = get_val(x)
    add2 = get_val(y)
    res_int = get_val(res)
    assert add1 + add2 == res_int, f'fp_add not working properly\n {hex(add1)} +\n {hex(add2)} should be\n {hex(add1 + add2)} but is \n {hex(res_int)}'

def fp_sub(x, y, res):
    add1 = get_val(x)
    add2 = get_val(y)
    res_int = get_val(res)
    assert add1 - add2 == res_int, f'fp_sub not working properly\n {hex(add1)} -\n {hex(add2)} should be\n {hex(add1 - add2)} but is \n {hex(res_int)}'

def fp_mul(x,y,res):
    add1 = get_val(x)
    add2 = get_val(y)
    res_int = get_val(res)
    assert add1 * add2 == res_int, f'fp_mul not working properly\n {hex(add1)} *\n {hex(add2)} should be\n {hex(add1 * add2)} but is \n {hex(res_int)}'

def fp_div(x, y, res):
    add1 = get_val(x)
    add2 = get_val(y)
    change = 1
    if add1 < 0 and add2 < 0:
        add1 = -add1
        add2 = -add2
    elif add1 < 0:
        add1 = -add1
        change = -1
    elif add2 < 0:
        add2 = -add2
        change = -1
    res_int = get_val(res)
    assert change * (add1 // add2) == res_int, f'fp_div not working properly\n {hex(add1)} /\n {hex(add2)} should be\n {hex(add1 // add2)} but is \n {hex(res_int)}'

def fp_mod(x, y, res):
    add1 = get_val(x)
    mod = get_val(y)
    res_int = get_val(res)
    assert add1 % mod == res_int, f'fp_mod not working properly\n {hex(add1)} %\n {hex(mod)} should be\n {hex(add1 % mod)} but is \n {hex(res_int)}'

def fp_mod_2k(x, y, res):
    add1 = get_val(x)
    mod = get_val(y)
    res_int = get_val(res)
    assert (add1 % mod) == res_int, f'fp_mod_2k not working properly\n {hex(add1)} %\n {hex(mod)} should be\n {hex(add1 % mod)} but is \n {hex(res_int)}'

def fp_addm(x, y, mod, res):
    add1 = get_val(x)
    add2 = get_val(y)
    m = get_val(mod)
    res_int = get_val(res)
    assert (add1 + add2) % m == res_int, f'fp_addm not working properly\n {hex(add1)} +\n {hex(add2)} should be\n {hex(add1 + add2)} but is \n {hex(res_int)}'

def fp_subm(x, y, mod, res):
    add1 = get_val(x)
    add2 = get_val(y)
    m = get_val(mod)
    res_int = get_val(res)
    assert (add1 - add2) % m == res_int % m, f'fp_subm not working properly\n {hex(add1)} -\n {hex(add2)} should be\n {hex(add1 + add2)} but is \n {hex(res_int)}'

def fp_get_len(x, res):
    add1 = get_val(x)
    r = add1
    Tab = []
    while (r!= 1):
        Tab.append(r%2)
        r=r//2
    Tab.append(1)
    res_int = len(Tab)
    res_int = int(res)
    assert add1.bit_length() == res_int, f'fp_get_len not working properly\n should be\n {add1.bit_length()} but is \n {res_int}'

def f2p_add(x, y, res):
    add1 = get_val(x)
    add2 = get_val(y)
    res_int = get_val(res)
    assert add1 + add2 == res_int, f'f2p_add not working properly\n {hex(add1)} +\n {hex(add2)} should be\n {hex(add1 + add2)} but is \n {hex(res_int)}'

def f2p_sub(x, y, res):
    add1 = get_val(x)
    add2 = get_val(y)
    res_int = get_val(res)
    assert add1 - add2 == res_int, f'f2p_sub not working properly\n {hex(add1)} -\n {hex(add2)} should be\n {hex(add1 - add2)} but is \n {hex(res_int)}'

def f2p_mul(x,y,res):
    add1 = get_val(x)
    add2 = get_val(y)
    res_int = get_val(res)
    assert add1 * add2 == res_int, f'f2p_mul not working properly\n {hex(add1)} +\n {hex(add2)} should be\n {hex(add1 * add2)} but is \n {hex(res_int)}'

def f2p_div(x, y, res):
    add1 = get_val(x)
    add2 = get_val(y)
    res_int = get_val(res)
    change = 1
    if add1 < 0 and add2 < 0:
        add1 = -add1
        add2 = -add2
    elif add1 < 0:
        add1 = -add1
        change = -1
    elif add2 < 0:
        add2 = -add2
        change = -1
    res_int = get_val(res)
    assert change * (add1 // add2) == res_int,f'f2p_div not working properly\n {hex(add1)} /\n {hex(add2)} should be\n {hex(add1 // add2)} but is \n {hex(res_int)}'

def f2p_mod(x, y, res):
    add1 = get_val(x)
    mod = get_val(y)
    res_int = get_val(res)
    assert add1 % mod == res_int, f'f2p_mod not working properly\n {hex(add1)} %\n {hex(mod)} should be\n {hex(add1 % mod)} but is \n {hex(res_int)}'

def f2p_div_2k(x, y, res):
    add1 = get_val(x)
    add2 = get_val(y)
    res_int = get_val(res)
    change = 1
    if add1 < 0 and add2 < 0:
        add1 = -add1
        add2 = -add2
    elif add1 < 0:
        add1 = -add1
        change = -1
    elif add2 < 0:
        add2 = -add2
        change = -1
    assert change * (add1 // add2) == res_int, f'f2p_div_2k not working properly\n {hex(add1)} /\n {hex(add2)} should be\n {hex(add1 // add2)} but is \n {hex(res_int)}'

def f2p_mod_2k(x, y, res):
    add1 = get_val(x)
    mod = get_val(y)
    res_int = get_val(res)
    assert add1 % mod == res_int, f'f2p_mod_2k not working properly\n {hex(add1)} %\n {hex(mod)} should be\n {hex(add1 % mod)} but is \n {hex(res_int)}'

def fp2_add(x1, x2, y1, y2, res1, res2):
    resx = get_val(res1) 
    resy = get_val(res2) 
    p1x = get_val(x1) 
    p1y = get_val(x2) 
    p2x = get_val(y1) 
    p2y = get_val(y2) 
    
    res = fp2(resx, resy)
    add1 = fp2(p1x, p1y)
    add2 = fp2(p2x, p2y)
    res = fp2(resx, resy)
    res_calc = fp2(0,0)
    res_calc.add(add1, add2)
    assert res == res_calc, f'fp2_add not working properly\n should be\n {hex(res_calc.a0)} \n {hex(res_calc.a1)} but is \n {hex(res.a0)} \n {hex(res.a1)}'

def fp2_sub(x1, x2, y1, y2, res1, res2):
    resx = get_val(res1) 
    resy = get_val(res2) 
    p1x = get_val(x1) 
    p1y = get_val(x2) 
    p2x = get_val(y1) 
    p2y = get_val(y2) 
    
    res = fp2(resx, resy)
    add1 = fp2(p1x, p1y)
    add2 = fp2(p2x, p2y)
    res = fp2(resx, resy)
    res_calc = fp2(0,0)
    res_calc.sous(add1, add2)
    assert res == res_calc, f'fp2_sub not working properly\n should be\n {hex(res_calc.a0)} \n {hex(res_calc.a1)} but is \n {hex(res.a0)} \n {hex(res.a1)}'

def fp2_mul(x1, x2, y1, y2, res1, res2):
    
    p1x = get_val(x1) 
    p1y = get_val(x2) 
    p2x = get_val(y1) 
    p2y = get_val(y2)
    resx = get_val(res1) 
    resy = get_val(res2) 
    
    res = fp2(resx, resy)
    add1 = fp2(p1x, p1y)
    add2 = fp2(p2x, p2y)
    res = fp2(resx, resy)
    res_calc = fp2(0,0)
    res_calc.mul_test(add1, add2)
    assert res == res_calc, f'fp2_mul not working properly\n should be\n {hex(res_calc.a0)} \n {hex(res_calc.a1)} but is \n {hex(res.a0)} \n {hex(res.a1)}'

# mont
# ogomory reduction fúnctions 

def REDCL_SANITY(x1, m, r):
    val = get_val(x1)
    mod = get_val(m)
    res = get_val(r)
    res_calc = val % mod
    assert res == res_calc, f'REDCL not working properly\n should be\n {hex(res_calc)} but is \n {hex(res)} \n'

def MODMUL(x1, x2, m, r):
    mul1 = get_val(x1)
    mul2 = get_val(x2)
    mod = get_val(m)
    res = get_val(r)
    res_calc = (mul1 * mul2) % mod
    assert res == res_calc, f'MODMUL not working properly\n should be\n {hex(res_calc)} but is \n {hex(res)} \n'

def REDCLADD(x1, x2, m, r):
    mul1 = get_val(x1)
    mul2 = get_val(x2)
    mod = get_val(m)
    res = get_val(r)
    res_calc = (mul1 + mul2) % mod
    assert res == res_calc, f'REDCL + not working properly\n should be\n {hex(res_calc)} but is \n {hex(res)} \n'

def REDCLSUB(x1, x2, m, r):
    mul1 = get_val(x1)
    mul2 = get_val(x2)
    mod = get_val(m)
    res = get_val(r)
    res_calc = (mul1 - mul2) % mod
    assert res == res_calc, f'REDCL - not working properly\n should be\n {hex(res_calc)} but is \n {hex(res)} \n'

def REDCL2_SANITY(x1, y1, m, rx1, ry1):
    px = get_val(x1)
    py = get_val(y1)
    mod = get_val(m)
    resx = get_val(rx1)
    resy = get_val(ry1)
    res_calc_x = px % mod
    res_calc_y = py % mod
    assert resx == res_calc_x and resy == res_calc_y, f'REDCL2 not working properly\n should be\n {hex(res_calc_x)} \n {hex(res_calc_y)} but is \n {hex(resx)} \n {hex(resy)}'

def REDCL2ADD(x1, y1, x2, y2, m, rx1, ry1):
    px1 = get_val(x1)
    py1 = get_val(y1)
    px2 = get_val(x2)
    py2 = get_val(y2)
    mod = get_val(m)
    resx = get_val(rx1)
    resy = get_val(ry1)
    res_calc_x = (px1 + px2) % mod
    res_calc_y = (py1 + py2) % mod
    assert resx == res_calc_x and resy == res_calc_y, f'REDCL2ADD not working properly\n should be\n {hex(res_calc_x)} \n {hex(res_calc_y)} but is \n {hex(resx)} \n {hex(resy)}'

def REDCL2MUL(x1, y1, x2, y2, m, rx1, ry1):
    p1 = fp2(get_val(x1), get_val(y1))
    p2 = fp2(get_val(x2), get_val(y2))
    resx = get_val(rx1)
    resy = get_val(ry1)
    mod = get_val(m)
    res_calc = fp2(0,0)
    res_calc.mul(p1, p2)
    res_calc_x = res_calc.a0 % mod
    res_calc_y = res_calc.a1 % mod
    assert resx == res_calc_x and resy == res_calc_y, f'REDCL2MUL not working properly\n should be\n {hex(res_calc_x)} \n {hex(res_calc_y)} but is \n {hex(resx)} \n {hex(resy)}'

def fp2_mont_mul(x1, y1, x2, y2, m, rx1, ry1):
    p1 = fp2(get_val(x1), get_val(y1))
    p2 = fp2(get_val(x2), get_val(y2))
    resx = get_val(rx1)
    resy = get_val(ry1)
    mod = get_val(m)
    res_calc = fp2(0,0)
    res_calc.mul(p1, p2)
    res_calc_x = res_calc.a0 % mod
    res_calc_y = res_calc.a1 % mod
    assert resx == res_calc_x and resy == res_calc_y, f'fp2_mont_mul not working properly\n should be\n {hex(res_calc_x)} \n {hex(res_calc_y)} but is \n {hex(resx)} \n {hex(resy)}'

def xDBL(x1, y1, x2, y2, A241x, A241y, A242x, A242y, m, rx1, ry1, rx2, ry2):
    p1 = fp2(get_val(x1), get_val(y1))
    p2 = fp2(get_val(x2), get_val(y2))
    p = point(p1, p2)
    A241 = fp2(get_val(A241x), get_val(A241y))
    A242 = fp2(get_val(A242x), get_val(A242y))
    A24 = point(A241, A242)
    mod = get_val(m)
    res1 = fp2(get_val(rx1), get_val(ry1))
    res2 = fp2(get_val(rx2), get_val(ry2))
    res = point(res1, res2)
    res_calc = point(fp2(0,0), fp2(0,0))
    res_calc.double(p, A24)
    assert res.X == res_calc.X and res.Z == res_calc.Z, f'xDBL not working properly\n should be\n {hex(res_calc.X.a0)} \n {hex(res_calc.X.a1)} \n {hex(res_calc.Z.a0)} \n {hex(res_calc.Z.a1)} \n but is \n {hex(res.X.a0)} \n {hex(res.X.a1)} \n {hex(res.Z.a0)} \n {hex(res.Z.a1)}'

def xDBLe(x1, y1, x2, y2, A241x, A241y, A242x, A242y, e, m, rx1, ry1, rx2, ry2):
    p1 = fp2(get_val(x1), get_val(y1))
    p2 = fp2(get_val(x2), get_val(y2))
    p = point(p1, p2)
    A241 = fp2(get_val(A241x), get_val(A241y))
    A242 = fp2(get_val(A242x), get_val(A242y))
    A24 = point(A241, A242)
    mod = get_val(m)
    e = int(e)
    res1 = fp2(get_val(rx1), get_val(ry1))
    res2 = fp2(get_val(rx2), get_val(ry2))
    res = point(res1, res2)
    res_calc = point()
    res_calc.puissance2_k(e, p, A24)
    assert res.X == res_calc.X and res.Z == res_calc.Z, f'xDBLe not working properly\n should be\n {hex(res_calc.X.a0)} \n {hex(res_calc.X.a1)} \n {hex(res_calc.Z.a0)} \n {hex(res_calc.Z.a1)} \n but is \n {hex(res.X.a0)} \n {hex(res.X.a1)} \n {hex(res.Z.a0)} \n {hex(res.Z.a1)}'

def xTPL(x1, y1, x2, y2, A241x, A241y, A242x, A242y, m, rx1, ry1, rx2, ry2):
    p1 = fp2(get_val(x1), get_val(y1))
    p2 = fp2(get_val(x2), get_val(y2))
    p = point(p1, p2)
    A241 = fp2(get_val(A241x), get_val(A241y))
    A242 = fp2(get_val(A242x), get_val(A242y))
    A24 = point(A241, A242)
    mod = get_val(m)
    res1 = fp2(get_val(rx1), get_val(ry1))
    res2 = fp2(get_val(rx2), get_val(ry2))
    res = point(res1, res2)
    res_calc = point(fp2(0,0), fp2(0,0))
    res_calc.TPL(p, A24)
    assert res.X == res_calc.X and res.Z == res_calc.Z, f'xTPL not working properly\n should be\n {hex(res_calc.X.a0)} \n {hex(res_calc.X.a1)} \n {hex(res_calc.Z.a0)} \n {hex(res_calc.Z.a1)} \n but is \n {hex(res.X.a0)} \n {hex(res.X.a1)} \n {hex(res.Z.a0)} \n {hex(res.Z.a1)}'

def xTPLe(x1, y1, x2, y2, A241x, A241y, A242x, A242y, e, m, rx1, ry1, rx2, ry2):
    p1 = fp2(get_val(x1), get_val(y1))
    p2 = fp2(get_val(x2), get_val(y2))
    p = point(p1, p2)
    A241 = fp2(get_val(A241x), get_val(A241y))
    A242 = fp2(get_val(A242x), get_val(A242y))
    A24 = point(A241, A242)
    mod = get_val(m)
    e = int(e)
    res1 = fp2(get_val(rx1), get_val(ry1))
    res2 = fp2(get_val(rx2), get_val(ry2))
    res = point(res1, res2)
    res_calc = point()
    res_calc.TPLe(e, p, A24)
    assert res.X == res_calc.X and res.Z == res_calc.Z, f'xTPLe not working properly\n should be\n {hex(res_calc.X.a0)} \n {hex(res_calc.X.a1)} \n {hex(res_calc.Z.a0)} \n {hex(res_calc.Z.a1)} \n but is \n {hex(res.X.a0)} \n {hex(res.X.a1)} \n {hex(res.Z.a0)} \n {hex(res.Z.a1)}'

def xDBLADD(xp1, yp1, xq1, yq1, xp2, yp2, xq2, yq2, xp3, yp3, xq3, yq3, A241x, A241y, A242x, A242y, m, rxp1, ryp1, rxq1, ryq1, rxp2, ryp2, rxq2, ryq2):
    p1 = fp2(get_val(xp1), get_val(yp1))
    p2 = fp2(get_val(xq1), get_val(yq1))
    add1 = point(p1, p2)
    
    p1 = fp2(get_val(xp2), get_val(yp2))
    p2 = fp2(get_val(xq2), get_val(yq2))
    add2 = point(p1, p2)

    p1 = fp2(get_val(xp3), get_val(yp3))
    p2 = fp2(get_val(xq3), get_val(yq3))
    add3 = point(p1, p2)


    A241 = fp2(get_val(A241x), get_val(A241y))
    A242 = fp2(get_val(A242x), get_val(A242y))
    A24 = point(A241, A242)
    mod = get_val(m)

    resp1 = fp2(get_val(rxp1), get_val(ryp1))
    resp2 = fp2(get_val(rxq1), get_val(ryq1))
    res1 = point(resp1, resp2)

    resp1 = fp2(get_val(rxp2), get_val(ryp2))
    resp2 = fp2(get_val(rxq2), get_val(ryq2))
    res2 = point(resp1, resp2)
    res_calc1, res_calc2 = DBLADD(add1, add2, add3, A24,)
    assert res1.X == res_calc1.X and res1.Z == res_calc1.Z, f'xDBLADD 1 not working properly\n should be\n {hex(res_calc1.X.a0)} \n {hex(res_calc1.X.a1)} \n {hex(res_calc1.Z.a0)} \n {hex(res_calc1.Z.a1)} \n but is \n {hex(res1.X.a0)} \n {hex(res1.X.a1)} \n {hex(res1.Z.a0)} \n {hex(res1.Z.a1)}'
    assert res2.X == res_calc2.X and res2.Z == res_calc2.Z, f'xDBLADD 2 not working properly\n should be\n {hex(res_calc2.X.a0)} \n {hex(res_calc2.X.a1)} \n {hex(res_calc2.Z.a0)} \n {hex(res_calc2.Z.a1)} \n but is \n {hex(res2.X.a0)} \n {hex(res2.X.a1)} \n {hex(res2.Z.a0)} \n {hex(res2.Z.a1)}'

def jinvariant(x1, y1, x2, y2, m, rx1, ry1):
    A = fp2(get_val(x1), get_val(y1))
    C = fp2(get_val(x2), get_val(y2))
    AC = point(A, C)
    mod = get_val(m)
    res = fp2(get_val(rx1), get_val(ry1))
    res_calc = j_invariance(AC)
    assert res == res_calc , f'jinvariant not working properly\n should be\n {hex(res_calc.a0)} \n {hex(res_calc.a1)} but is \n {hex(res.a0)} \n {hex(res.a1)}'

def iso_2_curve_(x1, y1, x2, y2, m, rx1, ry1, rx2, ry2):
    p1 = fp2(get_val(x1), get_val(y1))
    p2 = fp2(get_val(x2), get_val(y2))
    p = point(p1, p2)
    mod = get_val(m)
    res1 = fp2(get_val(rx1), get_val(ry1))
    res2 = fp2(get_val(rx2), get_val(ry2))
    res = point(res1, res2)
    #res.mod()
    res_calc = iso_2_curve(p)
    res_calc.mod()
    assert res.X == res_calc.X and res.Z == res_calc.Z, f'iso_2_curve not working properly\n should be\n {hex(res_calc.X.a0)} \n {hex(res_calc.X.a1)} \n {hex(res_calc.Z.a0)} \n {hex(res_calc.Z.a1)} \n but is \n {hex(res.X.a0)} \n {hex(res.X.a1)} \n {hex(res.Z.a0)} \n {hex(res.Z.a1)}'

def mont_mult_inv(x1, y1, m, r1, r2):
    val1 = get_val(x1)
    val2 = get_val(y1)
    val = fp2(val1, val2)
    mod = get_val(m)
    res = fp2(get_val(r1), get_val(r2))
    res_calc = fp2()
    res_calc.inv(val)
    assert res == res_calc, f'mont_mult_inv not working properly\n should be\n {hex(res_calc.a0)} \n {hex(res_calc.a1)} but is \n {hex(res.a0)} \n {hex(res.a1)}\n'

def iso_2_eval_(x1, y1, x2, y2, x3, y3, x4, y4, m, rx1, ry1, rx2, ry2):
    px1 = fp2(get_val(x1), get_val(y1))
    py1 = fp2(get_val(x2), get_val(y2))
    p1 = point(px1, py1)
    px2 = fp2(get_val(x3), get_val(y3))
    py2 = fp2(get_val(x4), get_val(y4))
    p2 = point(px2, py2)
    mod = get_val(m)
    res1 = fp2(get_val(rx1), get_val(ry1))
    res2 = fp2(get_val(rx2), get_val(ry2))
    res = point(res1, res2)
    res_calc = point(fp2(0,0), fp2(0,0))
    res_calc.isogeny_2_point(p2, p1)
    assert res.X == res_calc.X and res.Z == res_calc.Z, f'iso_2_eval not working properly\n should be\n {hex(res_calc.X.a0)} \n {hex(res_calc.X.a1)} \n {hex(res_calc.Z.a0)} \n {hex(res_calc.Z.a1)} \n but is \n {hex(res.X.a0)} \n {hex(res.X.a1)} \n {hex(res.Z.a0)} \n {hex(res.Z.a1)}'

def iso_4_curve_(x1, y1, x2, y2, K1x, K1y, K2x, K2y, K3x, K3y, m, rx1, ry1, rx2, ry2):
    p1 = fp2(get_val(x1), get_val(y1))
    p2 = fp2(get_val(x2), get_val(y2))
    P4 = point(p1, p2)
    K1 = fp2(get_val(K1x), get_val(K1y))
    K2 = fp2(get_val(K2x), get_val(K2y))
    K3 = fp2(get_val(K3x), get_val(K3y))
    mod = get_val(m)
    res1 = fp2(get_val(rx1), get_val(ry1))
    res2 = fp2(get_val(rx2), get_val(ry2))
    res = point(res1, res2)
    res_calc = point()
    K1_calc = fp2()
    K2_calc = fp2()
    K3_calc = fp2()
    (res_calc, K1_calc, K2_calc, K3_calc) = iso_4_curve(P4)
    K1.mod()
    K2.mod()
    K3.mod()
    K1_calc.mod()
    K2_calc.mod()
    K3_calc.mod()
    assert res.X == res_calc.X and res.Z == res_calc.Z, f'iso_4_curve not working properly\n AC point should be\n {hex(res_calc.X.a0)} \n {hex(res_calc.X.a1)} \n {hex(res_calc.Z.a0)} \n {hex(res_calc.Z.a1)} \n but is \n {hex(res.X.a0)} \n {hex(res.X.a1)} \n {hex(res.Z.a0)} \n {hex(res.Z.a1)}'
    assert (K1 == K1_calc) , f'iso_4_curve not working properly\n K1 point should be\n {hex(K1_calc.a0)} \n {hex(K1_calc.a1)} \n but is \n {hex(K1.a0)} \n {hex(K1.a1)}'
    assert (K2 == K2_calc) , f'iso_4_curve not working properly\n K2 point should be\n {hex(K2_calc.a0)} \n {hex(K2_calc.a1)} \n but is \n {hex(K2.a0)} \n {hex(K2.a1)}'
    assert (K3 == K3_calc) , f'iso_4_curve not working properly\n K3 point should be\n {hex(K3_calc.a0)} \n {hex(K3_calc.a1)} \n but is \n {hex(K3.a0)} \n {hex(K3.a1)}'

def iso_4_eval_(x1, y1, x2, y2, K1x, K1y, K2x, K2y, K3x, K3y, m, rx1, ry1, rx2, ry2):
    p1 = fp2(get_val(x1), get_val(y1))
    p2 = fp2(get_val(x2), get_val(y2))
    P = point(p1, p2)
    K1 = fp2(get_val(K1x), get_val(K1y))
    K2 = fp2(get_val(K2x), get_val(K2y))
    K3 = fp2(get_val(K3x), get_val(K3y))
    mod = get_val(m)
    res1 = fp2(get_val(rx1), get_val(ry1))
    res2 = fp2(get_val(rx2), get_val(ry2))
    res = point(res1, res2)
    P4 = point()
    res_calc = point()
    res_calc.isogeny_4_point(P, P4, K1, K2, K3)
    assert res.X == res_calc.X and res.Z == res_calc.Z, f'iso_4_eval not working properly\nshould be\n {hex(res_calc.X.a0)} \n {hex(res_calc.X.a1)} \n {hex(res_calc.Z.a0)} \n {hex(res_calc.Z.a1)} \n but is \n {hex(res.X.a0)} \n {hex(res.X.a1)} \n {hex(res.Z.a0)} \n {hex(res.Z.a1)}'

def iso_3_curve_(x1, y1, x2, y2, K1x, K1y, K2x, K2y, m, rx1, ry1, rx2, ry2):
    p1 = fp2(get_val(x1), get_val(y1))
    p2 = fp2(get_val(x2), get_val(y2))
    P3 = point(p1, p2)
    K1 = fp2(get_val(K1x), get_val(K1y))
    K2 = fp2(get_val(K2x), get_val(K2y))
    mod = get_val(m)
    res1 = fp2(get_val(rx1), get_val(ry1))
    res2 = fp2(get_val(rx2), get_val(ry2))
    res = point(res1, res2)
    res_calc = point()
    K1_calc = fp2()
    K2_calc = fp2()
    (res_calc, K1_calc, K2_calc) = iso_3_curve(P3)
    K1.mod()
    K2.mod()
    K1_calc.mod()
    K2_calc.mod()
    assert res.X == res_calc.X and res.Z == res_calc.Z, f'iso_4_curve not working properly\n AC point should be\n {hex(res_calc.X.a0)} \n {hex(res_calc.X.a1)} \n {hex(res_calc.Z.a0)} \n {hex(res_calc.Z.a1)} \n but is \n {hex(res.X.a0)} \n {hex(res.X.a1)} \n {hex(res.Z.a0)} \n {hex(res.Z.a1)}'
    assert (K1 == K1_calc) , f'iso_3_curve not working properly\n K1 point should be\n {hex(K1_calc.a0)} \n {hex(K1_calc.a1)} \n but is \n {hex(K1.a0)} \n {hex(K1.a1)}'
    assert (K2 == K2_calc) , f'iso_3_curve not working properly\n K2 point should be\n {hex(K2_calc.a0)} \n {hex(K2_calc.a1)} \n but is \n {hex(K2.a0)} \n {hex(K2.a1)}'

def iso_3_eval_(x1, y1, x2, y2, K1x, K1y, K2x, K2y, m, rx1, ry1, rx2, ry2):
    p1 = fp2(get_val(x1), get_val(y1))
    p2 = fp2(get_val(x2), get_val(y2))
    P = point(p1, p2)
    K1 = fp2(get_val(K1x), get_val(K1y))
    K2 = fp2(get_val(K2x), get_val(K2y))
    mod = get_val(m)
    res1 = fp2(get_val(rx1), get_val(ry1))
    res2 = fp2(get_val(rx2), get_val(ry2))
    res = point(res1, res2)
    P3 = point()
    res_calc = point()
    res_calc.isogeny_3_point(P, P3, K1, K2)
    assert res.X == res_calc.X and res.Z == res_calc.Z, f'iso_3_eval not working properly\nshould be\n {hex(res_calc.X.a0)} \n {hex(res_calc.X.a1)} \n {hex(res_calc.Z.a0)} \n {hex(res_calc.Z.a1)} \n but is \n {hex(res.X.a0)} \n {hex(res.X.a1)} \n {hex(res.Z.a0)} \n {hex(res.Z.a1)}'

def ladder3p_(x1, y1, x2, y2, x3, y3, secred, A241x, A241y, A242x, A242y, m, rx1, ry1, rx2, ry2):
    x1 = fp2(get_val(x1), get_val(y1))
    x2 = fp2(get_val(x2), get_val(y2))
    x3 = fp2(get_val(x3), get_val(y3))
    secred_key = get_val(secred)
    P  = point(x1, fp2(1, 0))
    Q  = point(x2, fp2(1, 0))
    PQ = point(x3, fp2(1, 0))
    A241 = fp2(get_val(A241x), get_val(A241y))
    A242 = fp2(get_val(A242x), get_val(A242y))
    A24 = point(A241, A242)
    mod = get_val(m)
    res1 = fp2(get_val(rx1), get_val(ry1))
    res2 = fp2(get_val(rx2), get_val(ry2))
    res = point(res1, res2)
    res_calc = point()
    res_calc.ladder3pt(secred_key, P, Q, PQ, A24)
    assert res.X == res_calc.X and res.Z == res_calc.Z, f'ladder3p not working properly\nshould be\n {hex(res_calc.X.a0)} \n {hex(res_calc.X.a1)} \n {hex(res_calc.Z.a0)} \n {hex(res_calc.Z.a1)} \n but is \n {hex(res.X.a0)} \n {hex(res.X.a1)} \n {hex(res.Z.a0)} \n {hex(res.Z.a1)}'

def e_2_iso_(xp1, yp1, xq1, yq1, A241x, A241y, A242x, A242y, m, A24n1x, A24n1y, A24n2x, A24n2y, K1_1x, K1_1y, K1_2x, K1_2y, K2_1x, K2_1y, K2_2x, K2_2y, K3_1x, K3_1y, K3_2x, K3_2y):
    p1 = fp2(get_val(xp1), get_val(yp1))
    p2 = fp2(get_val(xq1), get_val(yq1))
    S = point(p1, p2)

    A241 = fp2(get_val(A241x), get_val(A241y))
    A242 = fp2(get_val(A242x), get_val(A242y))
    A24 = point(A241, A242)

    mod = get_val(m)

    A24res1 = fp2(get_val(A24n1x), get_val(A24n1y))
    A24res2 = fp2(get_val(A24n2x), get_val(A24n2y))

    A24res = point(A24res1, A24res2)

    K1x = fp2(get_val(K1_1x), get_val(K1_1y))
    K1y = fp2(get_val(K1_2x), get_val(K1_2y))

    K1 = point(K1x, K1y)

    K2x = fp2(get_val(K2_1x), get_val(K2_1y))
    K2y = fp2(get_val(K2_2x), get_val(K2_2y))

    K2 = point(K2x, K2y)

    K3x = fp2(get_val(K3_1x), get_val(K3_1y))
    K3y = fp2(get_val(K3_2x), get_val(K3_2y))

    K3 = point(K3x, K3y)

    A24_calc = point()
    K1_calc = point()
    K2_calc = point()
    K3_calc = point()

    A24_calc2 = point()
    K1_calc2 = point()
    K2_calc2 = point()
    K3_calc2 = point()

    #(A24_calc2, K1_calc2, K2_calc2, K3_calc2) = e_2_iso(S, A24)

    #(A24_calc2, K1_calc2, K2_calc2, K3_calc2) = e_2_iso_strat(S, A24)

    #
    #(A24_calc, K1_calc, K2_calc, K3_calc) = e_2_iso_my_strat(S, A24)
    (A24_calc, K1_calc, K2_calc, K3_calc) = e_2_iso_x(S, A24, point(), point(), point())

    assert A24res.X == A24_calc.X and A24res.Z == A24_calc.Z, f'e_2_iso 1 not working properly\n should be\n {hex(A24_calc.X.a0)} \n {hex(A24_calc.X.a1)} \n {hex(A24_calc.Z.a0)} \n {hex(A24_calc.Z.a1)} \n but is \n {hex(A24res.X.a0)} \n {hex(A24res.X.a1)} \n {hex(A24res.Z.a0)} \n {hex(A24res.Z.a1)}'
    assert K1.X == K1_calc.X and K1.Z == K1_calc.Z, f'e_2_iso 2 not working properly\n should be\n {hex(K1_calc.X.a0)} \n {hex(K1_calc.X.a1)} \n {hex(K1_calc.Z.a0)} \n {hex(K1_calc.Z.a1)} \n but is \n {hex(K1.X.a0)} \n {hex(K1.X.a1)} \n {hex(K1.Z.a0)} \n {hex(K1.Z.a1)}'
    assert K2.X == K2_calc.X and K2.Z == K2_calc.Z, f'e_2_iso 2 not working properly\n should be\n {hex(K2_calc.X.a0)} \n {hex(K2_calc.X.a1)} \n {hex(K2_calc.Z.a0)} \n {hex(K2_calc.Z.a1)} \n but is \n {hex(K2.X.a0)} \n {hex(K2.X.a1)} \n {hex(K2.Z.a0)} \n {hex(K2.Z.a1)}'
    assert K3.X == K3_calc.X and K3.Z == K3_calc.Z, f'e_2_iso 2 not working properly\n should be\n {hex(K3_calc.X.a0)} \n {hex(K3_calc.X.a1)} \n {hex(K3_calc.Z.a0)} \n {hex(K3_calc.Z.a1)} \n but is \n {hex(K3.X.a0)} \n {hex(K3.X.a1)} \n {hex(K3.Z.a0)} \n {hex(K3.Z.a1)}'

    #assert A24_calc2.X == A24_calc.X and A24_calc2.Z == A24_calc.Z, f'e_2_iso 1 not working properly\n should be\n {hex(A24_calc.X.a0)} \n {hex(A24_calc.X.a1)} \n {hex(A24_calc.Z.a0)} \n {hex(A24_calc.Z.a1)} \n but is \n {hex(A24_calc2.X.a0)} \n {hex(A24_calc2.X.a1)} \n {hex(A24_calc2.Z.a0)} \n {hex(A24_calc2.Z.a1)}'
    #assert K1_calc2.X == K1_calc.X and K1_calc2.Z == K1_calc.Z, f'e_2_iso 2 not working properly\n should be\n {hex(K1_calc.X.a0)} \n {hex(K1_calc.X.a1)} \n {hex(K1_calc.Z.a0)} \n {hex(K1_calc.Z.a1)} \n but is \n {hex(K1_calc2.X.a0)} \n {hex(K1_calc2.X.a1)} \n {hex(K1_calc2.Z.a0)} \n {hex(K1_calc2.Z.a1)}'
    #assert K2_calc2.X == K2_calc.X and K2_calc2.Z == K2_calc.Z, f'e_2_iso 2 not working properly\n should be\n {hex(K2_calc.X.a0)} \n {hex(K2_calc.X.a1)} \n {hex(K2_calc.Z.a0)} \n {hex(K2_calc.Z.a1)} \n but is \n {hex(K2_calc2.X.a0)} \n {hex(K2_calc2.X.a1)} \n {hex(K2_calc2.Z.a0)} \n {hex(K2_calc2.Z.a1)}'
    #assert K3_calc2.X == K3_calc.X and K3_calc2.Z == K3_calc.Z, f'e_2_iso 2 not working properly\n should be\n {hex(K3_calc.X.a0)} \n {hex(K3_calc.X.a1)} \n {hex(K3_calc.Z.a0)} \n {hex(K3_calc.Z.a1)} \n but is \n {hex(K3_calc2.X.a0)} \n {hex(K3_calc2.X.a1)} \n {hex(K3_calc2.Z.a0)} \n {hex(K3_calc2.Z.a1)}'

def e_3_iso_(xp1, yp1, xq1, yq1, A241x, A241y, A242x, A242y, m, A24n1x, A24n1y, A24n2x, A24n2y, K1_1x, K1_1y, K1_2x, K1_2y, K2_1x, K2_1y, K2_2x, K2_2y, K3_1x, K3_1y, K3_2x, K3_2y):
    p1 = fp2(get_val(xp1), get_val(yp1))
    p2 = fp2(get_val(xq1), get_val(yq1))
    S = point(p1, p2)

    A241 = fp2(get_val(A241x), get_val(A241y))
    A242 = fp2(get_val(A242x), get_val(A242y))
    A24 = point(A241, A242)

    mod = get_val(m)

    A24res1 = fp2(get_val(A24n1x), get_val(A24n1y))
    A24res2 = fp2(get_val(A24n2x), get_val(A24n2y))

    A24res = point(A24res1, A24res2)

    K1x = fp2(get_val(K1_1x), get_val(K1_1y))
    K1y = fp2(get_val(K1_2x), get_val(K1_2y))

    K1 = point(K1x, K1y)

    K2x = fp2(get_val(K2_1x), get_val(K2_1y))
    K2y = fp2(get_val(K2_2x), get_val(K2_2y))

    K2 = point(K2x, K2y)

    K3x = fp2(get_val(K3_1x), get_val(K3_1y))
    K3y = fp2(get_val(K3_2x), get_val(K3_2y))

    K3 = point(K3x, K3y)

    A24_calc = point()
    K1_calc = point()
    K2_calc = point()
    K3_calc = point()

    A24_calc2 = point()
    K1_calc2 = point()
    K2_calc2 = point()
    K3_calc2 = point()

    #(A24_calc2, K1_calc2, K2_calc2, K3_calc2) = e_2_iso(S, A24)

    #(A24_calc2, K1_calc2, K2_calc2, K3_calc2) = e_2_iso_strat(S, A24)

    (A24_calc, K1_calc, K2_calc, K3_calc) = e_3_iso_strat(S, A24)
    #(A24_calc, K1_calc, K2_calc, K3_calc) = e_2_iso_strat(S, A24)

    assert A24res.X == A24_calc.X and A24res.Z == A24_calc.Z, f'e_3_iso 1 not working properly\n should be\n {hex(A24_calc.X.a0)} \n {hex(A24_calc.X.a1)} \n {hex(A24_calc.Z.a0)} \n {hex(A24_calc.Z.a1)} \n but is \n {hex(A24res.X.a0)} \n {hex(A24res.X.a1)} \n {hex(A24res.Z.a0)} \n {hex(A24res.Z.a1)}'
    assert K1.X == K1_calc.X and K1.Z == K1_calc.Z, f'e_3_iso 2 not working properly\n should be\n {hex(K1_calc.X.a0)} \n {hex(K1_calc.X.a1)} \n {hex(K1_calc.Z.a0)} \n {hex(K1_calc.Z.a1)} \n but is \n {hex(K1.X.a0)} \n {hex(K1.X.a1)} \n {hex(K1.Z.a0)} \n {hex(K1.Z.a1)}'
    assert K2.X == K2_calc.X and K2.Z == K2_calc.Z, f'e_3_iso 2 not working properly\n should be\n {hex(K2_calc.X.a0)} \n {hex(K2_calc.X.a1)} \n {hex(K2_calc.Z.a0)} \n {hex(K2_calc.Z.a1)} \n but is \n {hex(K2.X.a0)} \n {hex(K2.X.a1)} \n {hex(K2.Z.a0)} \n {hex(K2.Z.a1)}'
    assert K3.X == K3_calc.X and K3.Z == K3_calc.Z, f'e_3_iso 2 not working properly\n should be\n {hex(K3_calc.X.a0)} \n {hex(K3_calc.X.a1)} \n {hex(K3_calc.Z.a0)} \n {hex(K3_calc.Z.a1)} \n but is \n {hex(K3.X.a0)} \n {hex(K3.X.a1)} \n {hex(K3.Z.a0)} \n {hex(K3.Z.a1)}'

    #assert A24_calc2.X == A24_calc.X and A24_calc2.Z == A24_calc.Z, f'e_2_iso 1 not working properly\n should be\n {hex(A24_calc.X.a0)} \n {hex(A24_calc.X.a1)} \n {hex(A24_calc.Z.a0)} \n {hex(A24_calc.Z.a1)} \n but is \n {hex(A24_calc2.X.a0)} \n {hex(A24_calc2.X.a1)} \n {hex(A24_calc2.Z.a0)} \n {hex(A24_calc2.Z.a1)}'
    #assert K1_calc2.X == K1_calc.X and K1_calc2.Z == K1_calc.Z, f'e_2_iso 2 not working properly\n should be\n {hex(K1_calc.X.a0)} \n {hex(K1_calc.X.a1)} \n {hex(K1_calc.Z.a0)} \n {hex(K1_calc.Z.a1)} \n but is \n {hex(K1_calc2.X.a0)} \n {hex(K1_calc2.X.a1)} \n {hex(K1_calc2.Z.a0)} \n {hex(K1_calc2.Z.a1)}'
    #assert K2_calc2.X == K2_calc.X and K2_calc2.Z == K2_calc.Z, f'e_2_iso 2 not working properly\n should be\n {hex(K2_calc.X.a0)} \n {hex(K2_calc.X.a1)} \n {hex(K2_calc.Z.a0)} \n {hex(K2_calc.Z.a1)} \n but is \n {hex(K2_calc2.X.a0)} \n {hex(K2_calc2.X.a1)} \n {hex(K2_calc2.Z.a0)} \n {hex(K2_calc2.Z.a1)}'
    #assert K3_calc2.X == K3_calc.X and K3_calc2.Z == K3_calc.Z, f'e_2_iso 2 not working properly\n should be\n {hex(K3_calc.X.a0)} \n {hex(K3_calc.X.a1)} \n {hex(K3_calc.Z.a0)} \n {hex(K3_calc.Z.a1)} \n but is \n {hex(K3_calc2.X.a0)} \n {hex(K3_calc2.X.a1)} \n {hex(K3_calc2.Z.a0)} \n {hex(K3_calc2.Z.a1)}'

def get_A_(x1, y1, x2, y2, x3, y3, rx, ry):
    xp = fp2(get_val(x1), get_val(y1))
    xq = fp2(get_val(x2), get_val(y2))
    xpq = fp2(get_val(x3), get_val(y3))
    res = fp2(get_val(rx), get_val(ry))
    res.mod()
    res_calc = fp2()
    res_calc = get_A(xp, xq, xpq)
    res_calc.mod()
    assert res == res_calc, f'get_A not working properly\n should be\n {hex(res_calc.a0)} \n {hex(res_calc.a1)} \n but is \n {hex(res.a0)} \n {hex(res.a1)} \n'

try:
    with open(path) as file:
        line = file.readlines()
        i = 1
        if line[0] == '434\n':
            bits = 64 * 8#6
            e1 = 0xD8
            e2 = 0x89
            change_params(e1, e2)
        elif line[0] == '503\n':
            bits = 64 * 9#7
            e1 = 0xFA
            e2 = 0x9F
            change_params(e1, e2)
        elif line[0] == '610\n':
            bits = 64 * 11#9
            e1 = 0x131
            e2 = 0xC0
            change_params(e1, e2)
        elif line[0] == '751\n':
            bits = 64 * 13#12
            e1 = 0x174
            e2 = 0xEF
            print(hex(p))
            change_params(e1, e2)
            print(hex(p))
        else:
            print(f'Wrong mode {line[0]}')
            exit()
        print(f'Mode: {line[0]}')
        print(hex(e1))
        print(hex(e2))
        print(hex(p))
        while i < len(line):
            #print(i)
            if i % 1000 == 0:
                print(i)
            teststr = line[i]
            #match line[i]:
            #    case '+\n':
            if teststr == '+\n':
                fp_add(line[i+1], line[i+2], line[i+3])
                i = i + 3
            #case '-\n':
            elif teststr == '-\n':
                fp_sub(line[i+1], line[i+2], line[i+3])
                i = i + 3
            #case '*\n':
            elif teststr == '*\n':
                fp_mul(line[i+1], line[i+2], line[i+3])
                i = i + 3
            #case "/\n":
            elif teststr == '/\n':
                fp_div(line[i+1], line[i+2], line[i+3])
                i = i + 3
            #case "%\n":
            elif teststr == '%\n':
                fp_mod(line[i+1], line[i+2], line[i+3])
                i = i + 3
            #case "%2k\n":
            elif teststr == '%2k\n':
                fp_mod_2k(line[i+1], line[i+2], line[i+3])
                i = i + 3
            elif teststr == '+%\n':
                fp_addm(line[i+1], line[i+2], line[i+3], line[i+4])
                i = i + 4
            elif teststr == '-%\n':
                fp_subm(line[i+1], line[i+2], line[i+3], line[i+4])
                i = i + 4
            #case '++\n':
            elif teststr == '++\n':
                f2p_add(line[i+1], line[i+2], line[i+3])
                i = i + 3
            #case '--\n':
            elif teststr == '--\n':
                f2p_sub(line[i+1], line[i+2], line[i+3])
                i = i + 3
            #case '**\n':
            elif teststr == '**\n':
                f2p_mul(line[i+1], line[i+2], line[i+3])
                i = i + 3  
            #case "//\n":
            elif teststr == '//\n':
                f2p_div(line[i+1], line[i+2], line[i+3])
                i = i + 3
            #case "%%\n":
            elif teststr == '%%\n':
                f2p_mod(line[i+1], line[i+2], line[i+3])
                i = i + 3
            #case "//2k\n":
            elif teststr == '//2k\n':
                f2p_div_2k(line[i+1], line[i+2], line[i+3])
                i = i + 3
            #case "%%2k\n":
            elif teststr == '%%2k\n':
                f2p_mod_2k(line[i+1], line[i+2], line[i+3])
                i = i + 3
            if teststr == 'fp_get_len\n':
                fp_get_len(line[i+1], line[i+2])
                i = i + 2
            elif teststr == '+fp2\n':
                fp2_add(line[i+1], line[i+2], line[i+3], line[i+4], line[i+5], line[i+6])
                i = i + 6
            elif teststr == '-fp2\n':
                fp2_sub(line[i+1], line[i+2], line[i+3], line[i+4], line[i+5], line[i+6])
                i = i + 6
            elif teststr == '*fp2\n':
                fp2_mul(line[i+1], line[i+2], line[i+3], line[i+4], line[i+5], line[i+6])
                i = i + 6
            elif teststr == 'REDECL-Sanity\n':
                REDCL_SANITY(line[i+1], line[i+2], line[i+3])
                i = i + 3
            elif teststr == 'MODMUL\n':
                MODMUL(line[i+1], line[i+2], line[i+3], line[i+4])
                i = i + 4
            elif teststr == 'REDECL+\n':
                REDCLADD(line[i+1], line[i+2], line[i+3], line[i+4])
                i = i + 4
            elif teststr == 'REDECL-\n':
                REDCLSUB(line[i+1], line[i+2], line[i+3], line[i+4])
                i = i + 4
            elif teststr == 'REDECL2-Sanity\n':
                REDCL2_SANITY(line[i+1], line[i+2], line[i+3], line[i+4], line[i + 5])
                i = i + 5
            elif teststr == 'REDECL2+\n':
                REDCL2ADD(line[i+1], line[i+2], line[i+3], line[i+4], line[i + 5], line[i+6], line[i + 7])
                i = i + 7
            elif teststr == 'REDECL2*\n':
                REDCL2MUL(line[i+1], line[i+2], line[i+3], line[i+4], line[i + 5], line[i+6], line[i + 7])
                i = i + 7
            elif teststr == 'fp2_mul_mont\n':
                fp2_mont_mul(line[i+1], line[i+2], line[i+3], line[i+4], line[i + 5], line[i+6], line[i + 7])
                i = i + 7
            elif teststr == 'xDBL\n':
                xDBL(line[i+1], line[i+2], line[i+3], line[i+4], line[i + 5], line[i+6], line[i + 7], line[i+8], line[i+9], line[i+10], line[i+11], line[i + 12], line[i+13])
                i = i + 13
            elif teststr == 'xDBLe\n':
                xDBLe(line[i+1], line[i+2], line[i+3], line[i+4], line[i + 5], line[i+6], line[i + 7], line[i+8], line[i+9], line[i+10], line[i+11], line[i + 12], line[i+13], line[i+14])
                i = i + 14
            elif teststr == 'xDBLADD\n':
                xDBLADD(line[i+1], line[i+2], line[i+3], line[i+4], line[i + 5], line[i+6], line[i + 7], line[i+8], line[i+9], line[i+10], line[i+11], line[i + 12], line[i+13], line[i+14], line[i + 15], line[i + 16], line[i + 17], line[i + 18], line[i+19], line[i+20], line[i+21], line[i + 22], line[i+23], line[i+24], line[i+25])
                i = i + 25
            elif teststr == 'xTPL\n':
                xTPL(line[i+1], line[i+2], line[i+3], line[i+4], line[i + 5], line[i+6], line[i + 7], line[i+8], line[i+9], line[i+10], line[i+11], line[i + 12], line[i+13])
                i = i + 13
            elif teststr == 'xTPLe\n':
                xTPLe(line[i+1], line[i+2], line[i+3], line[i+4], line[i + 5], line[i+6], line[i + 7], line[i+8], line[i+9], line[i+10], line[i+11], line[i + 12], line[i+13], line[i+14])
                i = i + 14
            elif teststr == 'jinvariant\n':
                jinvariant(line[i+1], line[i+2], line[i+3], line[i+4], line[i + 5], line[i+6], line[i + 7])
                i = i + 7
            elif teststr == 'iso_2_curve\n':
                iso_2_curve_(line[i+1], line[i+2], line[i+3], line[i+4], line[i + 5], line[i+6], line[i + 7], line[i+8], line[i+9])
                i = i + 9
            elif teststr == 'Mont_mult_inv\n':
                mont_mult_inv(line[i+1], line[i+2], line[i+3], line[i+4], line[i+5])
                i = i + 5
            elif teststr == '2_iso_eval\n':
                iso_2_eval_(line[i+1], line[i+2], line[i+3], line[i+4], line[i + 5], line[i+6], line[i + 7], line[i+8], line[i+9], line[i+10], line[i+11], line[i + 12], line[i+13])
                i = i + 13
            elif teststr == 'iso_4_curve\n':
                iso_4_curve_(line[i+1], line[i+2], line[i+3], line[i+4], line[i + 5], line[i+6], line[i + 7], line[i+8], line[i+9], line[i+10], line[i+11], line[i + 12], line[i+13], line[i+14], line[i + 15])
                i = i + 15
            elif teststr == 'iso_4_eval\n':
                iso_4_eval_(line[i+1], line[i+2], line[i+3], line[i+4], line[i + 5], line[i+6], line[i + 7], line[i+8], line[i+9], line[i+10], line[i+11], line[i + 12], line[i+13], line[i+14], line[i + 15])
                i = i + 15
            elif teststr == 'iso_3_curve\n':
                iso_3_curve_(line[i+1], line[i+2], line[i+3], line[i+4], line[i + 5], line[i+6], line[i + 7], line[i+8], line[i+9], line[i+10], line[i+11], line[i + 12], line[i+13])
                i = i + 13
            elif teststr == 'iso_3_eval\n':
                iso_3_eval_(line[i+1], line[i+2], line[i+3], line[i+4], line[i + 5], line[i+6], line[i + 7], line[i+8], line[i+9], line[i+10], line[i+11], line[i + 12], line[i+13])
                i = i + 13
            elif teststr == 'ladder3p\n':
                ladder3p_(line[i+1], line[i+2], line[i+3], line[i+4], line[i + 5], line[i+6], line[i + 7], line[i+8], line[i+9], line[i+10], line[i+11], line[i + 12], line[i+13], line[i+14], line[i + 15], line[i + 16])
                i = i + 16
            elif teststr == 'e_2_iso\n':
                e_2_iso_(line[i+1], line[i+2], line[i+3], line[i+4], line[i + 5], line[i+6], line[i + 7], line[i+8], line[i+9], line[i+10], line[i+11], line[i + 12], line[i+13], line[i+14], line[i + 15], line[i + 16], line[i + 17], line[i + 18], line[i+19], line[i+20], line[i+21], line[i + 22], line[i+23], line[i+24], line[i+25])
                i = i + 25
            elif teststr == 'e_3_iso\n':
                e_3_iso_(line[i+1], line[i+2], line[i+3], line[i+4], line[i + 5], line[i+6], line[i + 7], line[i+8], line[i+9], line[i+10], line[i+11], line[i + 12], line[i+13], line[i+14], line[i + 15], line[i + 16], line[i + 17], line[i + 18], line[i+19], line[i+20], line[i+21], line[i + 22], line[i+23], line[i+24], line[i+25])
                i = i + 25
            elif teststr == 'get_A\n':
                get_A_(line[i+1], line[i+2], line[i+3], line[i+4], line[i + 5], line[i+6], line[i + 7], line[i+8])
                i = i + 9
            #case other:
            else:
                i = i + 1



except FileNotFoundError:
    print(f'File: {path} could not be found.')
print("end of test")