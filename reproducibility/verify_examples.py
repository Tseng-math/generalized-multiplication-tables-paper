#!/usr/bin/env python3
"""Exact rational-interval checks for the manuscript's worked examples.

The logarithm and exponential bounds follow Appendix C. Exact zeros are
analytic identities, not numerical assertions. Requires only Python 3.
Run without optimization: python3 reproducibility/verify_examples.py
"""
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product
from pathlib import Path
import json

GRID = 10**40


def down(x):
    return Q((x * GRID).__floor__(), GRID)


def up(x):
    return Q((x * GRID).__ceil__(), GRID)


class Interval:
    def __init__(self, lo, hi=None):
        self.lo = down(Q(lo))
        self.hi = up(Q(lo if hi is None else hi))
        assert self.lo <= self.hi

    def __add__(self, other):
        other = iv(other)
        return Interval(self.lo + other.lo, self.hi + other.hi)

    __radd__ = __add__

    def __neg__(self):
        return Interval(-self.hi, -self.lo)

    def __sub__(self, other):
        return self + -iv(other)

    def __rsub__(self, other):
        return iv(other) + -self

    def __mul__(self, other):
        other = iv(other)
        candidates = [a*b for a in (self.lo, self.hi)
                      for b in (other.lo, other.hi)]
        return Interval(min(candidates), max(candidates))

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = iv(other)
        assert other.lo > 0 or other.hi < 0
        return self * Interval(1/other.hi, 1/other.lo)

    def __rtruediv__(self, other):
        return iv(other) / self

    def bounds(self, digits=8):
        grid = 10**digits
        return [str(Q((self.lo*grid).__floor__(), grid)),
                str(Q((self.hi*grid).__ceil__(), grid))]


def iv(x):
    return x if isinstance(x, Interval) else Interval(x)


@lru_cache(None)
def log_reduced(x):
    """Exact series on [1,2], with a geometric upper remainder."""
    assert 1 <= x <= 2
    z = (x-1)/(x+1)
    n = 48
    value = 2*sum((z**(2*j+1)/Q(2*j+1) for j in range(n)), Q(0))
    remainder = 2*z**(2*n+1)/((2*n+1)*(1-z*z))
    return Interval(value, value+remainder)


@lru_cache(None)
def log_point(x):
    assert x > 0
    if x < 1:
        return -log_point(1/x)
    power2 = 0
    while x > 2:
        x /= 2
        power2 += 1
    return log_reduced(x) + power2*log_reduced(Q(2))


def log(x):
    x = iv(x)
    return Interval(log_point(x.lo).lo, log_point(x.hi).hi)


@lru_cache(None)
def exp_point(x):
    if x < 0:
        value = exp_point(-x)
        return 1/value
    squares = 0
    while x > Q(1,4):
        x /= 2
        squares += 1
    n = 40
    term = total = Q(1)
    for j in range(1, n+1):
        term *= x/j
        total += term
    remainder = term*x/(n+1)/(1-x/(n+2))
    result = Interval(total, total+remainder)
    for _ in range(squares):
        result = result*result
    return result


def exp(x):
    x = iv(x)
    return Interval(exp_point(x.lo).lo, exp_point(x.hi).hi)


def power(x, y):
    return exp(iv(y)*log(x))


def alpha(h):
    return log((h-1)/log(h))/log(h)


def E(r, z):
    if r == 0:
        return iv(0)
    return (r+z)*log(r+z)-z*log(z)


def secant(r, q, z):
    return 1-log((E(r,z)-E(q,z))/(r-q))/log(r+z)


def secant_derivative(r, q, z):
    x, y = r+z, q+z
    lx, ly = log(x), log(y)
    D = (x*lx-y*ly)/(r-q)
    return log(D)/(x*lx*lx)-(lx-ly)/((r-q)*D*lx)


def positive(value, bound=0):
    assert value.lo > Q(bound), value.bounds()


def negative(value, bound=0):
    assert value.hi < Q(bound), value.bounds()


def minimum(values):
    return Interval(min(v.lo for v in values), min(v.hi for v in values))


def main():
    output = {'method': 'independent exact rational, range-reduced log/exp series',
              'rounding_grid': str(GRID), 'log_terms': 48, 'exp_degree': 40}
    a2 = alpha(2)
    output['alpha2'] = a2.bounds(15)

    # The terminal critical family. All KKT signs and all remaining edge
    # signs reduce analytically to the following finite scalar inequalities.
    s = secant(6, 1, iv(1))
    e = {h: power(h,s)*log(h)-(h-1) for h in range(2,8)}
    positive(e[2])
    for h in range(3,8):
        negative(e[h])
    margins = {}
    for h in range(3,8):
        for qweight in range(2,h):
            if (h,qweight) == (7,2):
                continue  # Exact identity from the definition of s.
            d = power(h,s-1)*(h*log(h)-qweight*log(qweight))-(h-qweight)
            positive(d)
            margins[f'{h},{qweight}'] = d
    a, b = -e[7]/e[2], -sum((e[h] for h in range(3,7)),iv(0))/e[2]
    positive(a-54); negative(a-55)
    positive(b-66); negative(b-67)
    # tau is at most its explicit logarithmic component.
    tau_upper = 1/(8*7**3*log(7)*log(7))
    positive(s-a2-tau_upper); positive(alpha(3)-s-tau_upper)
    low_thresholds = [secant(H-1,1,iv(1))-a2 for H in range(3,7)]
    for v in low_thresholds:
        negative(v)
    output['critical_six'] = {
        'sigma':s.bounds(12), 'a':a.bounds(8), 'b':b.bounds(8),
        'e2':e[2].bounds(8),
        'min_strict_secant_margin':minimum(list(margins.values())).bounds(8),
        'lower_dimensional_sigma_minus_alpha': [v.bounds(8) for v in low_thresholds],
        'exact_zero':'Delta_7(2)=0 by defining identity'}

    # The strict, widely separated six-dimensional family.
    def first_split(H):
        return power(H,a2-1)*(H*log(H)-2*log(2))-(H-2)
    positive(first_split(6),Q(25,1000)); negative(first_split(6),Q(26,1000))
    positive(first_split(7),Q(-110,1000)); negative(first_split(7),Q(-109,1000))
    positive(a2*(1-a2)*log(2)-(2*a2-1))
    def split_equation(product_s):
        z = power(2,a2/product_s)
        return power(5+z,product_s-1)*E(5,z)-5
    negative(split_equation(iv('0.54070')))
    positive(split_equation(iv('0.54072')))
    sstar = Interval('0.54070','0.54072')
    zstar = power(2,a2/sstar)
    cstar = 1-log(2)*zstar*power(5+zstar,sstar-1)
    positive(cstar)
    output['split_six'] = {'f6':first_split(6).bounds(8),
                           'f7':first_split(7).bounds(8),
                           'sstar':['54070/100000','54072/100000'],
                           'incoming_surplus':cstar.bounds(8)}

    # Positive ray. KKT derivatives are rederived by differentiating the
    # common owner in initial coordinates. Templates use integer deadline
    # indices: cut index d retains name c in group g exactly when d>=g.
    s, theta = iv(Q(107,200)), iv(Q(267,500))
    z = power(2,theta/s)
    kappa = {r:power(r+z,s-1) for r in range(1,6)}
    f = {r:kappa[r]*E(r,z)-r for r in range(1,6)}
    c = {r:1-log(2)*z*kappa[r] for r in range(1,6)}
    a = -sum((f[r] for r in range(1,5)),iv(0))/f[5]
    b = (a*c[5]+sum((c[r] for r in range(1,5)),iv(0)))/(power(2,theta)*log(2)-1)
    positive(a);positive(b)
    lengths = [None,a,iv(1),iv(1),iv(1),iv(1)]
    ray_multipliers = []
    for j in range(1,5):
        q=5-j
        lam=sum((lengths[i]*(kappa[6-i]*(E(6-i,z)-E(q,z))-(6-i-q))
                 for i in range(1,j+1)),iv(0))
        positive(lam)
        ray_multipliers.append(lam.bounds(5))
    native_cost = {(g,q):q-kappa[6-g]*E(q,z)
                   for g in range(1,6) for q in range(7-g)}
    ray_edges = []
    for i in range(1,6):
        r=6-i
        for q in range(r):
            earlier=list(range(1,i)); retained=list(range(6-q,6))
            outside=earlier+retained
            options=[range(cname+1) for cname in earlier]+[range(i,cname+1) for cname in retained]
            costs=[]
            side=1 if (i,q)==(1,0) else 0
            count=0
            for choices in product(*options):
                cuts=dict(zip(outside,choices))
                for cname in range(i,6-q):
                    cuts[cname]=i-1+side
                cost=sum((lengths[g]*native_cost[g,sum(cuts[cname]>=g for cname in range(1,6))]
                          for g in range(1,6)),iv(0))
                costs.append(cost);count+=1
            min_cost=minimum(costs)
            positive(min_cost,Q(1,2) if side else Q(1,10))
            ray_edges.append({'i':i,'q':q,'side':'right' if side else 'left',
                              'templates':count,'min_C_over_t':min_cost.bounds(6)})
    output['positive_ray']={'a':a.bounds(8),'b':b.bounds(8),
                            'multipliers':ray_multipliers,'edges':ray_edges}

    # Same owner, two critical comparisons. Enclose a root, and prove its
    # uniqueness on the enclosing interval by a derivative bound.
    z=Interval('1.98455','1.98457')
    difference=lambda zz:secant(14,6,zz)-secant(6,0,zz)
    negative(difference(iv(z.lo)));positive(difference(iv(z.hi)))
    derivative=secant_derivative(14,6,z)-secant_derivative(6,0,z)
    positive(derivative)
    s=secant(6,0,z);theta=s*log(z)/log(2)
    positive(theta-a2);positive(s-theta)
    kap={r:power(r+z,s-1) for r in range(1,15)}
    g={r:r-kap[r]*E(r,z) for r in range(1,15)}
    negative(g[1]);positive(g[14])
    deltas={r:kap[r]*(E(r,z)-E(6,z))-(r-6) for r in range(7,14)}
    for d in deltas.values():positive(d)
    for q in range(1,6):positive(q-kap[6]*E(q,z))
    for q in range(7,14):positive(kap[14]*(E(14,z)-E(q,z))-(14-q))
    W_T=-g[14]/g[1]
    W_C=-sum((g[r] for r in range(2,14) if r!=6),iv(0))/g[1]
    positive(W_T);positive(W_C)
    ingress={r:1-z*log(2)*kap[r] for r in range(1,15)}
    for v in ingress.values():positive(v)
    denominator=power(2,theta)*log(2)-1;positive(denominator)
    output['same_owner']={'z':z.bounds(5),'root_left':difference(iv(z.lo)).bounds(12),
                          'root_right':difference(iv(z.hi)).bounds(12),
                          'derivative':derivative.bounds(8),'sigma':s.bounds(8),
                          'theta':theta.bounds(8),'g1':g[1].bounds(8),'g14':g[14].bounds(8),
                          'W_T':W_T.bounds(4),'W_C':W_C.bounds(4),
                          'R_U':(ingress[6]/denominator).bounds(4),
                          'min_ingress':minimum(list(ingress.values())).bounds(6),
                          'middle_margins':{str(r):d.bounds(6) for r,d in deltas.items()},
                          'exact_zeros':['g6=0','f14(6)-f14(14)=0']}

    # Different owners. Each finite check is performed for the precise
    # complete chain and the fixed geometric coefficient vector 1,2,... .
    node_threshold=lambda x:1-log(x*log(x)-(x-1)*log(x-1))/log(x)
    negative(node_threshold(17)-a2);positive(node_threshold(18)-a2)
    def drift(ss, child):
        zz=exp(child/ss)
        return power(1+zz,ss-1)*E(1,zz)-1
    v=Interval('0.531542356069','0.531542356071')
    child=a2*log(17)
    negative(drift(iv(v.lo),child));positive(drift(iv(v.hi),child))
    G2=v*log(1+exp(child/v))
    u=Interval('0.5343537358','0.5343537362')
    negative(drift(iv(u.lo),G2));positive(drift(iv(u.hi),G2))
    positive(u-v);positive(v-a2)
    terminal_lam={}
    for j in range(3,18):
        qweight=19-j
        lam=sum((2**(i-3)*(power(20-i,a2-1)*((20-i)*log(20-i)-qweight*log(qweight))-(j-i+1))
                 for i in range(3,j+1)),iv(0))
        positive(lam,Q(7,10000))
        terminal_lam[str(j)]=lam.bounds(6)
    edge_bounds={};edge_count=0
    for h in range(3,18):
        absolute=[]
        for q in range(1,h-1):
            d=(h-1-q)-power(h,a2-1)*(h*log(h)-(q+1)*log(q+1))
            if d.lo>0:
                absolute.append(d)
            else:
                negative(d);absolute.append(-d)
            positive(absolute[-1],Q(7,10000));edge_count+=1
        edge_bounds[str(h)]=minimum(absolute).bounds(8)
    assert edge_count==120
    output['different_owners']={'u':u.bounds(10),'v':v.bounds(12),
                                'threshold17_minus_alpha':(node_threshold(17)-a2).bounds(8),
                                'threshold18_minus_alpha':(node_threshold(18)-a2).bounds(8),
                                'terminal_multipliers':terminal_lam,
                                'minimum_absolute_edge_drift_by_h':edge_bounds,
                                'number_terminal_edges':edge_count}
    output['result']='PASS'
    destination=Path(__file__).resolve().parent/'examples_results.json'
    destination.write_text(json.dumps(output,indent=2)+'\n')
    print('PASS: independent rational signs, root isolations, KKT and all 15 ray / 120 terminal edges')
    print(destination)


if __name__=='__main__':
    main()
