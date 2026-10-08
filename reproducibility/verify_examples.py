#!/usr/bin/env python3
"""Reproduce the arithmetic in Section 3.7, ``Two worked profiles``.

Run with Python 3.8+ from any directory; only the standard library is used.
All numerical enclosures use Fraction endpoints, outward rational rounding,
and explicit series remainders. This script writes no files.

``profile:example-two``: certify the scalar saddle, ordered KKT sign,
endpoint coefficients, edge limit, and ordinary-suffix selection.
``profile:example-split``: certify the separation threshold and limiting
two-node saddle, its internal KKT signs, and positive terminal drift.
Also enumerate every lawful endpoint-context pair at T/lambda=1000,
certifying reference nonnegativity and the exceptional edge's identities.
The latter uses normalized lengths (1000,1,1,1,1,100000); homogeneity
extends its saddle and reference signs to positive rescalings. It does
not claim lambda=1 meets the manuscript's arithmetic lift threshold.

These checks concern the worked examples' deterministic geometry. They do
not compute chi, certify its strict annealing inequality, reprove the
asymptotic theorems, or test the full manuscript. The finite face-root
certificate uses continuity and the manuscript's strict-concavity theorem
to identify the stationary point; it does not reprove that theorem.
Optimized Python (-O) is deliberately rejected, so checks cannot vanish.
"""

from fractions import Fraction as F
from itertools import accumulate, product
import sys


class VerificationError(Exception):
    """A certificate or its stated domain failed."""


def require(condition, message):
    if not condition:
        raise VerificationError(message)


# A fixed outward grid keeps exact rational denominators manageable.
SCALE = 10**24


def rounded(lower, upper):
    require(lower <= upper, "reversed interval")
    lo = lower.numerator * SCALE // lower.denominator
    hi = -((-upper.numerator * SCALE) // upper.denominator)
    return F(lo, SCALE), F(hi, SCALE)


def interval(value):
    return value if isinstance(value, tuple) else (F(value), F(value))


def add(x, y):
    x, y = interval(x), interval(y)
    return rounded(x[0] + y[0], x[1] + y[1])


def sub(x, y):
    x, y = interval(x), interval(y)
    return rounded(x[0] - y[1], x[1] - y[0])


def mul(x, y):
    x, y = interval(x), interval(y)
    values = [a * b for a in x for b in y]
    return rounded(min(values), max(values))


def div(x, y):
    y = interval(y)
    require(y[0] > 0, "division needs a positive denominator")
    return mul(x, rounded(1 / y[1], 1 / y[0]))


def total(values):
    result = interval(0)
    for value in values:
        result = add(result, value)
    return result


def log_positive(x):
    x = interval(x)
    require(x[0] >= 1, "positive log series needs x >= 1")
    # log(x)=2 sum z^(2j+1)/(2j+1), z=(x-1)/(x+1).
    z = div(sub(x, 1), add(x, 1))
    z2 = mul(z, z)
    require(z2[1] < 1, "logarithm tail ratio must be below one")
    power, answer = z, interval(0)
    terms = 120
    for j in range(terms):
        answer = add(answer, div(power, 2*j + 1))
        power = mul(power, z2)
    lower, upper = mul(2, answer)
    tail = 2*power[1] / ((2*terms + 1)*(1 - z2[1]))
    return rounded(lower, upper + tail)


def log_value(x):
    x = interval(x)
    require(x[0] > 0, "logarithm needs a positive argument")
    if x[0] >= 1:
        return log_positive(x)
    require(x[1] <= 1, "logarithm interval straddles one")
    return sub(0, log_positive(div(1, x)))


def exp_value(x):
    x = interval(x)
    if x == interval(0):
        return interval(1)
    if x[1] <= 0:
        return div(1, exp_value(sub(0, x)))
    require(0 <= x[0] <= x[1] < 2, "exponential series domain exceeded")
    term, answer = interval(1), interval(1)
    terms = 40
    for j in range(1, terms + 1):
        term = div(mul(term, x), j)
        answer = add(answer, term)
    first_tail = term[1]*x[1] / (terms + 1)
    tail = first_tail / (1 - x[1]/(terms + 2))
    return rounded(answer[0], answer[1] + tail)


def within(bounds, lower, upper, name):
    require(F(lower) < bounds[0] <= bounds[1] < F(upper), name)


LOG2 = log_value(2)
LOG3 = log_value(3)
LOG7 = log_value(7)
ALPHA2 = div(sub(0, log_value(LOG2)), LOG2)
ALPHA3 = div(log_value(div(2, LOG3)), LOG3)
ALPHA7 = div(log_value(div(6, LOG7)), LOG7)


def power(base_log, exponent):
    return exp_value(mul(base_log, exponent))


def comparable_example():
    """The k=2 equality-face saddle and its exact two contexts."""
    def stationarity(s):
        return sub(add(mul(power(LOG3, s), LOG3),
                       mul(power(LOG2, s), LOG2)), 3)

    sigma = F(54133, 100000), F(54134, 100000)
    require(stationarity(sigma[0])[1] < 0 < stationarity(sigma[1])[0],
            "k=2 scalar saddle bracket")
    within(ALPHA2, F(52876, 100000), F(52877, 100000), "alpha_2 bracket")
    require(ALPHA3[0] > F(54532, 100000), "alpha_3 lower bound")
    require(ALPHA2[1] < sigma[0] < sigma[1] < ALPHA3[0],
            "k=2 saddle lies between the one-node thresholds")
    z, w = power(LOG3, sigma), power(LOG2, sigma)
    m = sub(LOG3, mul(F(2, 3), LOG2))
    nu = sub(mul(z, m), 1)
    c1 = sub(1, mul(F(2, 3), mul(z, LOG2)))
    c2 = sub(mul(w, LOG2), 1)
    # nu=c1-c2 at the exact stationary root.
    require(nu[0] > F(1536, 10000), "k=2 proper prefix KKT sign")
    require(sub(mul(power(LOG3, ALPHA2), m), 1)[0] > 0,
            "k=2 KKT sign already holds at alpha_2")
    within(c1, F(1624, 10000), F(1625, 10000), "k=2 c1 bracket")
    within(c2, F(87, 10000), F(88, 10000), "k=2 c2 bracket")
    require(c1[0] > c2[1] > 0, "k=2 endpoint reserve ordering")
    rho = div(mul(m, m), add(mul(m, m), mul(F(2, 9), mul(LOG2, LOG2))))
    within(rho, F(7914, 10000), F(7915, 10000), "k=2 rho bracket")
    rate = sub(add(2, mul(3, sigma)), add(z, w))
    within(rate, F(3561, 10000), F(3563, 10000), "k=2 entropy rate")
    edge_limit = div(mul(2, nu), z)
    within(edge_limit, F(1695, 10000), F(1697, 10000), "k=2 edge limit")
    # For H=3, tau <= 1/(216 log(3)^2). This upper bound alone
    # separates sigma from both bands, without evaluating their minimum.
    tau_bound = div(1, mul(216, mul(LOG3, LOG3)))
    require(tau_bound[1] < F(384, 100000), "k=2 suffix band width")
    require(sigma[0] > ALPHA2[1] + tau_bound[1] and
            sigma[1] < ALPHA3[0] - tau_bound[1], "k=2 suffix selection")
    # Direct reference costs for cuts (0,1),(1,1),(0,2),(1,2).
    # The right, late full-deadline cost is exactly zero by stationarity.
    require(c1[0] > 0 and c2[0] > 0 and nu[0] > 0,
            "k=2 early/late endpoint context costs")


def split_example():
    """Threshold and limiting two-node saddle for the k=6 example."""
    E5_2 = sub(mul(7, LOG7), mul(2, LOG2))
    theta = sub(1, div(log_value(div(E5_2, 5)), LOG7))
    numerator = sub(6, mul(power(LOG7, theta), LOG7))
    denominator = sub(mul(power(LOG2, theta), LOG2), 1)
    require(ALPHA2[1] < theta[0] < theta[1] < ALPHA7[0],
            "split threshold lies between alpha_2 and alpha_7")
    require(numerator[0] > 0 and denominator[0] > 0,
            "split threshold quotient signs")
    within(div(numerator, denominator), 54, 55, "54 < c_* < 55")

    def at_z(z):
        logz, logS = log_value(z), log_value(add(5, z))
        entropy = sub(mul(add(5, z), logS), mul(z, logz))
        sigma = sub(1, div(log_value(div(entropy, 5)), logS))
        tau = div(mul(sigma, logz), LOG2)
        kappa = div(5, entropy)
        residual = sub(mul(LOG2, add(mul(100, power(LOG2, tau)),
                                    mul(z, kappa))), 101)
        return sigma, tau, residual

    z = F(198635663, 100000000), F(198635665, 100000000)
    require(at_z(z[0])[2][1] < 0 < at_z(z[1])[2][0],
            "split limiting root bracket")
    sigma, tau, _ = at_z(z)
    within(sigma, F(54038, 100000), F(54040, 100000), "split sigma_0 bracket")
    within(tau, F(53504, 100000), F(53506, 100000), "split tau_0 bracket")
    require(ALPHA2[1] < tau[0] < tau[1] < sigma[0] < sigma[1] < 1,
            "split limiting products and terminal drift")
    logz, logS = log_value(z), log_value(add(5, z))
    kappa = exp_value(mul(sub(sigma, 1), logS))
    for q in range(1, 5):
        entropy = sub(mul(add(q, z), log_value(add(q, z))), mul(z, logz))
        require(sub(q, mul(kappa, entropy))[0] > 0,
                "split limiting internal KKT sign at q=" + str(q))
    require(sub(mul(power(LOG2, tau), LOG2), 1)[0] > 0,
            "split positive terminal ordinary drift")
    p = div(z, add(5, z))
    m = sub(logS, mul(p, logz))
    private_variance = mul(mul(logz, logz), mul(p, sub(1, p)))
    require(m[0] > 0 and private_variance[0] > 0,
            "split common/private variances are both positive")


def finite_face_data(z, lengths):
    """Certify the unique top-face product at a specified z interval."""
    z = interval(z)
    logz = log_value(z)
    ranks = range(5, 0, -1)
    logs = [log_value(add(r, z)) for r in ranks]
    entropies = [sub(mul(add(r, z), lr), mul(z, logz))
                 for r, lr in zip(ranks, logs)]

    def coefficients(sigma):
        return [exp_value(mul(sub(sigma, 1), lr)) for lr in logs]

    def first(sigma):
        return total(mul(L, sub(r, mul(kap, ent)))
                     for L, r, kap, ent in zip(lengths[:5], ranks,
                                             coefficients(sigma), entropies))

    lower, upper = F(1, 2), F(7, 10)
    require(first(lower)[0] > 0 > first(upper)[1],
            "finite top-face root endpoints")
    # Uniform bracket for every z in the supplied interval. Retain a
    # positive lower and negative upper sign; an undecided midpoint stops
    # refinement and leaves the already certified enclosing bracket.
    for _ in range(44):
        midpoint = (lower + upper)/2
        value = first(midpoint)
        if value[0] > 0:
            lower = midpoint
        elif value[1] < 0:
            upper = midpoint
        else:
            break
    sigma = lower, upper
    tau = div(mul(sigma, logz), LOG2)
    kappas = coefficients(sigma)
    incoming = total(mul(L, kap) for L, kap in zip(lengths[:5], kappas))
    residual = sub(sum(lengths), mul(LOG2, add(mul(lengths[5], power(LOG2, tau)),
                                             mul(z, incoming))))
    return sigma, tau, kappas, residual


def lawful_contexts(i, q):
    """Boundary indices in the definition of profile:scores, a=1,b=5."""
    variable = list(range(1, i)) + list(range(6-q, 6))
    options = [range(c+1) if c < i else range(i, c+1) for c in variable]
    for selected in product(*options):
        left, right = [None]*5, [None]*5
        for c, cut in zip(variable, selected):
            left[c-1] = right[c-1] = cut
        for c in range(i, 6-q):
            left[c-1], right[c-1] = i-1, i
        require(all(left[c-1] is not None and right[c-1] is not None
                    and 0 <= left[c-1] <= c and 0 <= right[c-1] <= c
                    for c in range(1, 6)), "lawful-context deadlines")
        yield left, right


def finite_context_example():
    """Exact finite geometry; no use of floating-point zero tolerances."""
    lengths = (1000, 1, 1, 1, 1, 100000)
    # Establish the continuous sigma(z) branch on the entire z interval,
    # not only at sampled bisection points. E_r(z) increases in z because
    # E_r'(z)=log(1+r/z)>0; (r+z)^(s-1) decreases when s<1.
    # These separated endpoint bounds therefore enclose the first equation
    # uniformly for 1<=z<=2, at s=1/2 and s=7/10 respectively.
    lower_first, upper_first = interval(0), interval(0)
    for r, L in zip(range(5, 0, -1), lengths[:5]):
        E_at_1 = mul(r+1, log_value(r+1))
        E_at_2 = sub(mul(r+2, log_value(r+2)), mul(2, LOG2))
        lower_first = add(lower_first, mul(L, sub(r, mul(
            exp_value(mul(F(-1, 2), log_value(r+1))), E_at_2))))
        upper_first = add(upper_first, mul(L, sub(r, mul(
            exp_value(mul(F(-3, 10), log_value(r+2))), E_at_1))))
    require(lower_first[0] > 0 > upper_first[1],
            "finite sigma(z) branch exists throughout 1 <= z <= 2")
    # The first equation has a strictly negative sigma derivative, so its
    # unique root is continuous in z. Opposite signs of the remaining
    # equation now give a face-critical point by the intermediate value
    # theorem. Every bisection step retains that certified sign bracket.
    lower, upper = F(1), F(2)
    require(finite_face_data(lower, lengths)[3][0] > 0 >
            finite_face_data(upper, lengths)[3][1], "finite face-root endpoints")
    for _ in range(36):
        midpoint = (lower + upper)/2
        residual = finite_face_data(midpoint, lengths)[3]
        if residual[0] > 0:
            lower = midpoint
        elif residual[1] < 0:
            upper = midpoint
        else:
            raise VerificationError("finite z root needs a tighter sigma enclosure")
    z = lower, upper
    sigma, tau, kappas, _ = finite_face_data(z, lengths)
    require(ALPHA2[1] < tau[0] < tau[1] < sigma[0] < sigma[1] < 1,
            "finite saddle has a strict plateau separation")
    V = [0] + list(accumulate(lengths))
    logz = log_value(z)
    entropies = [sub(mul(add(q, z), log_value(add(q, z))), mul(z, logz))
                 if q else interval(0) for q in range(6)]
    costs = [[sub(q, mul(kap, entropies[q])) for q in range(6)] for kap in kappas]
    # Full-plateau stationarity is an exact identity at the certified root.
    require(total(mul(L, costs[g][5-g]) for g, L in enumerate(lengths[:5]))[0] <= 0 <=
            total(mul(L, costs[g][5-g]) for g, L in enumerate(lengths[:5]))[1],
            "finite top-plateau balance enclosure")
    prefix = interval(0)
    incoming = interval(0)
    for j in range(1, 5):
        incoming = add(incoming, mul(lengths[j-1], kappas[j-1]))
        r = 6-j
        derivative = sub(V[j], mul(incoming, sub(
            mul(add(r, z), log_value(add(r, z))),
            mul(add(r-1, z), log_value(add(r-1, z))))))
        prefix = add(prefix, derivative)
        require(prefix[1] < 0, "finite ordered KKT prefix " + str(j))

    def reference(cuts):
        counts = tuple(sum(cut >= g for cut in cuts) for g in range(1, 6))
        require(all(0 <= count <= 6-g for g, count in enumerate(counts, 1)),
                "lawful context has an ineligible exterior name")
        if counts == (0, 0, 0, 0, 0) or counts == (5, 4, 3, 2, 1):
            return interval(0), counts  # exact zero / exact stationarity
        C = total(mul(L, costs[g][counts[g]]) for g, L in enumerate(lengths[:5]))
        require(C[0] > 0, "finite context reference nonnegativity")
        return C, counts

    pairs, edges = 0, 0
    for i in range(1, 6):
        for q in range(6-i):
            edge_pairs = 0
            for left, right in lawful_contexts(i, q):
                CL, _ = reference(left)
                CR, _ = reference(right)
                pairs += 1
                edge_pairs += 1
                if i == 1 and q == 0:
                    require(left == [0]*5 and right == [1]*5,
                            "critical split context pair")
                    require(2*V[max(left)] <= V[5] < 2*V[max(right)],
                            "critical left/right early/late classification")
                    require(CL == interval(0) and CR[0] > 0,
                            "critical split reference scores")
                    # Y=1+T*g_1(5)=1-lambda*sum_{r=1}^4 g_r(r).
                    balance = sub(0, total(costs[g][5-g] for g in range(1, 5)))
                    require(CR[0] <= balance[1] and balance[0] <= CR[1],
                            "critical split plateau-balance identity")
                    require(not any(0 < sum(v >= g for v in right) < 6-g
                                    and cut < g for g, cut in enumerate(right, 1)),
                            "critical right context has no good groups")
            require(edge_pairs > 0, "canonical edge has no lawful contexts")
            if i == 1 and q == 0:
                require(edge_pairs == 1, "critical split has one context pair")
            edges += 1
    require(edges == 15 and pairs == 540, "finite canonical-context inventory")
    require(sub(mul(power(LOG2, tau), LOG2), 1)[0] > 0,
            "finite terminal ordinary drift")
    return edges, pairs


def main():
    require(__debug__, "Run without -O: optimized Python is deliberately rejected.")
    comparable_example()
    print("PASS profile:example-two: exact saddle, endpoints, edge and suffix")
    split_example()
    print("PASS profile:example-split: exact separation and limiting saddle")
    edges, pairs = finite_context_example()
    print("PASS finite split geometry: " + str(edges) + " edges, " + str(pairs) +
          " lawful context pairs")


if __name__ == "__main__":
    try:
        main()
    except VerificationError as error:
        print("FAIL: " + str(error), file=sys.stderr)
        sys.exit(1)
