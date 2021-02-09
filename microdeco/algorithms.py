# =============================================================================

# 

# =============================================================================



def eq_gf_limit(gf, p_n2, p_he, a_n2, b_n2, a_he, b_he):

    """

    Calculate ascent ceiling limit of a tissue compartment using Buhlmann

    equation extended with gradient factors by Erik Baker.



    The returned value is absolute pressure of depth of the ascent ceiling.



    :param gf: Gradient factor value.

    :param p_n2: Current tissue pressure for nitrogen.

    :param p_he: Current tissue pressure for helium.

    :param a_n2: Nitrox Buhlmann coefficient A.

    :param b_n2: Nitrox Buhlmann coefficient B.

    :param a_he: Helium Buhlmann coefficient A.

    :param b_he: Helium Buhlmann coefficient B.

    """

    assert gf > 0 and gf <= 1.5

    p = p_n2 + p_he

    a = (a_n2 * p_n2 + a_he * p_he) / p

    b = (b_n2 * p_n2 + b_he * p_he) / p

    return (p - a * gf) / (gf / b + 1 - gf)





def recurse_while(predicate, f, *args):

    """

    Accumulate value by executing recursively function `f`.



    The function `f` is executed with starting arguments. While the

    predicate for the result is true, the result is fed into function `f`.



    If predicate is never true then starting arguments are returned.



    :param predicate: Predicate function guarding execution.

    :param f: Function to execute.

    :param *args: Starting arguments.

    """

    result = f(*args)

    result = result if type(result) == tuple else (result, )

    while predicate(*result):

        args = result # predicate(args) is always true

        result = f(*args)

        result = result if type(result) == tuple else (result, )



    return args if len(args) > 1 else args[0]





def bisect_find(n, f, *args, **kw):

    """

    Find largest `k` for which `f(k)` is true.



    The k is integer in range 1 <= k <= n.  If there is no `k` for which

    `f(k)` is true, then return `0`.



    :param n: Range for `k`, so :math:`1 <= k <= n`.

    :param f: Invariant function accepting `k`.

    :param *args: Additional positional parameters of `f`.

    :param **kw: Additional named parameters of `f`.

    """

    lo = 1

    hi = n + 1



    while lo < hi:

        k = (lo + hi) // 2



        if f(k, *args, **kw):

            lo = k + 1

        else:

            hi = k



    return hi - 1 # hi is first k for which f(k) is not true, so f(hi - 1) is true