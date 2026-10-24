import config

def classical_factor(N:int)->dict :
    """
    Trial division for finding smallest prime factor of N
    """
    start = time.time()
    if N<2:
        raise ValueError("N must be >= 2")
    if N%2==0:
        return Result(
                N, 'classical',
                True, (2,N//2), None,
                time.time()-start, None,
                None, None, None,
                None, None,
                None, None,
                None
            ) #{"factors": (2, N//2), "method":"classical"}

    i = 3
    while i*i <= N:
        if N%i == 0:
            return Result (
                    N, 'classical',
                    True, (i,N//i), None,
                    time.time()-start, None,
                    None, None, None,
                    None, None,
                    None, None,
                    None
                ) #{"factors":(i, N//i), "method":"classical:trial_division"}
        i += 2
    return Result(
            N, 'classical',
            True, (1,N), None,
            time.time()-start, None,
            None, None, None,
            None, None,
            None, None,
            None
        ) #{"factors":(1, N), "method":"classical_prime"}
