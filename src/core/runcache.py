"""
A per-run cache for the long experiments: with the environment variable EXP_CACHE set to a directory, a job decorated with `cached` stores
its result there under a key made of its file, name and arguments, and returns the stored result when the same job is asked again —
so an interrupted experiment restarts where it stopped. Without EXP_CACHE the jobs run as usual. The results are the same either way (the
jobs are deterministic in their arguments).
"""
import functools
import hashlib
import os
import pickle


def cached(fn):
    @functools.wraps(fn)
    def wrapper(args):
        d = os.environ.get('EXP_CACHE')
        if not d:
            return fn(args)
        src = os.path.splitext(os.path.basename(fn.__code__.co_filename))[0]      # the same whether run as a script or imported
        key = hashlib.sha1(repr((src, fn.__name__, args)).encode()).hexdigest()
        path = os.path.join(d, f'{src}.{fn.__name__}_{key}.pkl')
        if os.path.exists(path):
            with open(path, 'rb') as f:
                return pickle.load(f)
        out = fn(args)
        os.makedirs(d, exist_ok=True)
        tmp = f'{path}.tmp{os.getpid()}'
        with open(tmp, 'wb') as f:
            pickle.dump(out, f)
        os.replace(tmp, path)
        return out
    return wrapper
