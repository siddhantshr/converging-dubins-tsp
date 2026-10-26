"""Optional process pool, for `use_all_cores: true` in config.yml.

LKH (elkai) holds the GIL and keeps global state in C, so the parallel work
runs in worker processes, not threads. Two things use the pool:

* edge costs: the uncached Dubins lengths (upper bound) and interval
  problems (lower bound) are computed in chunks, one chunk per worker task;
* LKH: the upper-bound solve runs in a worker while this process builds the
  lower-bound graph (see bounds.compute_bounds).

The results are the same with or without the pool; only the time changes.
With `pool=None` every helper here runs inline, in this process.
"""

from __future__ import annotations

import contextlib
import multiprocessing
import os
from concurrent.futures import Future, ProcessPoolExecutor

WORKERS = os.cpu_count() or 1
MIN_PARALLEL_ITEMS = 2000  # below this, sending work to the workers costs more than it saves
TASKS_PER_WORKER = 4  # chunks per worker, so a slow chunk doesn't hold up the others


def pool(enabled: bool):
    """Context manager: a process pool with one worker per core, or None."""
    if not enabled or WORKERS < 2:
        return contextlib.nullcontext(None)
    # spawn, not fork: forking a process that has numpy / LKH loaded is unsafe on macOS
    return ProcessPoolExecutor(WORKERS, mp_context=multiprocessing.get_context("spawn"))


def submit(executor: ProcessPoolExecutor | None, fn, *args) -> Future:
    """fn(*args) in a worker, or right now in this process if there is no pool."""
    if executor is not None:
        return executor.submit(fn, *args)
    done: Future = Future()
    try:
        done.set_result(fn(*args))
    except Exception as err:  # re-raised by done.result(), as a pool would
        done.set_exception(err)
    return done


def map_chunked(executor: ProcessPoolExecutor | None, fn, items: list, *args) -> list:
    """fn(*args, chunk) returns one result per item of chunk; results in item order."""
    if executor is None or len(items) < MIN_PARALLEL_ITEMS:
        return fn(*args, items)
    size = -(-len(items) // (WORKERS * TASKS_PER_WORKER))
    futures = [executor.submit(fn, *args, items[k : k + size]) for k in range(0, len(items), size)]
    return [result for future in futures for result in future.result()]
