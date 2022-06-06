import os, time
from multiprocessing_on_dill import Process, Pipe  # noqa
from multiprocessing_on_dill.connection import Connection
from typing import Dict, Any, Callable

def subprocess_task(fun: Callable, child_conn: Connection, options: Dict[str, Any]):  # Pipe <- Union[dict, None]:
    try:
        os.nice(19)
        data = fun(**options)
        child_conn.send(data)
    except BaseException as error:
        child_conn.send({'error': f'{error.__class__.__name__}:{error}'})


def run_subprocess(fun:Callable, timeout:float=60*10, **options) -> Any:
    """
    Run a function (fun) and wait for it.

    :param fun: unbound function.
    :return: whatever fun returns.
    """
    parent_conn, child_conn = Pipe()
    print(type(child_conn))
    p = Process(target=subprocess_task, args=(fun, child_conn, options))
    p.start()
    tick = time.time()
    while 1:
        if time.time() - tick > timeout:
            p.terminate()
            raise Exception(f'Timeout ({timeout} s) for {fun.__qualname__}')
        elif parent_conn.poll():
            break
        elif not p.is_alive():
            break
        else:
            pass
    return parent_conn.recv()
