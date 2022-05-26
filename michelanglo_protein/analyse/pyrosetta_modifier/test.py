# these are not real test...


    # m.output()


def paratest():
    import requests
    import time
    from multiprocessing_on_dill import (Pipe, Process)  # noqa
    # 1SFT/A/A/HIS`166
    pdbblock = requests.get('https://files.rcsb.org/download/1SFT.pdb').text
    kwargs = dict(pdbblock=pdbblock, target_resi=166, target_chain='A', cycles=1, radius=3)

    def subpro(child_conn, **kwargs):  # Pipe <- Union[dict, None]: # noqa
        try:
            print('started child')
            Mutator.reinit()
            mut = Mutator(**kwargs)
            data = mut.analyse_mutation('W')  # {ddG: float, scores: Dict[str, float], native:str, mutant:str, rmsd:int}
            print('done', len(data))
            child_conn.send(data)
            print('completed child')
        except BaseException as error:
            print('error child')
            child_conn.send({'error': f'{error.__class__.__name__}:{error}'})

    parent_conn, child_conn = Pipe()
    p = Process(target=subpro, args=(child_conn,), kwargs=kwargs, name='pyrosetta')
    p.start()

    while 1:
        time.sleep(5)
        print(parent_conn.poll())
        if parent_conn.poll():
            # p.terminate()
            break
        elif not p.is_alive():
            child_conn.send({'error': 'segmentation fault'})
            break
    msg = parent_conn.recv()
    print(msg)
    print('DONE!')
