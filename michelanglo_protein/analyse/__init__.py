try:
    from .Pymol_StructureAnalyser import StructureAnalyser
    from .pyrosetta_modifier import Mutator
except ModuleNotFoundError as err:
    import warnings
    warnings.warn(f'**{err.__class__.__name__}** {err}.')
    print('If you are not running Michelanglo in full this is fine.')
    print('Without pyrosetta VENUS will not complete. Without PyMOL editing of coordinates will not work.')
    StructureAnalyser = None
    Mutator = None