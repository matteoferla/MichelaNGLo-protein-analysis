"""
This script got run once.. To make the images of the amino acids.

>>> AA().write('/Users/matteo/Coding/michelanglo/app/michelanglo_app/static/aa')

"""
from rdkit import Chem
from rdkit.Chem import AllChem, rdFMCS, rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Geometry.rdGeometry import Point2D
from typing import *

class AA:
    def __init__(self):
        self.aminoacids = self.read()

    def read(self):
        #found randomly online: https://binfalse.de/software/parseable-biodata/amino-acids/
        table = r'''Alanine Ala     A       C3H7NO2 O=C(O)C(N)C     nonpolar        neutral 1.8     nonessential    6.01    2.35    9.87
Arginine        Arg     R       C6H14N4O2       O=C(O)C(N)CCC/[NH+]=C(\[NH2])[NH2]  polar   positive        −4.5    essential       10.76   1.82    8.99
Asparagine      Asn     N       C4H8N2O3        O=C(N)CC(N)C(=O)O   polar   neutral −3.5    nonessential    5.41    2.14    8.72
Asparticacid    Asp     D       C4H7NO4 O=C(O)CC(N)C(=O)[O-]       polar   negative        −3.5    nonessential    2.85    1.99    9.90
Cysteine        Cys     C       C3H7NO2S        C(C(C(=O)O)N)S     nonpolar        neutral 2.5     nonessential    5.05    1.92    10.70
Glutamicacid    Glu     E       C5H9NO4 C(CC(=O)[O-])C(C(=O)O)N    polar   negative        −3.5    nonessential    3.15    2.10    9.47
Glutamine       Gln     Q       C5H10N2O3       O=C(N)CCC(N)C(=O)O      polar   neutral −3.5    nonessential    5.65    2.17    9.13
Glycine Gly     G       C2H5NO2 C(C(=O)O)N      nonpolar        neutral −0.4    nonessential    6.06    2.35    9.78
Histidine       His     H       C6H9N3O2        O=C(C(CC1=CNC=N1)N)O polar   neutral(90%)−3.2        essential       7.60    1.80    9.33
Isoleucine      Ile     I       C6H13NO2        CC[C@H](C)C(C(=O)O)N       nonpolar        neutral 4.5     essential       6.05    2.32    9.76
Leucine Leu     L       C6H13NO2        CC(C)CC(C(=O)O)N   nonpolar        neutral 3.8     essential       6.01    2.33    9.74
Lysine  Lys     K       C6H14N2O2       C(CC[NH3+])CC(C(=O)O)N       polar   positive        −3.9    essential       9.60    2.16    9.06
Methionine      Met     M       C5H11NO2S       CSCCC(C(=O)O)N  nonpolar        neutral 1.9     essential       5.74    2.13    9.28
Phenylalanine   Phe     F       C9H11NO2        c1ccc(cc1)CC(C(=O)O)N      nonpolar        neutral 2.8     essential       5.49    2.20    9.31
Proline Pro     P       C5H9NO2 C1CC(NC1)C(=O)O nonpolar        neutral −1.6    nonessential    6.30    1.95    10.64
Serine  Ser     S       C3H7NO3 C(C(C(=O)O)N)O     polar   neutral −0.8    nonessential    5.68    2.19    9.21
Threonine       Thr     T       C4H9NO3 C[C@H](C(C(=O)O)N)O        polar   neutral −0.7    essential       2.09    2.09    9.10
Tryptophan      Trp     W       C11H12N2O2      c1ccc2c(c1)c(c[nH]2)CC(C(=O)O)N    nonpolar        neutral −0.9    essential       5.89    2.46    9.41
Tyrosine        Tyr     Y       C9H11NO3        NC(Cc1ccc(O)cc1)C(O)=O     polar   neutral −1.3    nonessential    5.64    2.20    9.21
Valine  Val     V       C5H11NO2        CC(C)C(C(=O)O)N    nonpolar        neutral 4.2     essential       6.00    2.39    9.74'''
        aminoacids = {} # single letter to mol.
        for aa in table.split('\n'):
            a = aa.split()
            mol = Chem.MolFromSmiles(a[4])
            AllChem.Compute2DCoords(mol)
            mol.SetProp('_Name', a[0])
            if a[0] == 'Proline':
                self.fix_carboxy(mol)
            else:
                self.fix_backbone(mol)
            aminoacids[a[2]] = mol
        ref = Chem.MolFromSmiles('NCC(=O)O')
        AllChem.Compute2DCoords(ref, coordMap={0: Point2D(0,0), 1: Point2D(0,1.34),})
        for name1, mol in aminoacids.items():
            AllChem.GenerateDepictionMatching2DStructure(mol, ref)
        return aminoacids

    def fix_carboxy(self, mol: Chem.Mol):
        carbo = Chem.MolFromSmiles('C(=O)O')
        if mol.HasSubstructMatch(carbo):
            O = mol.GetSubstructMatch(carbo)[2]
            mol.GetAtomWithIdx(O).SetFormalCharge(-1)

    def fix_amine(self, mol: Chem.Mol):
        amine = Chem.MolFromSmiles('CN([H])([H])')
        if mol.HasSubstructMatch(amine):
            N = mol.GetSubstructMatch(amine)[1]
            mol.GetAtomWithIdx(N).SetFormalCharge(+1)

    def fix_backbone(self, mol: Chem.Mol):
        N, Ca, C, O1, O2 = self.get_backbone(mol)
        mol.GetAtomWithIdx(N).SetFormalCharge(+1)
        mol.GetAtomWithIdx(O2).SetFormalCharge(-1)

    def get_backbone(self, mol: Chem.Mol) -> Tuple[int, int, int, int, int]:
        gly = Chem.MolFromSmiles('NCC(=O)O')
        return mol.GetSubstructMatch(gly)

    def get_Calpha(self, mol: Chem.Mol) -> int:
        return self.get_backbone(mol)[1]

    def get_svg(self, ab: str, ad: str) -> str:
        def mod(letter: str):
            # self.aminoacids[k] is a Chem.Mol, but let\s copy it.
            mol = Chem.Mol(self.aminoacids[letter])
            C = self.get_Calpha(mol)
            # to cheat with selection let's make Calpha isotope 14.
            mol.GetAtomWithIdx(C).SetAtomicNum(0) #.SetIsotope(14)
            return mol
        #mols = [mod(k) for k in (ab, ad)]
        mols = [self.aminoacids[k] for k in (ab, ad)]
        # mod atoms
        if 'P' not in (ab, ad):
            res=Chem.rdFMCS.FindMCS(mols,
                                    #ringMatchesRingOnly=True,
                                    #completeRingsOnly=True,
                                    matchValences=True,
                                    #matchChiralTag=True,
                                    bondCompare=Chem.rdFMCS.BondCompare.CompareOrderExact
                                    #atomCompare=Chem.rdFMCS.AtomCompare.CompareIsotopes
                                   )
        else:
            res=Chem.rdFMCS.FindMCS(mols,
                                    ringMatchesRingOnly=True,
                                    completeRingsOnly=True,
                                    matchValences=True,
                                    #matchChiralTag=True,
                                    bondCompare=Chem.rdFMCS.BondCompare.CompareOrderExact
                                    #atomCompare=Chem.rdFMCS.AtomCompare.CompareIsotopes
                                   )
        common = Chem.MolFromSmarts(res.smartsString)
        drawer = rdMolDraw2D.MolDraw2DSVG(400,200)
        # remove the Calpha C-14 dodginess
        #mols[0].GetAtomsMatchingQuery(Chem.rdqueries.IsotopeEqualsQueryAtom(14))[0].SetIsotope(0)
        #mols[0].GetAtomsMatchingQuery(Chem.rdqueries.AtomNumEqualsQueryAtom(0))[0].SetAtomicNum(6)
        #rdDepictor.Compute2DCoords(mols[0])
        mol = mols[0]
        match = mol.GetSubstructMatch(common)
        backbone = self.get_backbone(mol)
        unconserved = [i for i in range(mol.GetNumAtoms()) if i not in match and i not in backbone]
        coral = (1, 0.498, 0.314)
        robin_egg = (0.122, 0.808, 0.796)
        pale_turquoise = (0.65, 0.90, 0.88)
        azure_mist = (0.941, 1., 1.)
        atom_highlights = {**{atom_idx: [coral] for atom_idx in unconserved},
                           **{atom_idx: [pale_turquoise] for atom_idx in backbone}}
        bond_filter = lambda bond, atoms: bond.GetBeginAtomIdx() in atoms and bond.GetEndAtomIdx() in atoms
        bond_highlights = {**{bond.GetIdx(): [coral] for bond in mol.GetBonds() if bond_filter(bond, unconserved)},
                          **{bond.GetIdx(): [pale_turquoise] for bond in mol.GetBonds() if bond_filter(bond,   backbone)}}
        drawer.DrawMoleculeWithHighlights(mol,
                                          '', atom_highlights, bond_highlights,
                                          {},{})
        drawer.FinishDrawing()
        return drawer.GetDrawingText()

    def draw(self, ab: str, ad: str):
        from IPython.display import display, SVG
        display(SVG(self.get_svg(ab, ad)))

    def write(self, folder:str):
        """
        makes and writes the images of the amino acids.
        The format is '{ab}{ad}.svg' where ab is from/whence, and ad is to/whither.
        The match matches only complete rings ``FindMCS(mols, ringMatchesRingOnly=True, completeRingsOnly=True)``,
        because it was matching any bonding and making a weird open match with Y and W.

        :return: None
        """
        for ab in self.aminoacids:
            for ad in self.aminoacids:
                svg = self.get_svg(ab, ad)
                with open(f'{folder}/{ab}{ad}.svg', 'w') as fh:
                    fh.write(svg)

    def sizes(self):
        """
        Gets the number of heavy atoms that a AA ought to have.

        :return:
        """
        return {aa: self.aminoacids[aa].GetNumAtoms() for aa in self.aminoacids}