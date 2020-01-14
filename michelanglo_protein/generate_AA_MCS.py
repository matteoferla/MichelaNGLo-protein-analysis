"""
This script got run once.. To make the images of the amino acids.

"""
from rdkit import Chem
from rdkit.Chem import AllChem, rdFMCS, rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D

class AA:
    def __init__(self):
        self.aminoacids = self.read()

    def read(self):
        #found randomly online: https://binfalse.de/software/parseable-biodata/amino-acids/
        table = r'''Alanine Ala     A       C3H7NO2 O=C(O)C(N)C     nonpolar        neutral 1.8     nonessential    6.01    2.35    9.87
Arginine        Arg     R       C6H14N4O2       O=C(O)C(N)CCC/N=C(\N)N  polar   positive        −4.5    essential       10.76   1.82    8.99
Asparagine      Asn     N       C4H8N2O3        O=C(N)C[C@H](N)C(=O)O   polar   neutral −3.5    nonessential    5.41    2.14    8.72
Asparticacid    Asp     D       C4H7NO4 O=C(O)CC(N)C(=O)O       polar   negative        −3.5    nonessential    2.85    1.99    9.90
Cysteine        Cys     C       C3H7NO2S        C([C@@H](C(=O)O)N)S     nonpolar        neutral 2.5     nonessential    5.05    1.92    10.70
Glutamicacid    Glu     E       C5H9NO4 C(CC(=O)O)C(C(=O)O)N    polar   negative        −3.5    nonessential    3.15    2.10    9.47
Glutamine       Gln     Q       C5H10N2O3       O=C(N)CCC(N)C(=O)O      polar   neutral −3.5    nonessential    5.65    2.17    9.13
Glycine Gly     G       C2H5NO2 C(C(=O)O)N      nonpolar        neutral −0.4    nonessential    6.06    2.35    9.78
Histidine       His     H       C6H9N3O2        O=C([C@H](CC1=CNC=N1)N)O polar   neutral(90%)−3.2        essential       7.60    1.80    9.33
Isoleucine      Ile     I       C6H13NO2        CC[C@H](C)[C@@H](C(=O)O)N       nonpolar        neutral 4.5     essential       6.05    2.32    9.76
Leucine Leu     L       C6H13NO2        CC(C)C[C@@H](C(=O)O)N   nonpolar        neutral 3.8     essential       6.01    2.33    9.74
Lysine  Lys     K       C6H14N2O2       C(CCN)CC(C(=O)O)N       polar   positive        −3.9    essential       9.60    2.16    9.06
Methionine      Met     M       C5H11NO2S       CSCCC(C(=O)O)N  nonpolar        neutral 1.9     essential       5.74    2.13    9.28
Phenylalanine   Phe     F       C9H11NO2        c1ccc(cc1)C[C@@H](C(=O)O)N      nonpolar        neutral 2.8     essential       5.49    2.20    9.31
Proline Pro     P       C5H9NO2 C1CC(NC1)C(=O)O nonpolar        neutral −1.6    nonessential    6.30    1.95    10.64
Serine  Ser     S       C3H7NO3 C([C@@H](C(=O)O)N)O     polar   neutral −0.8    nonessential    5.68    2.19    9.21
Threonine       Thr     T       C4H9NO3 C[C@H]([C@@H](C(=O)O)N)O        polar   neutral −0.7    essential       2.09    2.09    9.10
Tryptophan      Trp     W       C11H12N2O2      c1ccc2c(c1)c(c[nH]2)C[C@@H](C(=O)O)N    nonpolar        neutral −0.9    essential       5.89    2.46    9.41
Tyrosine        Tyr     Y       C9H11NO3        N[C@@H](Cc1ccc(O)cc1)C(O)=O     polar   neutral −1.3    nonessential    5.64    2.20    9.21
Valine  Val     V       C5H11NO2        CC(C)[C@@H](C(=O)O)N    nonpolar        neutral 4.2     essential       6.00    2.39    9.74'''
        aminoacids = {} # single letter to mol.
        for aa in table.split('\n'):
            a = aa.split()
            mol = Chem.MolFromSmiles(a[4])
            aminoacids[a[2]] = mol
        return aminoacids

    def draw(self):
        for ab in self.aminoacids:
            for ad in self.aminoacids:
                mols = [self.aminoacids[k] for k in (ab, ad)]
                res=Chem.rdFMCS.FindMCS(mols)
                common = Chem.MolFromSmarts(res.smartsString)
                drawer = rdMolDraw2D.MolDraw2DSVG(400,200)
                rdDepictor.Compute2DCoords(mols[0])
                inv = [i for i in range(mols[0].GetNumAtoms()) if i not in mols[0].GetSubstructMatch(common)]
                drawer.DrawMolecule(mols[0], highlightAtoms=inv)
                drawer.FinishDrawing()
                svg = drawer.GetDrawingText()
                open(f'{ab}{ad}.svg','w').write(svg)

    def sizes(self):
        return {aa: self.aminoacids[aa].GetNumAtoms() for aa in self.aminoacids}



if __name__ == '__main__':
    AA().draw()
    #print(AA().sizes())