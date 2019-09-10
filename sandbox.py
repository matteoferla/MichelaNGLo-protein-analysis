from protein import ProteinAnalyser, ProteinCore, Mutation
from protein.settings_handler import global_settings
from protein.generate import ProteinGatherer, ProteomeGatherer
from protein.protein_analysis import StructureAnalyser

import pickle


def test_ProteinAnalyser():
    p = ProteinAnalyser(uniprot = 'O75015').gload()
    print(p)
    #p.mutation = Mutation('p.I124P')

    structure = '''HEADER    IMMUNE SYSTEM                           07-JUL-00   1E4J              
    TITLE     CRYSTAL STRUCTURE OF THE SOLUBLE HUMAN FC-GAMMA                       
    TITLE    2 RECEPTOR III                                                         
    COMPND    MOL_ID: 1;                                                            
    COMPND   2 MOLECULE: LOW AFFINITY IMMUNOGLOBULIN GAMMA FC RECEPTOR              
    COMPND   3  III;                                                                
    COMPND   4 CHAIN: A;                                                            
    COMPND   5 FRAGMENT: EXTRACELLULAR DOMAIN;                                      
    COMPND   6 SYNONYM: CD16;                                                       
    COMPND   7 ENGINEERED: YES                                                      
    SOURCE    MOL_ID: 1;                                                            
    SOURCE   2 ORGANISM_SCIENTIFIC: HOMO SAPIENS;                                   
    SOURCE   3 ORGANISM_COMMON: HUMAN;                                              
    SOURCE   4 ORGANISM_TAXID: 9606;                                                
    SOURCE   5 CELL: LEUKOCYTE;                                                     
    SOURCE   6 EXPRESSION_SYSTEM: ESCHERICHIA COLI;                                 
    SOURCE   7 EXPRESSION_SYSTEM_TAXID: 469008;                                     
    SOURCE   8 EXPRESSION_SYSTEM_STRAIN: BL21(DE3);                                 
    SOURCE   9 EXPRESSION_SYSTEM_CELLULAR_LOCATION: INCLUSION BODIES;               
    SOURCE  10 EXPRESSION_SYSTEM_PLASMID: PET21                                     
    KEYWDS    IMMUNE SYSTEM, IGG, FC, RECEPTOR, CD16, GAMMA                         
    EXPDTA    X-RAY DIFFRACTION                                                     
    AUTHOR    P.SONDERMANN,R.HUBER,U.JACOB                                          
    REVDAT   2   24-FEB-09 1E4J    1       VERSN                                    
    REVDAT   1   04-AUG-00 1E4J    0                                                
    JRNL        AUTH   P.SONDERMANN,R.HUBER,V.OOSTHUIZEN,U.JACOB                    
    JRNL        TITL   THE 3.2-A CRYSTAL STRUCTURE OF THE HUMAN IGG1 FC             
    JRNL        TITL 2 FRAGMENT-FC GAMMARIII COMPLEX.                               
    JRNL        REF    NATURE                        V. 406   267 2000              
    JRNL        REFN                   ISSN 0028-0836                               
    JRNL        PMID   10917521                                                     
    JRNL        DOI    10.1038/35018508                                             
    REMARK   2                                                                      
    REMARK   2 RESOLUTION.    2.5  ANGSTROMS.                                       
    REMARK   3                                                                      
    REMARK   3 REFINEMENT.                                                          
    REMARK   3   PROGRAM     : CNS 0.9                                              
    REMARK   3   AUTHORS     : BRUNGER,ADAMS,CLORE,DELANO,GROS,GROSSE-              
    REMARK   3               : KUNSTLEVE,JIANG,KUSZEWSKI,NILGES, PANNU,             
    REMARK   3               : READ,RICE,SIMONSON,WARREN                            
    REMARK   3                                                                      
    REMARK   3  REFINEMENT TARGET : NULL                                            
    REMARK   3                                                                      
    REMARK   3  DATA USED IN REFINEMENT.                                            
    REMARK   3   RESOLUTION RANGE HIGH (ANGSTROMS) : 2.5                            
    REMARK   3   RESOLUTION RANGE LOW  (ANGSTROMS) : 50                             
    REMARK   3   DATA CUTOFF            (SIGMA(F)) : NULL                           
    REMARK   3   DATA CUTOFF HIGH         (ABS(F)) : NULL                           
    REMARK   3   DATA CUTOFF LOW          (ABS(F)) : NULL                           
    REMARK   3   COMPLETENESS (WORKING+TEST)   (%) : 94.7                           
    REMARK   3   NUMBER OF REFLECTIONS             : 6725                           
    REMARK   3                                                                      
    REMARK   3  FIT TO DATA USED IN REFINEMENT.                                     
    REMARK   3   CROSS-VALIDATION METHOD          : THROUGHOUT                      
    REMARK   3   FREE R VALUE TEST SET SELECTION  : RANDOM                          
    REMARK   3   R VALUE            (WORKING SET) : 0.1951                          
    REMARK   3   FREE R VALUE                     : 0.2614                          
    REMARK   3   FREE R VALUE TEST SET SIZE   (%) : 5.0                             
    REMARK   3   FREE R VALUE TEST SET COUNT      : NULL                            
    REMARK   3   ESTIMATED ERROR OF FREE R VALUE  : NULL                            
    REMARK   3                                                                      
    REMARK   3  FIT IN THE HIGHEST RESOLUTION BIN.                                  
    REMARK   3   TOTAL NUMBER OF BINS USED           : NULL                         
    REMARK   3   BIN RESOLUTION RANGE HIGH       (A) : NULL                         
    REMARK   3   BIN RESOLUTION RANGE LOW        (A) : NULL                         
    REMARK   3   BIN COMPLETENESS (WORKING+TEST) (%) : 72.3                         
    REMARK   3   REFLECTIONS IN BIN    (WORKING SET) : NULL                         
    REMARK   3   BIN R VALUE           (WORKING SET) : NULL                         
    REMARK   3   BIN FREE R VALUE                    : NULL                         
    REMARK   3   BIN FREE R VALUE TEST SET SIZE  (%) : NULL                         
    REMARK   3   BIN FREE R VALUE TEST SET COUNT     : NULL                         
    REMARK   3   ESTIMATED ERROR OF BIN FREE R VALUE : NULL                         
    REMARK   3                                                                      
    REMARK   3  NUMBER OF NON-HYDROGEN ATOMS USED IN REFINEMENT.                    
    REMARK   3   PROTEIN ATOMS            : 1376                                    
    REMARK   3   NUCLEIC ACID ATOMS       : 0                                       
    REMARK   3   HETEROGEN ATOMS          : 0                                       
    REMARK   3   SOLVENT ATOMS            : 67                                      
    REMARK   3                                                                      
    REMARK   3  B VALUES.                                                           
    REMARK   3   FROM WILSON PLOT           (A**2) : NULL                           
    REMARK   3   MEAN B VALUE      (OVERALL, A**2) : 29.7                           
    REMARK   3   OVERALL ANISOTROPIC B VALUE.                                       
    REMARK   3    B11 (A**2) : 0.521                                                
    REMARK   3    B22 (A**2) : -6.305                                               
    REMARK   3    B33 (A**2) : 5.784                                                
    REMARK   3    B12 (A**2) : 0                                                    
    REMARK   3    B13 (A**2) : 0                                                    
    REMARK   3    B23 (A**2) : 0                                                    
    REMARK   3                                                                      
    REMARK   3  ESTIMATED COORDINATE ERROR.                                         
    REMARK   3   ESD FROM LUZZATI PLOT        (A) : NULL                            
    REMARK   3   ESD FROM SIGMAA              (A) : NULL                            
    REMARK   3   LOW RESOLUTION CUTOFF        (A) : NULL                            
    REMARK   3                                                                      
    REMARK   3  CROSS-VALIDATED ESTIMATED COORDINATE ERROR.                         
    REMARK   3   ESD FROM C-V LUZZATI PLOT    (A) : NULL                            
    REMARK   3   ESD FROM C-V SIGMAA          (A) : NULL                            
    REMARK   3                                                                      
    REMARK   3  RMS DEVIATIONS FROM IDEAL VALUES.                                   
    REMARK   3   BOND LENGTHS                 (A) : 0.012                           
    REMARK   3   BOND ANGLES            (DEGREES) : 1.64                            
    REMARK   3   DIHEDRAL ANGLES        (DEGREES) : NULL                            
    REMARK   3   IMPROPER ANGLES        (DEGREES) : NULL                            
    REMARK   3                                                                      
    REMARK   3  ISOTROPIC THERMAL MODEL : NULL                                      
    REMARK   3                                                                      
    REMARK   3  ISOTROPIC THERMAL FACTOR RESTRAINTS.    RMS    SIGMA                
    REMARK   3   MAIN-CHAIN BOND              (A**2) : NULL  ; NULL                 
    REMARK   3   MAIN-CHAIN ANGLE             (A**2) : NULL  ; NULL                 
    REMARK   3   SIDE-CHAIN BOND              (A**2) : NULL  ; NULL                 
    REMARK   3   SIDE-CHAIN ANGLE             (A**2) : NULL  ; NULL                 
    REMARK   3                                                                      
    REMARK   3  BULK SOLVENT MODELING.                                              
    REMARK   3   METHOD USED : DENSITY MODIFICATION                                 
    REMARK   3   KSOL        : NULL                                                 
    REMARK   3   BSOL        : NULL                                                 
    REMARK   3                                                                      
    REMARK   3  NCS MODEL : NULL                                                    
    REMARK   3                                                                      
    REMARK   3  NCS RESTRAINTS.                         RMS   SIGMA/WEIGHT          
    REMARK   3   GROUP  1  POSITIONAL            (A) : NULL  ; NULL                 
    REMARK   3   GROUP  1  B-FACTOR           (A**2) : NULL  ; NULL                 
    REMARK   3                                                                      
    REMARK   3  PARAMETER FILE  1  : PROTEIN.PARAM                                  
    REMARK   3  PARAMETER FILE  2  : WATER.PARAM                                    
    REMARK   3  PARAMETER FILE  3  : NULL                                           
    REMARK   3  TOPOLOGY FILE  1   : NULL                                           
    REMARK   3  TOPOLOGY FILE  2   : NULL                                           
    REMARK   3  TOPOLOGY FILE  3   : NULL                                           
    REMARK   3                                                                      
    REMARK   3  OTHER REFINEMENT REMARKS: NULL                                      
    REMARK   4                                                                      
    REMARK   4 1E4J COMPLIES WITH FORMAT V. 3.15, 01-DEC-08                         
    REMARK 100                                                                      
    REMARK 100 THIS ENTRY HAS BEEN PROCESSED BY PDBE ON 10-JUL-00.                  
    REMARK 100 THE PDBE ID CODE IS EBI-5144.                                        
    REMARK 200                                                                      
    REMARK 200 EXPERIMENTAL DETAILS                                                 
    REMARK 200  EXPERIMENT TYPE                : X-RAY DIFFRACTION                  
    REMARK 200  DATE OF DATA COLLECTION        : 15-SEP-99                          
    REMARK 200  TEMPERATURE           (KELVIN) : 291                                
    REMARK 200  PH                             : 8.0                                
    REMARK 200  NUMBER OF CRYSTALS USED        : 1                                  
    REMARK 200                                                                      
    REMARK 200  SYNCHROTRON              (Y/N) : N                                  
    REMARK 200  RADIATION SOURCE               : ROTATING ANODE                     
    REMARK 200  BEAMLINE                       : NULL                               
    REMARK 200  X-RAY GENERATOR MODEL          : RIGAKU RU200                       
    REMARK 200  MONOCHROMATIC OR LAUE    (M/L) : M                                  
    REMARK 200  WAVELENGTH OR RANGE        (A) : 1.5418                             
    REMARK 200  MONOCHROMATOR                  : NI FILTER                          
    REMARK 200  OPTICS                         : NULL                               
    REMARK 200                                                                      
    REMARK 200  DETECTOR TYPE                  : IMAGE PLATE                        
    REMARK 200  DETECTOR MANUFACTURER          : MARRESEARCH                        
    REMARK 200  INTENSITY-INTEGRATION SOFTWARE : MOSFLM                             
    REMARK 200  DATA SCALING SOFTWARE          : SCALA                              
    REMARK 200                                                                      
    REMARK 200  NUMBER OF UNIQUE REFLECTIONS   : 6725                               
    REMARK 200  RESOLUTION RANGE HIGH      (A) : 2.5                                
    REMARK 200  RESOLUTION RANGE LOW       (A) : 40                                 
    REMARK 200  REJECTION CRITERIA  (SIGMA(I)) : 2                                  
    REMARK 200                                                                      
    REMARK 200 OVERALL.                                                             
    REMARK 200  COMPLETENESS FOR RANGE     (%) : 94.7                               
    REMARK 200  DATA REDUNDANCY                : 3.0                                
    REMARK 200  R MERGE                    (I) : 0.114                              
    REMARK 200  R SYM                      (I) : NONE                               
    REMARK 200  <I/SIGMA(I)> FOR THE DATA SET  : NULL                               
    REMARK 200                                                                      
    REMARK 200 IN THE HIGHEST RESOLUTION SHELL.                                     
    REMARK 200  HIGHEST RESOLUTION SHELL, RANGE HIGH (A) : 2.5                      
    REMARK 200  HIGHEST RESOLUTION SHELL, RANGE LOW  (A) : 2.3                      
    REMARK 200  COMPLETENESS FOR SHELL     (%) : 50.0                               
    REMARK 200  DATA REDUNDANCY IN SHELL       : NULL                               
    REMARK 200  R MERGE FOR SHELL          (I) : NULL                               
    REMARK 200  R SYM FOR SHELL            (I) : NULL                               
    REMARK 200  <I/SIGMA(I)> FOR SHELL         : NULL                               
    REMARK 200                                                                      
    REMARK 200 DIFFRACTION PROTOCOL: SINGLE WAVELENGTH                              
    REMARK 200 METHOD USED TO DETERMINE THE STRUCTURE: MOLECULAR REPLACEMENT        
    REMARK 200 SOFTWARE USED: AMORE                                                 
    REMARK 200 STARTING MODEL: 2FCB                                                 
    REMARK 200                                                                      
    REMARK 200 REMARK: NULL                                                         
    REMARK 280                                                                      
    REMARK 280 CRYSTAL                                                              
    REMARK 280 SOLVENT CONTENT, VS  (%): 48                                         
    REMARK 280 MATTHEWS COEFFICIENT, VM (ANGSTROMS**3/DA): 2.37                     
    REMARK 280                                                                      
    REMARK 280 CRYSTALLIZATION CONDITIONS: 0.1M MES/TRIS PH 7.8, 22.0% PEG8000      
    REMARK 290                                                                      
    REMARK 290 CRYSTALLOGRAPHIC SYMMETRY                                            
    REMARK 290 SYMMETRY OPERATORS FOR SPACE GROUP: P 2 21 21                        
    REMARK 290                                                                      
    REMARK 290      SYMOP   SYMMETRY                                                
    REMARK 290     NNNMMM   OPERATOR                                                
    REMARK 290       1555   X,Y,Z                                                   
    REMARK 290       2555   X,-Y,-Z                                                 
    REMARK 290       3555   -X,Y+1/2,-Z+1/2                                         
    REMARK 290       4555   -X,-Y+1/2,Z+1/2                                         
    REMARK 290                                                                      
    REMARK 290     WHERE NNN -> OPERATOR NUMBER                                     
    REMARK 290           MMM -> TRANSLATION VECTOR                                  
    REMARK 290                                                                      
    REMARK 290 CRYSTALLOGRAPHIC SYMMETRY TRANSFORMATIONS                            
    REMARK 290 THE FOLLOWING TRANSFORMATIONS OPERATE ON THE ATOM/HETATM             
    REMARK 290 RECORDS IN THIS ENTRY TO PRODUCE CRYSTALLOGRAPHICALLY                
    REMARK 290 RELATED MOLECULES.                                                   
    REMARK 290   SMTRY1   1  1.000000  0.000000  0.000000        0.00000            
    REMARK 290   SMTRY2   1  0.000000  1.000000  0.000000        0.00000            
    REMARK 290   SMTRY3   1  0.000000  0.000000  1.000000        0.00000            
    REMARK 290   SMTRY1   2  1.000000  0.000000  0.000000        0.00000            
    REMARK 290   SMTRY2   2  0.000000 -1.000000  0.000000        0.00000            
    REMARK 290   SMTRY3   2  0.000000  0.000000 -1.000000        0.00000            
    REMARK 290   SMTRY1   3 -1.000000  0.000000  0.000000        0.00000            
    REMARK 290   SMTRY2   3  0.000000  1.000000  0.000000       30.14500            
    REMARK 290   SMTRY3   3  0.000000  0.000000 -1.000000       42.79000            
    REMARK 290   SMTRY1   4 -1.000000  0.000000  0.000000        0.00000            
    REMARK 290   SMTRY2   4  0.000000 -1.000000  0.000000       30.14500            
    REMARK 290   SMTRY3   4  0.000000  0.000000  1.000000       42.79000            
    REMARK 290                                                                      
    REMARK 290 REMARK: NULL                                                         
    REMARK 300                                                                      
    REMARK 300 BIOMOLECULE: 1                                                       
    REMARK 300 SEE REMARK 350 FOR THE AUTHOR PROVIDED AND/OR PROGRAM                
    REMARK 300 GENERATED ASSEMBLY INFORMATION FOR THE STRUCTURE IN                  
    REMARK 300 THIS ENTRY.  THE REMARK MAY ALSO PROVIDE INFORMATION ON              
    REMARK 300 BURIED SURFACE AREA.                                                 
    REMARK 350                                                                      
    REMARK 350 GENERATING THE BIOMOLECULE                                           
    REMARK 350 COORDINATES FOR A COMPLETE MULTIMER REPRESENTING THE KNOWN           
    REMARK 350 BIOLOGICALLY SIGNIFICANT OLIGOMERIZATION STATE OF THE                
    REMARK 350 MOLECULE CAN BE GENERATED BY APPLYING BIOMT TRANSFORMATIONS          
    REMARK 350 GIVEN BELOW.  BOTH NON-CRYSTALLOGRAPHIC AND                          
    REMARK 350 CRYSTALLOGRAPHIC OPERATIONS ARE GIVEN.                               
    REMARK 350                                                                      
    REMARK 350 BIOMOLECULE:  1                                                      
    REMARK 350 AUTHOR DETERMINED BIOLOGICAL UNIT: MONOMERIC                         
    REMARK 350 SOFTWARE DETERMINED QUATERNARY STRUCTURE: MONOMERIC                  
    REMARK 350 SOFTWARE USED: PQS                                                   
    REMARK 350 APPLY THE FOLLOWING TO CHAINS: A                                     
    REMARK 350   BIOMT1   1  1.000000  0.000000  0.000000        0.00000            
    REMARK 350   BIOMT2   1  0.000000  1.000000  0.000000        0.00000            
    REMARK 350   BIOMT3   1  0.000000  0.000000  1.000000        0.00000            
    REMARK 465                                                                      
    REMARK 465 MISSING RESIDUES                                                     
    REMARK 465 THE FOLLOWING RESIDUES WERE NOT LOCATED IN THE                       
    REMARK 465 EXPERIMENT. (M=MODEL NUMBER; RES=RESIDUE NAME; C=CHAIN               
    REMARK 465 IDENTIFIER; SSSEQ=SEQUENCE NUMBER; I=INSERTION CODE):                
    REMARK 465                                                                      
    REMARK 465   M RES C SSSEQI                                                     
    REMARK 465     MET A    -4                                                      
    REMARK 465     ARG A    -3                                                      
    REMARK 465     THR A    -2                                                      
    REMARK 465     GLU A    -1                                                      
    REMARK 465     ASP A     1                                                      
    REMARK 480                                                                      
    REMARK 480 ZERO OCCUPANCY ATOM                                                  
    REMARK 480 THE FOLLOWING RESIDUES HAVE ATOMS MODELED WITH ZERO                  
    REMARK 480 OCCUPANCY. THE LOCATION AND PROPERTIES OF THESE ATOMS                
    REMARK 480 MAY NOT BE RELIABLE. (M=MODEL NUMBER;                                
    REMARK 480 RES=RESIDUE NAME; C=CHAIN IDENTIFIER; SSEQ=SEQUENCE NUMBER;          
    REMARK 480 I=INSERTION CODE):                                                   
    REMARK 480   M RES CSSEQI  ATOMS                                                
    REMARK 480     LEU A   2    CG  CD1  CD2                                        
    REMARK 500                                                                      
    REMARK 500 GEOMETRY AND STEREOCHEMISTRY                                         
    REMARK 500 SUBTOPIC: TORSION ANGLES                                             
    REMARK 500                                                                      
    REMARK 500 TORSION ANGLES OUTSIDE THE EXPECTED RAMACHANDRAN REGIONS:            
    REMARK 500 (M=MODEL NUMBER; RES=RESIDUE NAME; C=CHAIN IDENTIFIER;               
    REMARK 500 SSEQ=SEQUENCE NUMBER; I=INSERTION CODE).                             
    REMARK 500 STANDARD TABLE:                                                      
    REMARK 500 FORMAT:(10X,I3,1X,A3,1X,A1,I4,A1,4X,F7.2,3X,F7.2)                    
    REMARK 500                                                                      
    REMARK 500 EXPECTED VALUES: GJ KLEYWEGT AND TA JONES (1996). PHI/PSI-           
    REMARK 500 CHOLOGY: RAMACHANDRAN REVISITED. STRUCTURE 4, 1395 - 1400            
    REMARK 500                                                                      
    REMARK 500  M RES CSSEQI        PSI       PHI                                   
    REMARK 500    LYS A  19       -3.45     88.18                                   
    REMARK 500    PRO A  32      -21.16    -37.17                                   
    REMARK 500    ASN A  42      -90.78     64.81                                   
    REMARK 500                                                                      
    REMARK 500 REMARK: NULL                                                         
    REMARK 525                                                                      
    REMARK 525 SOLVENT                                                              
    REMARK 525                                                                      
    REMARK 525 THE SOLVENT MOLECULES HAVE CHAIN IDENTIFIERS THAT                    
    REMARK 525 INDICATE THE POLYMER CHAIN WITH WHICH THEY ARE MOST                  
    REMARK 525 CLOSELY ASSOCIATED. THE REMARK LISTS ALL THE SOLVENT                 
    REMARK 525 MOLECULES WHICH ARE MORE THAN 5A AWAY FROM THE                       
    REMARK 525 NEAREST POLYMER CHAIN (M = MODEL NUMBER;                             
    REMARK 525 RES=RESIDUE NAME; C=CHAIN IDENTIFIER; SSEQ=SEQUENCE                  
    REMARK 525 NUMBER; I=INSERTION CODE):                                           
    REMARK 900                                                                      
    REMARK 900 RELATED ID: 2FCB   RELATED DB: PDB                                   
    REMARK 900  HUMAN FC GAMMA RECEPTOR IIB ECTODOMAIN (CD32)                       
    REMARK 900 RELATED ID: 1E4K   RELATED DB: PDB                                   
    REMARK 900  CRYSTAL STRUCTURE OF SOLUBLE HUMAN IGG1 FC                          
    REMARK 900  FRAGMENT-FC-GAMMA RECEPTOR III COMPLEX                              
    DBREF  1E4J A   -3   172  UNP    O75015   FCG3B_HUMAN     18    193             
    SEQRES   1 A  176  MET ARG THR GLU ASP LEU PRO LYS ALA VAL VAL PHE LEU          
    SEQRES   2 A  176  GLU PRO GLN TRP TYR SER VAL LEU GLU LYS ASP SER VAL          
    SEQRES   3 A  176  THR LEU LYS CYS GLN GLY ALA TYR SER PRO GLU ASP ASN          
    SEQRES   4 A  176  SER THR GLN TRP PHE HIS ASN GLU SER LEU ILE SER SER          
    SEQRES   5 A  176  GLN ALA SER SER TYR PHE ILE ASP ALA ALA THR VAL ASN          
    SEQRES   6 A  176  ASP SER GLY GLU TYR ARG CYS GLN THR ASN LEU SER THR          
    SEQRES   7 A  176  LEU SER ASP PRO VAL GLN LEU GLU VAL HIS ILE GLY TRP          
    SEQRES   8 A  176  LEU LEU LEU GLN ALA PRO ARG TRP VAL PHE LYS GLU GLU          
    SEQRES   9 A  176  ASP PRO ILE HIS LEU ARG CYS HIS SER TRP LYS ASN THR          
    SEQRES  10 A  176  ALA LEU HIS LYS VAL THR TYR LEU GLN ASN GLY LYS ASP          
    SEQRES  11 A  176  ARG LYS TYR PHE HIS HIS ASN SER ASP PHE HIS ILE PRO          
    SEQRES  12 A  176  LYS ALA THR LEU LYS ASP SER GLY SER TYR PHE CYS ARG          
    SEQRES  13 A  176  GLY LEU VAL GLY SER LYS ASN VAL SER SER GLU THR VAL          
    SEQRES  14 A  176  ASN ILE THR ILE THR GLN GLY                                  
    FORMUL   2  HOH   *67(H2 O1)                                                    
    HELIX    1   1 THR A   59  SER A   63  5                                   5    
    HELIX    2   2 THR A  142  SER A  146  5                                   5    
    SHEET    1   A 3 VAL A   6  GLU A  10  0                                        
    SHEET    2   A 3 VAL A  22  GLN A  27 -1  N  GLN A  27   O  VAL A   6           
    SHEET    3   A 3 SER A  52  ILE A  55 -1  N  ILE A  55   O  VAL A  22           
    SHEET    1   B 2 SER A  15  LEU A  17  0                                        
    SHEET    2   B 2 GLU A  82  HIS A  84  1  N  GLU A  82   O  VAL A  16           
    SHEET    1   C 3 GLN A  38  HIS A  41  0                                        
    SHEET    2   C 3 GLY A  64  GLN A  69 -1  N  GLN A  69   O  GLN A  38           
    SHEET    3   C 3 VAL A  79  LEU A  81 -1  N  LEU A  81   O  GLY A  64           
    SHEET    1   D 2 LEU A  88  GLN A  91  0                                        
    SHEET    2   D 2 ARG A 106  SER A 109 -1  N  HIS A 108   O  LEU A  89           
    SHEET    1   E 2 VAL A  96  LYS A  98  0                                        
    SHEET    2   E 2 THR A 168  THR A 170  1  N  THR A 168   O  PHE A  97           
    SHEET    1   F 2 ILE A 103  LEU A 105  0                                        
    SHEET    2   F 2 PHE A 136  ILE A 138 -1  N  ILE A 138   O  ILE A 103           
    SHEET    1   G 4 VAL A 165  ILE A 167  0                                        
    SHEET    2   G 4 GLY A 147  ARG A 152 -1  N  TYR A 149   O  VAL A 165           
    SHEET    3   G 4 VAL A 118  GLN A 122 -1  N  LEU A 121   O  PHE A 150           
    SHEET    4   G 4 LYS A 125  PHE A 130 -1  N  PHE A 130   O  VAL A 118           
    SHEET    1   H 2 ARG A 152  VAL A 155  0                                        
    SHEET    2   H 2 LYS A 158  SER A 161 -1  N  VAL A 160   O  GLY A 153           
    SSBOND   1 CYS A   26    CYS A   68                          1555   1555  2.03  
    SSBOND   2 CYS A  107    CYS A  151                          1555   1555  2.02  
    CISPEP   1 GLU A   10    PRO A   11          0        -0.56                     
    CRYST1   36.700   60.290   85.580  90.00  90.00  90.00 P 2 21 21     4          
    ORIGX1      1.000000  0.000000  0.000000        0.00000                         
    ORIGX2      0.000000  1.000000  0.000000        0.00000                         
    ORIGX3      0.000000  0.000000  1.000000        0.00000                         
    SCALE1      0.027248  0.000000  0.000000        0.00000                         
    SCALE2      0.000000  0.016586  0.000000        0.00000                         
    SCALE3      0.000000  0.000000  0.011685        0.00000                         
    ATOM      1  N   LEU A   2      -7.331  27.189  41.296  1.00 70.01           N  
    ATOM      2  CA  LEU A   2      -5.891  26.836  41.484  1.00 71.04           C  
    ATOM      3  C   LEU A   2      -5.173  26.208  40.263  1.00 69.04           C  
    ATOM      4  O   LEU A   2      -4.000  26.525  40.038  1.00 71.44           O  
    ATOM      5  CB  LEU A   2      -5.749  25.823  42.615  1.00 99.99           C  
    ATOM      6  CG  LEU A   2      -4.399  26.009  43.299  0.00 99.99           C  
    ATOM      7  CD1 LEU A   2      -4.549  25.754  44.795  0.00 99.99           C  
    ATOM      8  CD2 LEU A   2      -3.393  25.025  42.713  0.00 99.99           C  
    ATOM      9  N   PRO A   3      -5.838  25.303  39.479  1.00 64.95           N  
    ATOM     10  CA  PRO A   3      -5.178  24.698  38.307  1.00 58.24           C  
    ATOM     11  C   PRO A   3      -4.991  25.672  37.145  1.00 52.49           C  
    ATOM     12  O   PRO A   3      -5.806  26.572  36.930  1.00 47.15           O  
    ATOM     13  CB  PRO A   3      -6.143  23.578  37.898  1.00 59.87           C  
    ATOM     14  CG  PRO A   3      -6.822  23.234  39.158  1.00 65.35           C  
    ATOM     15  CD  PRO A   3      -7.107  24.593  39.719  1.00 65.11           C  
    ATOM     16  N   LYS A   4      -3.909  25.460  36.404  1.00 50.86           N  
    ATOM     17  CA  LYS A   4      -3.536  26.284  35.261  1.00 43.50           C  
    ATOM     18  C   LYS A   4      -3.863  25.522  33.975  1.00 42.04           C  
    ATOM     19  O   LYS A   4      -3.745  24.292  33.922  1.00 43.05           O  
    ATOM     20  CB  LYS A   4      -2.033  26.588  35.357  1.00 38.98           C  
    ATOM     21  CG  LYS A   4      -1.507  27.723  34.475  1.00 48.83           C  
    ATOM     22  CD  LYS A   4      -0.070  28.137  34.843  1.00 49.21           C  
    ATOM     23  CE  LYS A   4       0.953  27.030  34.575  1.00 56.53           C  
    ATOM     24  NZ  LYS A   4       2.336  27.415  34.984  1.00 58.94           N  
    ATOM     25  N   ALA A   5      -4.301  26.255  32.954  1.00 36.87           N  
    ATOM     26  CA  ALA A   5      -4.638  25.675  31.654  1.00 33.21           C  
    ATOM     27  C   ALA A   5      -3.368  25.371  30.845  1.00 32.44           C  
    ATOM     28  O   ALA A   5      -2.291  25.876  31.173  1.00 36.59           O  
    ATOM     29  CB  ALA A   5      -5.549  26.623  30.884  1.00 20.33           C  
    ATOM     30  N   VAL A   6      -3.470  24.488  29.848  1.00 30.09           N  
    ATOM     31  CA  VAL A   6      -2.320  24.134  28.998  1.00 30.21           C  
    ATOM     32  C   VAL A   6      -2.761  24.210  27.527  1.00 30.97           C  
    ATOM     33  O   VAL A   6      -3.900  23.853  27.201  1.00 25.97           O  
    ATOM     34  CB  VAL A   6      -1.769  22.668  29.275  1.00 31.08           C  
    ATOM     35  CG1 VAL A   6      -0.436  22.434  28.554  1.00 30.12           C  
    ATOM     36  CG2 VAL A   6      -1.568  22.408  30.757  1.00 39.68           C  
    ATOM     37  N   VAL A   7      -1.869  24.708  26.663  1.00 26.05           N  
    ATOM     38  CA  VAL A   7      -2.123  24.790  25.219  1.00 26.58           C  
    ATOM     39  C   VAL A   7      -1.266  23.724  24.519  1.00 27.05           C  
    ATOM     40  O   VAL A   7      -0.043  23.695  24.688  1.00 32.27           O  
    ATOM     41  CB  VAL A   7      -1.784  26.195  24.602  1.00 27.32           C  
    ATOM     42  CG1 VAL A   7      -2.144  26.241  23.113  1.00 24.82           C  
    ATOM     43  CG2 VAL A   7      -2.499  27.303  25.318  1.00 16.61           C  
    ATOM     44  N   PHE A   8      -1.925  22.839  23.770  1.00 29.51           N  
    ATOM     45  CA  PHE A   8      -1.262  21.770  23.009  1.00 28.31           C  
    ATOM     46  C   PHE A   8      -1.348  22.104  21.537  1.00 23.90           C  
    ATOM     47  O   PHE A   8      -2.323  22.699  21.086  1.00 18.64           O  
    ATOM     48  CB  PHE A   8      -1.967  20.411  23.171  1.00 32.79           C  
    ATOM     49  CG  PHE A   8      -2.017  19.895  24.575  1.00 48.33           C  
    ATOM     50  CD1 PHE A   8      -0.831  19.677  25.317  1.00 51.67           C  
    ATOM     51  CD2 PHE A   8      -3.257  19.584  25.159  1.00 50.28           C  
    ATOM     52  CE1 PHE A   8      -0.879  19.147  26.641  1.00 51.01           C  
    ATOM     53  CE2 PHE A   8      -3.330  19.053  26.476  1.00 56.38           C  
    ATOM     54  CZ  PHE A   8      -2.136  18.834  27.223  1.00 56.57           C  
    ATOM     55  N   LEU A   9      -0.341  21.681  20.787  1.00 23.39           N  
    ATOM     56  CA  LEU A   9      -0.320  21.884  19.349  1.00 26.66           C  
    ATOM     57  C   LEU A   9      -0.659  20.548  18.687  1.00 28.57           C  
    ATOM     58  O   LEU A   9      -0.127  19.508  19.087  1.00 27.52           O  
    ATOM     59  CB  LEU A   9       1.061  22.376  18.897  1.00 26.09           C  
    ATOM     60  CG  LEU A   9       1.397  23.865  19.063  1.00 22.19           C  
    ATOM     61  CD1 LEU A   9       2.885  24.091  18.908  1.00 11.75           C  
    ATOM     62  CD2 LEU A   9       0.642  24.695  18.049  1.00 17.70           C  
    ATOM     63  N   GLU A  10      -1.601  20.567  17.744  1.00 30.64           N  
    ATOM     64  CA  GLU A  10      -1.995  19.361  17.011  1.00 30.84           C  
    ATOM     65  C   GLU A  10      -1.923  19.635  15.496  1.00 25.57           C  
    ATOM     66  O   GLU A  10      -2.773  20.350  14.960  1.00 22.29           O  
    ATOM     67  CB  GLU A  10      -3.384  18.867  17.434  1.00 41.19           C  
    ATOM     68  CG  GLU A  10      -3.730  17.479  16.891  1.00 54.09           C  
    ATOM     69  CD  GLU A  10      -4.215  16.531  17.966  1.00 70.24           C  
    ATOM     70  OE1 GLU A  10      -3.419  16.186  18.865  1.00 77.33           O  
    ATOM     71  OE2 GLU A  10      -5.397  16.128  17.912  1.00 77.71           O  
    ATOM     72  N   PRO A  11      -0.901  19.082  14.788  1.00 24.30           N  
    ATOM     73  CA  PRO A  11       0.236  18.227  15.184  1.00 25.76           C  
    ATOM     74  C   PRO A  11       1.300  18.942  16.036  1.00 28.59           C  
    ATOM     75  O   PRO A  11       1.360  20.177  16.047  1.00 19.61           O  
    ATOM     76  CB  PRO A  11       0.817  17.783  13.836  1.00 24.64           C  
    ATOM     77  CG  PRO A  11      -0.306  17.933  12.891  1.00 30.14           C  
    ATOM     78  CD  PRO A  11      -0.900  19.232  13.323  1.00 26.52           C  
    ATOM     79  N   GLN A  12       2.090  18.150  16.766  1.00 29.82           N  
    ATOM     80  CA  GLN A  12       3.140  18.610  17.696  1.00 26.42           C  
    ATOM     81  C   GLN A  12       4.204  19.641  17.279  1.00 23.69           C  
    ATOM     82  O   GLN A  12       4.867  20.208  18.152  1.00 30.21           O  
    ATOM     83  CB  GLN A  12       3.834  17.386  18.329  1.00 33.08           C  
    ATOM     84  CG  GLN A  12       4.427  16.389  17.324  1.00 49.82           C  
    ATOM     85  CD  GLN A  12       4.853  15.069  17.953  1.00 66.38           C  
    ATOM     86  OE1 GLN A  12       6.017  14.673  17.863  1.00 77.43           O  
    ATOM     87  NE2 GLN A  12       3.903  14.365  18.563  1.00 75.67           N  
    ATOM     88  N   TRP A  13       4.304  19.942  15.980  1.00 19.90           N  
    ATOM     89  CA  TRP A  13       5.317  20.874  15.440  1.00 17.55           C  
    ATOM     90  C   TRP A  13       5.126  22.377  15.728  1.00 16.97           C  
    ATOM     91  O   TRP A  13       4.099  22.954  15.366  1.00 15.49           O  
    ATOM     92  CB  TRP A  13       5.495  20.647  13.927  1.00  6.09           C  
    ATOM     93  CG  TRP A  13       5.246  19.211  13.429  1.00 13.19           C  
    ATOM     94  CD1 TRP A  13       4.213  18.797  12.629  1.00  6.55           C  
    ATOM     95  CD2 TRP A  13       6.022  18.029  13.715  1.00  9.77           C  
    ATOM     96  NE1 TRP A  13       4.287  17.438  12.407  1.00  7.60           N  
    ATOM     97  CE2 TRP A  13       5.382  16.939  13.059  1.00  2.92           C  
    ATOM     98  CE3 TRP A  13       7.200  17.780  14.457  1.00 11.70           C  
    ATOM     99  CZ2 TRP A  13       5.874  15.612  13.129  1.00 10.14           C  
    ATOM    100  CZ3 TRP A  13       7.699  16.455  14.525  1.00 13.99           C  
    ATOM    101  CH2 TRP A  13       7.028  15.390  13.859  1.00 14.38           C  
    ATOM    102  N   TYR A  14       6.134  23.005  16.348  1.00 18.64           N  
    ATOM    103  CA  TYR A  14       6.087  24.438  16.687  1.00 25.69           C  
    ATOM    104  C   TYR A  14       6.254  25.335  15.448  1.00 28.52           C  
    ATOM    105  O   TYR A  14       5.712  26.446  15.380  1.00 26.89           O  
    ATOM    106  CB  TYR A  14       7.112  24.803  17.813  1.00 23.57           C  
    ATOM    107  CG  TYR A  14       8.600  24.829  17.442  1.00 18.23           C  
    ATOM    108  CD1 TYR A  14       9.212  26.008  16.953  1.00 22.67           C  
    ATOM    109  CD2 TYR A  14       9.392  23.673  17.546  1.00 24.33           C  
    ATOM    110  CE1 TYR A  14      10.576  26.026  16.562  1.00 28.15           C  
    ATOM    111  CE2 TYR A  14      10.767  23.679  17.161  1.00 25.61           C  
    ATOM    112  CZ  TYR A  14      11.345  24.858  16.669  1.00 29.92           C  
    ATOM    113  OH  TYR A  14      12.667  24.872  16.279  1.00 30.76           O  
    ATOM    114  N   SER A  15       7.008  24.818  14.479  1.00 28.78           N  
    ATOM    115  CA  SER A  15       7.309  25.513  13.238  1.00 21.42           C  
    ATOM    116  C   SER A  15       6.510  24.921  12.084  1.00 22.02           C  
    ATOM    117  O   SER A  15       6.580  23.721  11.815  1.00 22.04           O  
    ATOM    118  CB  SER A  15       8.810  25.414  12.947  1.00 17.25           C  
    ATOM    119  OG  SER A  15       9.201  26.288  11.914  1.00 30.48           O  
    ATOM    120  N   VAL A  16       5.693  25.765  11.459  1.00 20.28           N  
    ATOM    121  CA  VAL A  16       4.883  25.367  10.310  1.00 19.08           C  
    ATOM    122  C   VAL A  16       5.231  26.221   9.093  1.00 19.83           C  
    ATOM    123  O   VAL A  16       5.959  27.213   9.181  1.00 21.11           O  
    ATOM    124  CB  VAL A  16       3.349  25.419  10.588  1.00 17.85           C  
    ATOM    125  CG1 VAL A  16       2.957  24.339  11.581  1.00 25.41           C  
    ATOM    126  CG2 VAL A  16       2.922  26.798  11.097  1.00 27.71           C  
    ATOM    127  N   LEU A  17       4.747  25.778   7.948  1.00 17.80           N  
    ATOM    128  CA  LEU A  17       4.972  26.445   6.682  1.00 17.79           C  
    ATOM    129  C   LEU A  17       3.625  27.064   6.301  1.00 17.82           C  
    ATOM    130  O   LEU A  17       2.592  26.605   6.784  1.00 21.29           O  
    ATOM    131  CB  LEU A  17       5.388  25.368   5.678  1.00 11.59           C  
    ATOM    132  CG  LEU A  17       6.445  25.499   4.583  1.00 17.81           C  
    ATOM    133  CD1 LEU A  17       7.704  26.220   5.027  1.00  4.27           C  
    ATOM    134  CD2 LEU A  17       6.766  24.087   4.120  1.00  7.40           C  
    ATOM    135  N   GLU A  18       3.623  28.131   5.495  1.00 16.40           N  
    ATOM    136  CA  GLU A  18       2.362  28.747   5.054  1.00 12.98           C  
    ATOM    137  C   GLU A  18       1.539  27.727   4.278  1.00 14.44           C  
    ATOM    138  O   GLU A  18       2.100  26.921   3.525  1.00 17.01           O  
    ATOM    139  CB  GLU A  18       2.602  29.972   4.178  1.00 14.61           C  
    ATOM    140  CG  GLU A  18       2.534  31.284   4.916  1.00 13.81           C  
    ATOM    141  CD  GLU A  18       1.914  32.402   4.098  1.00 19.27           C  
    ATOM    142  OE1 GLU A  18       1.426  32.156   2.978  1.00 29.34           O  
    ATOM    143  OE2 GLU A  18       1.895  33.546   4.588  1.00 28.67           O  
    ATOM    144  N   LYS A  19       0.221  27.772   4.490  1.00 11.18           N  
    ATOM    145  CA  LYS A  19      -0.794  26.870   3.907  1.00 14.75           C  
    ATOM    146  C   LYS A  19      -1.036  25.599   4.741  1.00 17.96           C  
    ATOM    147  O   LYS A  19      -1.939  24.810   4.435  1.00 31.88           O  
    ATOM    148  CB  LYS A  19      -0.552  26.526   2.419  1.00 10.34           C  
    ATOM    149  CG  LYS A  19      -0.826  27.676   1.459  1.00 14.93           C  
    ATOM    150  CD  LYS A  19      -0.809  27.211   0.025  1.00 18.49           C  
    ATOM    151  CE  LYS A  19      -1.262  28.321  -0.917  1.00 30.24           C  
    ATOM    152  NZ  LYS A  19      -1.246  27.890  -2.350  1.00 31.94           N  
    ATOM    153  N   ASP A  20      -0.269  25.436   5.820  1.00 11.85           N  
    ATOM    154  CA  ASP A  20      -0.410  24.289   6.718  1.00 14.68           C  
    ATOM    155  C   ASP A  20      -1.560  24.485   7.701  1.00 15.74           C  
    ATOM    156  O   ASP A  20      -1.889  25.614   8.079  1.00 17.96           O  
    ATOM    157  CB  ASP A  20       0.865  24.059   7.545  1.00  9.66           C  
    ATOM    158  CG  ASP A  20       1.977  23.355   6.779  1.00 16.43           C  
    ATOM    159  OD1 ASP A  20       1.764  22.875   5.642  1.00 25.73           O  
    ATOM    160  OD2 ASP A  20       3.089  23.272   7.346  1.00 11.23           O  
    ATOM    161  N   SER A  21      -2.146  23.368   8.119  1.00 16.51           N  
    ATOM    162  CA  SER A  21      -3.231  23.360   9.087  1.00 17.75           C  
    ATOM    163  C   SER A  21      -2.652  23.136  10.475  1.00 18.89           C  
    ATOM    164  O   SER A  21      -1.742  22.326  10.644  1.00 20.02           O  
    ATOM    165  CB  SER A  21      -4.237  22.253   8.770  1.00 17.60           C  
    ATOM    166  OG  SER A  21      -5.031  22.599   7.650  1.00 28.63           O  
    ATOM    167  N   VAL A  22      -3.100  23.945  11.438  1.00 21.28           N  
    ATOM    168  CA  VAL A  22      -2.669  23.824  12.834  1.00 20.52           C  
    ATOM    169  C   VAL A  22      -3.905  23.952  13.729  1.00 21.26           C  
    ATOM    170  O   VAL A  22      -4.784  24.788  13.477  1.00 17.53           O  
    ATOM    171  CB  VAL A  22      -1.668  24.947  13.302  1.00 19.45           C  
    ATOM    172  CG1 VAL A  22      -0.874  24.457  14.499  1.00 16.75           C  
    ATOM    173  CG2 VAL A  22      -0.724  25.390  12.206  1.00 20.89           C  
    ATOM    174  N   THR A  23      -3.969  23.108  14.758  1.00 20.52           N  
    ATOM    175  CA  THR A  23      -5.060  23.148  15.731  1.00 16.24           C  
    ATOM    176  C   THR A  23      -4.429  23.339  17.104  1.00 17.54           C  
    ATOM    177  O   THR A  23      -3.572  22.554  17.525  1.00 11.93           O  
    ATOM    178  CB  THR A  23      -5.933  21.868  15.736  1.00 12.51           C  
    ATOM    179  OG1 THR A  23      -6.366  21.573  14.405  1.00 24.31           O  
    ATOM    180  CG2 THR A  23      -7.169  22.071  16.596  1.00 12.13           C  
    ATOM    181  N   LEU A  24      -4.813  24.439  17.748  1.00 18.09           N  
    ATOM    182  CA  LEU A  24      -4.340  24.795  19.075  1.00 13.49           C  
    ATOM    183  C   LEU A  24      -5.434  24.395  20.063  1.00 20.94           C  
    ATOM    184  O   LEU A  24      -6.548  24.902  19.984  1.00 16.37           O  
    ATOM    185  CB  LEU A  24      -4.054  26.299  19.137  1.00 22.12           C  
    ATOM    186  CG  LEU A  24      -2.814  26.780  18.374  1.00 16.19           C  
    ATOM    187  CD1 LEU A  24      -3.165  27.801  17.327  1.00 11.67           C  
    ATOM    188  CD2 LEU A  24      -1.838  27.352  19.355  1.00 21.85           C  
    ATOM    189  N   LYS A  25      -5.119  23.440  20.941  1.00 19.74           N  
    ATOM    190  CA  LYS A  25      -6.057  22.923  21.945  1.00 22.41           C  
    ATOM    191  C   LYS A  25      -5.827  23.506  23.331  1.00 20.73           C  
    ATOM    192  O   LYS A  25      -4.688  23.635  23.759  1.00 24.77           O  
    ATOM    193  CB  LYS A  25      -5.917  21.407  22.088  1.00 24.24           C  
    ATOM    194  CG  LYS A  25      -6.428  20.552  20.964  1.00 30.01           C  
    ATOM    195  CD  LYS A  25      -6.052  19.109  21.280  1.00 43.19           C  
    ATOM    196  CE  LYS A  25      -6.703  18.112  20.344  1.00 56.11           C  
    ATOM    197  NZ  LYS A  25      -8.174  18.016  20.533  1.00 58.67           N  
    ATOM    198  N   CYS A  26      -6.914  23.819  24.036  1.00 12.33           N  
    ATOM    199  CA  CYS A  26      -6.839  24.324  25.402  1.00  8.42           C  
    ATOM    200  C   CYS A  26      -7.405  23.246  26.321  1.00 19.33           C  
    ATOM    201  O   CYS A  26      -8.507  22.734  26.081  1.00 25.96           O  
    ATOM    202  CB  CYS A  26      -7.660  25.602  25.583  1.00 10.15           C  
    ATOM    203  SG  CYS A  26      -7.342  26.479  27.158  1.00 10.02           S  
    ATOM    204  N   GLN A  27      -6.634  22.884  27.345  1.00 23.99           N  
    ATOM    205  CA  GLN A  27      -7.051  21.888  28.335  1.00 24.25           C  
    ATOM    206  C   GLN A  27      -7.369  22.640  29.638  1.00 24.27           C  
    ATOM    207  O   GLN A  27      -6.462  23.121  30.323  1.00 20.76           O  
    ATOM    208  CB  GLN A  27      -5.936  20.852  28.545  1.00 27.28           C  
    ATOM    209  CG  GLN A  27      -6.290  19.685  29.468  1.00 42.34           C  
    ATOM    210  CD  GLN A  27      -5.150  18.685  29.631  1.00 52.18           C  
    ATOM    211  OE1 GLN A  27      -5.272  17.524  29.236  1.00 56.37           O  
    ATOM    212  NE2 GLN A  27      -4.036  19.134  30.207  1.00 41.26           N  
    ATOM    213  N   GLY A  28      -8.661  22.770  29.940  1.00 23.71           N  
    ATOM    214  CA  GLY A  28      -9.097  23.460  31.145  1.00 27.38           C  
    ATOM    215  C   GLY A  28     -10.391  22.909  31.718  1.00 31.02           C  
    ATOM    216  O   GLY A  28     -11.009  22.020  31.131  1.00 32.51           O  
    ATOM    217  N   ALA A  29     -10.799  23.448  32.867  1.00 34.16           N  
    ATOM    218  CA  ALA A  29     -12.020  23.025  33.552  1.00 37.97           C  
    ATOM    219  C   ALA A  29     -13.154  24.037  33.378  1.00 40.39           C  
    ATOM    220  O   ALA A  29     -12.924  25.249  33.369  1.00 31.64           O  
    ATOM    221  CB  ALA A  29     -11.733  22.792  35.038  1.00 43.62           C  
    ATOM    222  N   TYR A  30     -14.376  23.523  33.218  1.00 45.48           N  
    ATOM    223  CA  TYR A  30     -15.575  24.350  33.027  1.00 45.17           C  
    ATOM    224  C   TYR A  30     -16.693  24.016  34.019  1.00 52.02           C  
    ATOM    225  O   TYR A  30     -16.712  22.936  34.612  1.00 58.26           O  
    ATOM    226  CB  TYR A  30     -16.152  24.154  31.614  1.00 38.74           C  
    ATOM    227  CG  TYR A  30     -15.171  24.276  30.472  1.00 32.83           C  
    ATOM    228  CD1 TYR A  30     -14.404  23.166  30.053  1.00 26.56           C  
    ATOM    229  CD2 TYR A  30     -15.011  25.495  29.790  1.00 35.92           C  
    ATOM    230  CE1 TYR A  30     -13.489  23.270  28.968  1.00 31.38           C  
    ATOM    231  CE2 TYR A  30     -14.105  25.616  28.700  1.00 30.85           C  
    ATOM    232  CZ  TYR A  30     -13.349  24.501  28.298  1.00 32.18           C  
    ATOM    233  OH  TYR A  30     -12.472  24.622  27.243  1.00 28.20           O  
    ATOM    234  N   SER A  31     -17.613  24.965  34.186  1.00 55.71           N  
    ATOM    235  CA  SER A  31     -18.794  24.819  35.037  1.00 58.97           C  
    ATOM    236  C   SER A  31     -19.954  25.138  34.078  1.00 62.00           C  
    ATOM    237  O   SER A  31     -19.845  26.089  33.304  1.00 61.51           O  
    ATOM    238  CB  SER A  31     -18.747  25.796  36.222  1.00 58.02           C  
    ATOM    239  OG  SER A  31     -18.630  27.140  35.801  1.00 50.85           O  
    ATOM    240  N   PRO A  32     -21.051  24.336  34.081  1.00 68.64           N  
    ATOM    241  CA  PRO A  32     -22.241  24.494  33.220  1.00 70.97           C  
    ATOM    242  C   PRO A  32     -22.797  25.881  32.854  1.00 69.70           C  
    ATOM    243  O   PRO A  32     -23.488  26.024  31.839  1.00 68.31           O  
    ATOM    244  CB  PRO A  32     -23.277  23.636  33.932  1.00 73.56           C  
    ATOM    245  CG  PRO A  32     -22.451  22.472  34.350  1.00 75.54           C  
    ATOM    246  CD  PRO A  32     -21.222  23.144  34.943  1.00 72.69           C  
    ATOM    247  N   GLU A  33     -22.462  26.889  33.657  1.00 68.10           N  
    ATOM    248  CA  GLU A  33     -22.904  28.272  33.443  1.00 73.24           C  
    ATOM    249  C   GLU A  33     -22.106  29.037  32.365  1.00 73.56           C  
    ATOM    250  O   GLU A  33     -22.686  29.803  31.585  1.00 76.43           O  
    ATOM    251  CB  GLU A  33     -22.878  29.034  34.776  1.00 73.56           C  
    ATOM    252  CG  GLU A  33     -21.590  28.852  35.579  1.00 81.60           C  
    ATOM    253  CD  GLU A  33     -21.636  29.485  36.950  1.00 89.47           C  
    ATOM    254  OE1 GLU A  33     -22.632  29.273  37.681  1.00 94.77           O  
    ATOM    255  OE2 GLU A  33     -20.656  30.173  37.309  1.00 92.29           O  
    ATOM    256  N   ASP A  34     -20.790  28.808  32.330  1.00 69.94           N  
    ATOM    257  CA  ASP A  34     -19.868  29.447  31.383  1.00 60.85           C  
    ATOM    258  C   ASP A  34     -19.030  28.347  30.703  1.00 56.97           C  
    ATOM    259  O   ASP A  34     -18.096  27.796  31.304  1.00 54.60           O  
    ATOM    260  CB  ASP A  34     -18.969  30.455  32.140  1.00 61.33           C  
    ATOM    261  CG  ASP A  34     -18.194  31.403  31.211  1.00 58.15           C  
    ATOM    262  OD1 ASP A  34     -18.645  31.662  30.074  1.00 57.91           O  
    ATOM    263  OD2 ASP A  34     -17.132  31.909  31.638  1.00 57.44           O  
    ATOM    264  N   ASN A  35     -19.388  28.027  29.457  1.00 49.29           N  
    ATOM    265  CA  ASN A  35     -18.711  26.988  28.667  1.00 48.45           C  
    ATOM    266  C   ASN A  35     -17.645  27.487  27.681  1.00 46.83           C  
    ATOM    267  O   ASN A  35     -16.980  26.678  27.027  1.00 43.43           O  
    ATOM    268  CB  ASN A  35     -19.740  26.145  27.901  1.00 46.15           C  
    ATOM    269  CG  ASN A  35     -20.540  25.229  28.805  1.00 57.11           C  
    ATOM    270  OD1 ASN A  35     -19.988  24.341  29.469  1.00 55.84           O  
    ATOM    271  ND2 ASN A  35     -21.854  25.422  28.818  1.00 63.43           N  
    ATOM    272  N   SER A  36     -17.466  28.807  27.601  1.00 46.44           N  
    ATOM    273  CA  SER A  36     -16.502  29.423  26.684  1.00 40.17           C  
    ATOM    274  C   SER A  36     -15.037  29.376  27.135  1.00 34.79           C  
    ATOM    275  O   SER A  36     -14.742  29.182  28.311  1.00 34.80           O  
    ATOM    276  CB  SER A  36     -16.913  30.871  26.378  1.00 36.98           C  
    ATOM    277  OG  SER A  36     -16.854  31.683  27.539  1.00 45.10           O  
    ATOM    278  N   THR A  37     -14.137  29.526  26.165  1.00 25.55           N  
    ATOM    279  CA  THR A  37     -12.694  29.528  26.376  1.00 16.92           C  
    ATOM    280  C   THR A  37     -12.191  30.866  25.838  1.00 19.33           C  
    ATOM    281  O   THR A  37     -12.686  31.353  24.821  1.00 21.37           O  
    ATOM    282  CB  THR A  37     -12.016  28.385  25.563  1.00 16.42           C  
    ATOM    283  OG1 THR A  37     -12.612  27.135  25.910  1.00 17.63           O  
    ATOM    284  CG2 THR A  37     -10.515  28.306  25.827  1.00  4.24           C  
    ATOM    285  N   GLN A  38     -11.239  31.474  26.538  1.00 16.76           N  
    ATOM    286  CA  GLN A  38     -10.663  32.739  26.088  1.00 20.44           C  
    ATOM    287  C   GLN A  38      -9.299  32.472  25.462  1.00 16.43           C  
    ATOM    288  O   GLN A  38      -8.464  31.804  26.058  1.00 19.68           O  
    ATOM    289  CB  GLN A  38     -10.537  33.744  27.240  1.00 24.28           C  
    ATOM    290  CG  GLN A  38     -11.861  34.274  27.776  1.00 18.67           C  
    ATOM    291  CD  GLN A  38     -12.371  33.506  28.982  1.00 31.09           C  
    ATOM    292  OE1 GLN A  38     -11.723  33.466  30.034  1.00 39.33           O  
    ATOM    293  NE2 GLN A  38     -13.549  32.908  28.844  1.00 32.96           N  
    ATOM    294  N   TRP A  39      -9.110  32.935  24.230  1.00 16.32           N  
    ATOM    295  CA  TRP A  39      -7.847  32.754  23.521  1.00 12.31           C  
    ATOM    296  C   TRP A  39      -7.109  34.071  23.389  1.00 17.81           C  
    ATOM    297  O   TRP A  39      -7.727  35.127  23.212  1.00 17.71           O  
    ATOM    298  CB  TRP A  39      -8.077  32.155  22.134  1.00  5.57           C  
    ATOM    299  CG  TRP A  39      -8.345  30.698  22.131  1.00 13.51           C  
    ATOM    300  CD1 TRP A  39      -9.572  30.086  22.122  1.00  9.46           C  
    ATOM    301  CD2 TRP A  39      -7.371  29.652  22.182  1.00 15.62           C  
    ATOM    302  NE1 TRP A  39      -9.422  28.721  22.176  1.00 23.33           N  
    ATOM    303  CE2 TRP A  39      -8.085  28.421  22.214  1.00 16.82           C  
    ATOM    304  CE3 TRP A  39      -5.957  29.627  22.218  1.00  7.46           C  
    ATOM    305  CZ2 TRP A  39      -7.434  27.171  22.279  1.00 11.30           C  
    ATOM    306  CZ3 TRP A  39      -5.302  28.380  22.288  1.00 12.05           C  
    ATOM    307  CH2 TRP A  39      -6.049  27.168  22.318  1.00 12.78           C  
    ATOM    308  N   PHE A  40      -5.786  34.004  23.507  1.00 12.75           N  
    ATOM    309  CA  PHE A  40      -4.945  35.188  23.396  1.00 12.45           C  
    ATOM    310  C   PHE A  40      -3.837  34.977  22.385  1.00 13.90           C  
    ATOM    311  O   PHE A  40      -3.107  33.992  22.471  1.00 14.32           O  
    ATOM    312  CB  PHE A  40      -4.304  35.556  24.749  1.00  9.14           C  
    ATOM    313  CG  PHE A  40      -5.284  36.018  25.798  1.00 24.04           C  
    ATOM    314  CD1 PHE A  40      -5.894  35.092  26.671  1.00 25.76           C  
    ATOM    315  CD2 PHE A  40      -5.608  37.384  25.923  1.00 31.03           C  
    ATOM    316  CE1 PHE A  40      -6.822  35.518  27.659  1.00 22.62           C  
    ATOM    317  CE2 PHE A  40      -6.534  37.832  26.905  1.00 27.84           C  
    ATOM    318  CZ  PHE A  40      -7.142  36.892  27.776  1.00 25.37           C  
    ATOM    319  N   HIS A  41      -3.776  35.854  21.380  1.00 18.30           N  
    ATOM    320  CA  HIS A  41      -2.706  35.817  20.384  1.00 26.03           C  
    ATOM    321  C   HIS A  41      -1.788  36.955  20.816  1.00 33.54           C  
    ATOM    322  O   HIS A  41      -2.189  38.126  20.829  1.00 33.83           O  
    ATOM    323  CB  HIS A  41      -3.216  36.023  18.952  1.00 21.33           C  
    ATOM    324  CG  HIS A  41      -2.137  35.989  17.906  1.00 21.51           C  
    ATOM    325  ND1 HIS A  41      -2.281  36.594  16.678  1.00 21.07           N  
    ATOM    326  CD2 HIS A  41      -0.890  35.453  17.913  1.00 16.69           C  
    ATOM    327  CE1 HIS A  41      -1.174  36.438  15.978  1.00 12.79           C  
    ATOM    328  NE2 HIS A  41      -0.315  35.749  16.705  1.00 15.50           N  
    ATOM    329  N   ASN A  42      -0.550  36.574  21.142  1.00 36.46           N  
    ATOM    330  CA  ASN A  42       0.506  37.455  21.652  1.00 42.65           C  
    ATOM    331  C   ASN A  42       0.061  37.989  23.026  1.00 49.41           C  
    ATOM    332  O   ASN A  42       0.287  37.322  24.046  1.00 63.95           O  
    ATOM    333  CB  ASN A  42       0.922  38.540  20.632  1.00 36.24           C  
    ATOM    334  CG  ASN A  42       1.762  37.968  19.471  1.00 35.45           C  
    ATOM    335  OD1 ASN A  42       1.544  38.293  18.304  1.00 23.44           O  
    ATOM    336  ND2 ASN A  42       2.730  37.121  19.802  1.00 26.02           N  
    ATOM    337  N   GLU A  43      -0.619  39.132  23.059  1.00 42.61           N  
    ATOM    338  CA  GLU A  43      -1.128  39.694  24.316  1.00 45.24           C  
    ATOM    339  C   GLU A  43      -2.566  40.152  24.105  1.00 39.78           C  
    ATOM    340  O   GLU A  43      -3.245  40.551  25.050  1.00 44.76           O  
    ATOM    341  CB  GLU A  43      -0.275  40.882  24.789  1.00 53.80           C  
    ATOM    342  CG  GLU A  43       1.004  40.511  25.539  1.00 71.65           C  
    ATOM    343  CD  GLU A  43       1.835  41.729  25.916  1.00 84.96           C  
    ATOM    344  OE1 GLU A  43       2.541  42.270  25.033  1.00 87.54           O  
    ATOM    345  OE2 GLU A  43       1.785  42.146  27.095  1.00 90.72           O  
    ATOM    346  N   SER A  44      -3.024  40.055  22.860  1.00 31.78           N  
    ATOM    347  CA  SER A  44      -4.359  40.480  22.465  1.00 28.35           C  
    ATOM    348  C   SER A  44      -5.368  39.336  22.409  1.00 22.63           C  
    ATOM    349  O   SER A  44      -5.062  38.254  21.912  1.00 28.97           O  
    ATOM    350  CB  SER A  44      -4.282  41.183  21.103  1.00 25.69           C  
    ATOM    351  OG  SER A  44      -5.532  41.744  20.728  1.00 41.95           O  
    ATOM    352  N   LEU A  45      -6.578  39.616  22.896  1.00 19.38           N  
    ATOM    353  CA  LEU A  45      -7.701  38.679  22.920  1.00 10.79           C  
    ATOM    354  C   LEU A  45      -8.235  38.501  21.508  1.00 14.56           C  
    ATOM    355  O   LEU A  45      -8.349  39.479  20.769  1.00 26.92           O  
    ATOM    356  CB  LEU A  45      -8.823  39.232  23.818  1.00 15.42           C  
    ATOM    357  CG  LEU A  45     -10.164  38.488  23.957  1.00 19.20           C  
    ATOM    358  CD1 LEU A  45     -10.036  37.302  24.907  1.00 11.13           C  
    ATOM    359  CD2 LEU A  45     -11.239  39.443  24.435  1.00 13.72           C  
    ATOM    360  N   ILE A  46      -8.504  37.252  21.122  1.00 16.78           N  
    ATOM    361  CA  ILE A  46      -9.057  36.962  19.800  1.00 14.48           C  
    ATOM    362  C   ILE A  46     -10.500  36.441  19.864  1.00 18.46           C  
    ATOM    363  O   ILE A  46     -11.025  36.183  20.950  1.00 13.48           O  
    ATOM    364  CB  ILE A  46      -8.125  36.073  18.918  1.00 21.38           C  
    ATOM    365  CG1 ILE A  46      -7.799  34.749  19.599  1.00 15.57           C  
    ATOM    366  CG2 ILE A  46      -6.862  36.851  18.561  1.00  5.50           C  
    ATOM    367  CD1 ILE A  46      -7.129  33.726  18.702  1.00 19.34           C  
    ATOM    368  N   SER A  47     -11.134  36.335  18.694  1.00 24.47           N  
    ATOM    369  CA  SER A  47     -12.540  35.927  18.538  1.00 29.17           C  
    ATOM    370  C   SER A  47     -12.978  34.491  18.812  1.00 30.13           C  
    ATOM    371  O   SER A  47     -14.183  34.226  18.962  1.00 36.10           O  
    ATOM    372  CB  SER A  47     -13.041  36.331  17.148  1.00 25.09           C  
    ATOM    373  OG  SER A  47     -12.939  37.729  16.961  1.00 39.90           O  
    ATOM    374  N   SER A  48     -12.023  33.568  18.868  1.00 26.97           N  
    ATOM    375  CA  SER A  48     -12.340  32.166  19.105  1.00 29.69           C  
    ATOM    376  C   SER A  48     -12.734  31.901  20.550  1.00 26.75           C  
    ATOM    377  O   SER A  48     -12.024  32.292  21.477  1.00 28.33           O  
    ATOM    378  CB  SER A  48     -11.163  31.279  18.710  1.00 41.00           C  
    ATOM    379  OG  SER A  48     -11.508  29.897  18.751  1.00 57.32           O  
    ATOM    380  N   GLN A  49     -13.899  31.283  20.723  1.00 26.00           N  
    ATOM    381  CA  GLN A  49     -14.398  30.941  22.049  1.00 29.39           C  
    ATOM    382  C   GLN A  49     -14.547  29.427  22.274  1.00 28.85           C  
    ATOM    383  O   GLN A  49     -15.065  28.993  23.304  1.00 29.50           O  
    ATOM    384  CB  GLN A  49     -15.697  31.702  22.369  1.00 30.82           C  
    ATOM    385  CG  GLN A  49     -16.824  31.563  21.360  1.00 47.64           C  
    ATOM    386  CD  GLN A  49     -18.085  32.268  21.818  1.00 56.71           C  
    ATOM    387  OE1 GLN A  49     -18.860  31.723  22.606  1.00 63.29           O  
    ATOM    388  NE2 GLN A  49     -18.289  33.493  21.342  1.00 55.45           N  
    ATOM    389  N   ALA A  50     -14.037  28.637  21.326  1.00 24.65           N  
    ATOM    390  CA  ALA A  50     -14.095  27.178  21.391  1.00 21.59           C  
    ATOM    391  C   ALA A  50     -12.900  26.577  22.135  1.00 24.07           C  
    ATOM    392  O   ALA A  50     -11.896  27.258  22.372  1.00 24.90           O  
    ATOM    393  CB  ALA A  50     -14.184  26.604  19.991  1.00 22.51           C  
    ATOM    394  N   SER A  51     -13.022  25.294  22.486  1.00 26.43           N  
    ATOM    395  CA  SER A  51     -11.990  24.520  23.198  1.00 24.69           C  
    ATOM    396  C   SER A  51     -10.685  24.431  22.397  1.00 22.72           C  
    ATOM    397  O   SER A  51      -9.592  24.405  22.961  1.00 27.38           O  
    ATOM    398  CB  SER A  51     -12.513  23.108  23.483  1.00 16.25           C  
    ATOM    399  OG  SER A  51     -11.586  22.351  24.231  1.00 42.12           O  
    ATOM    400  N   SER A  52     -10.823  24.370  21.080  1.00 20.15           N  
    ATOM    401  CA  SER A  52      -9.682  24.308  20.192  1.00 23.68           C  
    ATOM    402  C   SER A  52      -9.797  25.358  19.094  1.00 24.31           C  
    ATOM    403  O   SER A  52     -10.879  25.570  18.538  1.00 28.94           O  
    ATOM    404  CB  SER A  52      -9.530  22.905  19.604  1.00 21.52           C  
    ATOM    405  OG  SER A  52     -10.750  22.434  19.076  1.00 35.30           O  
    ATOM    406  N   TYR A  53      -8.701  26.078  18.867  1.00 16.58           N  
    ATOM    407  CA  TYR A  53      -8.643  27.110  17.834  1.00 18.78           C  
    ATOM    408  C   TYR A  53      -7.986  26.486  16.605  1.00 17.73           C  
    ATOM    409  O   TYR A  53      -6.844  26.014  16.665  1.00 16.56           O  
    ATOM    410  CB  TYR A  53      -7.872  28.349  18.336  1.00 12.58           C  
    ATOM    411  CG  TYR A  53      -7.636  29.453  17.315  1.00 19.09           C  
    ATOM    412  CD1 TYR A  53      -8.680  29.935  16.492  1.00 19.34           C  
    ATOM    413  CD2 TYR A  53      -6.353  30.015  17.154  1.00 20.17           C  
    ATOM    414  CE1 TYR A  53      -8.448  30.954  15.524  1.00 16.81           C  
    ATOM    415  CE2 TYR A  53      -6.111  31.036  16.191  1.00 15.99           C  
    ATOM    416  CZ  TYR A  53      -7.163  31.492  15.385  1.00 13.48           C  
    ATOM    417  OH  TYR A  53      -6.930  32.472  14.448  1.00 24.41           O  
    ATOM    418  N   PHE A  54      -8.729  26.508  15.501  1.00 14.21           N  
    ATOM    419  CA  PHE A  54      -8.289  25.940  14.234  1.00 16.34           C  
    ATOM    420  C   PHE A  54      -7.919  26.990  13.193  1.00 20.37           C  
    ATOM    421  O   PHE A  54      -8.690  27.909  12.898  1.00 26.87           O  
    ATOM    422  CB  PHE A  54      -9.366  24.979  13.690  1.00 17.00           C  
    ATOM    423  CG  PHE A  54      -8.995  24.276  12.403  1.00 10.22           C  
    ATOM    424  CD1 PHE A  54      -7.809  23.511  12.301  1.00 11.17           C  
    ATOM    425  CD2 PHE A  54      -9.842  24.366  11.287  1.00  1.00           C  
    ATOM    426  CE1 PHE A  54      -7.469  22.840  11.092  1.00 13.42           C  
    ATOM    427  CE2 PHE A  54      -9.524  23.700  10.068  1.00 11.66           C  
    ATOM    428  CZ  PHE A  54      -8.334  22.936   9.967  1.00  2.62           C  
    ATOM    429  N   ILE A  55      -6.711  26.812  12.658  1.00 21.02           N  
    ATOM    430  CA  ILE A  55      -6.132  27.646  11.612  1.00 21.75           C  
    ATOM    431  C   ILE A  55      -6.047  26.705  10.394  1.00 28.31           C  
    ATOM    432  O   ILE A  55      -5.193  25.805  10.345  1.00 17.53           O  
    ATOM    433  CB  ILE A  55      -4.713  28.149  12.001  1.00 19.23           C  
    ATOM    434  CG1 ILE A  55      -4.756  28.895  13.340  1.00 12.95           C  
    ATOM    435  CG2 ILE A  55      -4.191  29.070  10.927  1.00 21.43           C  
    ATOM    436  CD1 ILE A  55      -3.409  29.337  13.874  1.00 14.64           C  
    ATOM    437  N   ASP A  56      -6.962  26.909   9.444  1.00 27.93           N  
    ATOM    438  CA  ASP A  56      -7.068  26.103   8.222  1.00 30.44           C  
    ATOM    439  C   ASP A  56      -5.890  26.175   7.251  1.00 28.95           C  
    ATOM    440  O   ASP A  56      -5.453  25.155   6.708  1.00 28.66           O  
    ATOM    441  CB  ASP A  56      -8.400  26.399   7.493  1.00 40.34           C  
    ATOM    442  CG  ASP A  56      -8.590  27.885   7.131  1.00 54.79           C  
    ATOM    443  OD1 ASP A  56      -7.922  28.775   7.713  1.00 55.82           O  
    ATOM    444  OD2 ASP A  56      -9.432  28.156   6.249  1.00 60.21           O  
    ATOM    445  N   ALA A  57      -5.387  27.391   7.059  1.00 21.21           N  
    ATOM    446  CA  ALA A  57      -4.273  27.665   6.168  1.00 15.71           C  
    ATOM    447  C   ALA A  57      -3.488  28.817   6.789  1.00 16.55           C  
    ATOM    448  O   ALA A  57      -3.846  29.997   6.653  1.00 22.89           O  
    ATOM    449  CB  ALA A  57      -4.791  28.028   4.769  1.00  3.19           C  
    ATOM    450  N   ALA A  58      -2.424  28.446   7.493  1.00 14.56           N  
    ATOM    451  CA  ALA A  58      -1.562  29.387   8.199  1.00 13.59           C  
    ATOM    452  C   ALA A  58      -0.868  30.437   7.347  1.00 12.36           C  
    ATOM    453  O   ALA A  58      -0.472  30.172   6.216  1.00 16.98           O  
    ATOM    454  CB  ALA A  58      -0.542  28.622   9.044  1.00  3.79           C  
    ATOM    455  N   THR A  59      -0.858  31.660   7.871  1.00 15.53           N  
    ATOM    456  CA  THR A  59      -0.200  32.802   7.239  1.00 15.75           C  
    ATOM    457  C   THR A  59       0.904  33.255   8.188  1.00 20.30           C  
    ATOM    458  O   THR A  59       1.004  32.728   9.297  1.00 22.78           O  
    ATOM    459  CB  THR A  59      -1.174  33.979   6.927  1.00  8.09           C  
    ATOM    460  OG1 THR A  59      -1.853  34.400   8.119  1.00 15.92           O  
    ATOM    461  CG2 THR A  59      -2.186  33.562   5.875  1.00 11.51           C  
    ATOM    462  N   VAL A  60       1.733  34.209   7.758  1.00 21.48           N  
    ATOM    463  CA  VAL A  60       2.835  34.743   8.577  1.00 23.82           C  
    ATOM    464  C   VAL A  60       2.290  35.535   9.775  1.00 21.77           C  
    ATOM    465  O   VAL A  60       2.951  35.639  10.814  1.00 27.13           O  
    ATOM    466  CB  VAL A  60       3.785  35.660   7.735  1.00 30.69           C  
    ATOM    467  CG1 VAL A  60       5.076  35.984   8.502  1.00 41.51           C  
    ATOM    468  CG2 VAL A  60       4.139  34.994   6.445  1.00 36.48           C  
    ATOM    469  N   ASN A  61       1.065  36.042   9.626  1.00 23.60           N  
    ATOM    470  CA  ASN A  61       0.388  36.824  10.663  1.00 29.53           C  
    ATOM    471  C   ASN A  61      -0.134  35.972  11.823  1.00 24.80           C  
    ATOM    472  O   ASN A  61      -0.579  36.514  12.838  1.00 32.35           O  
    ATOM    473  CB  ASN A  61      -0.752  37.659  10.060  1.00 32.51           C  
    ATOM    474  CG  ASN A  61      -0.271  38.635   8.993  1.00 38.90           C  
    ATOM    475  OD1 ASN A  61      -0.837  38.695   7.900  1.00 46.35           O  
    ATOM    476  ND2 ASN A  61       0.772  39.404   9.306  1.00 35.07           N  
    ATOM    477  N   ASP A  62      -0.074  34.648  11.666  1.00 17.72           N  
    ATOM    478  CA  ASP A  62      -0.509  33.708  12.704  1.00 20.09           C  
    ATOM    479  C   ASP A  62       0.638  33.316  13.641  1.00 18.71           C  
    ATOM    480  O   ASP A  62       0.425  32.583  14.609  1.00 21.36           O  
    ATOM    481  CB  ASP A  62      -1.149  32.448  12.098  1.00 20.64           C  
    ATOM    482  CG  ASP A  62      -2.478  32.726  11.424  1.00 20.71           C  
    ATOM    483  OD1 ASP A  62      -3.360  33.341  12.060  1.00 29.54           O  
    ATOM    484  OD2 ASP A  62      -2.642  32.320  10.254  1.00 23.06           O  
    ATOM    485  N   SER A  63       1.846  33.810  13.354  1.00 13.25           N  
    ATOM    486  CA  SER A  63       3.033  33.543  14.179  1.00 19.38           C  
    ATOM    487  C   SER A  63       2.916  34.290  15.503  1.00 23.29           C  
    ATOM    488  O   SER A  63       2.245  35.329  15.573  1.00 25.83           O  
    ATOM    489  CB  SER A  63       4.318  34.026  13.493  1.00 12.81           C  
    ATOM    490  OG  SER A  63       4.514  33.425  12.233  1.00 30.83           O  
    ATOM    491  N   GLY A  64       3.546  33.747  16.544  1.00 18.13           N  
    ATOM    492  CA  GLY A  64       3.528  34.393  17.844  1.00 16.01           C  
    ATOM    493  C   GLY A  64       3.221  33.474  18.998  1.00 22.54           C  
    ATOM    494  O   GLY A  64       3.258  32.253  18.842  1.00 21.67           O  
    ATOM    495  N   GLU A  65       2.953  34.079  20.159  1.00 22.86           N  
    ATOM    496  CA  GLU A  65       2.615  33.359  21.386  1.00 24.21           C  
    ATOM    497  C   GLU A  65       1.105  33.193  21.512  1.00 20.84           C  
    ATOM    498  O   GLU A  65       0.336  34.085  21.164  1.00 24.07           O  
    ATOM    499  CB  GLU A  65       3.150  34.082  22.629  1.00 29.65           C  
    ATOM    500  CG  GLU A  65       4.672  34.182  22.721  1.00 45.43           C  
    ATOM    501  CD  GLU A  65       5.157  34.774  24.044  1.00 55.37           C  
    ATOM    502  OE1 GLU A  65       4.530  35.734  24.558  1.00 55.50           O  
    ATOM    503  OE2 GLU A  65       6.175  34.270  24.572  1.00 56.96           O  
    ATOM    504  N   TYR A  66       0.690  32.010  21.948  1.00 19.33           N  
    ATOM    505  CA  TYR A  66      -0.719  31.704  22.135  1.00 17.05           C  
    ATOM    506  C   TYR A  66      -0.960  31.276  23.563  1.00 18.93           C  
    ATOM    507  O   TYR A  66      -0.205  30.474  24.112  1.00 19.95           O  
    ATOM    508  CB  TYR A  66      -1.191  30.611  21.164  1.00 15.04           C  
    ATOM    509  CG  TYR A  66      -1.429  31.114  19.754  1.00 24.92           C  
    ATOM    510  CD1 TYR A  66      -0.377  31.157  18.814  1.00 28.61           C  
    ATOM    511  CD2 TYR A  66      -2.696  31.593  19.359  1.00 21.19           C  
    ATOM    512  CE1 TYR A  66      -0.577  31.672  17.509  1.00 25.19           C  
    ATOM    513  CE2 TYR A  66      -2.909  32.112  18.051  1.00 19.46           C  
    ATOM    514  CZ  TYR A  66      -1.845  32.148  17.137  1.00 21.90           C  
    ATOM    515  OH  TYR A  66      -2.032  32.665  15.873  1.00 20.06           O  
    ATOM    516  N   ARG A  67      -1.972  31.882  24.182  1.00 17.07           N  
    ATOM    517  CA  ARG A  67      -2.361  31.575  25.559  1.00 16.73           C  
    ATOM    518  C   ARG A  67      -3.871  31.340  25.625  1.00 18.90           C  
    ATOM    519  O   ARG A  67      -4.619  31.778  24.743  1.00 18.97           O  
    ATOM    520  CB  ARG A  67      -1.969  32.720  26.518  1.00 20.21           C  
    ATOM    521  CG  ARG A  67      -0.470  32.909  26.731  1.00 21.71           C  
    ATOM    522  CD  ARG A  67      -0.147  34.133  27.563  1.00 36.78           C  
    ATOM    523  NE  ARG A  67       0.083  33.800  28.969  1.00 47.46           N  
    ATOM    524  CZ  ARG A  67       1.103  34.249  29.698  1.00 54.32           C  
    ATOM    525  NH1 ARG A  67       2.016  35.060  29.170  1.00 55.21           N  
    ATOM    526  NH2 ARG A  67       1.212  33.886  30.968  1.00 59.62           N  
    ATOM    527  N   CYS A  68      -4.306  30.609  26.650  1.00 15.31           N  
    ATOM    528  CA  CYS A  68      -5.725  30.338  26.855  1.00 15.94           C  
    ATOM    529  C   CYS A  68      -6.071  30.205  28.334  1.00 15.56           C  
    ATOM    530  O   CYS A  68      -5.209  29.898  29.159  1.00 16.49           O  
    ATOM    531  CB  CYS A  68      -6.188  29.088  26.084  1.00  8.46           C  
    ATOM    532  SG  CYS A  68      -5.628  27.493  26.752  1.00 18.47           S  
    ATOM    533  N   GLN A  69      -7.331  30.490  28.650  1.00 11.54           N  
    ATOM    534  CA  GLN A  69      -7.879  30.390  29.998  1.00 15.71           C  
    ATOM    535  C   GLN A  69      -9.377  30.118  29.937  1.00 21.44           C  
    ATOM    536  O   GLN A  69     -10.040  30.394  28.932  1.00 26.59           O  
    ATOM    537  CB  GLN A  69      -7.616  31.656  30.835  1.00 16.86           C  
    ATOM    538  CG  GLN A  69      -8.189  32.959  30.291  1.00 27.47           C  
    ATOM    539  CD  GLN A  69      -8.314  34.047  31.342  1.00 22.40           C  
    ATOM    540  OE1 GLN A  69      -7.433  34.237  32.181  1.00 36.37           O  
    ATOM    541  NE2 GLN A  69      -9.418  34.774  31.293  1.00 32.01           N  
    ATOM    542  N   THR A  70      -9.879  29.491  30.994  1.00 24.52           N  
    ATOM    543  CA  THR A  70     -11.296  29.178  31.142  1.00 24.11           C  
    ATOM    544  C   THR A  70     -11.744  29.883  32.437  1.00 29.32           C  
    ATOM    545  O   THR A  70     -10.956  30.629  33.034  1.00 31.16           O  
    ATOM    546  CB  THR A  70     -11.543  27.642  31.191  1.00 23.93           C  
    ATOM    547  OG1 THR A  70     -10.999  27.092  32.392  1.00 27.70           O  
    ATOM    548  CG2 THR A  70     -10.890  26.949  30.001  1.00 19.03           C  
    ATOM    549  N   ASN A  71     -12.993  29.681  32.862  1.00 33.49           N  
    ATOM    550  CA  ASN A  71     -13.492  30.333  34.081  1.00 36.67           C  
    ATOM    551  C   ASN A  71     -12.901  29.742  35.367  1.00 38.00           C  
    ATOM    552  O   ASN A  71     -12.691  30.466  36.344  1.00 38.18           O  
    ATOM    553  CB  ASN A  71     -15.036  30.361  34.119  1.00 43.15           C  
    ATOM    554  CG  ASN A  71     -15.673  28.965  34.147  1.00 49.91           C  
    ATOM    555  OD1 ASN A  71     -16.342  28.601  35.113  1.00 60.94           O  
    ATOM    556  ND2 ASN A  71     -15.484  28.198  33.082  1.00 52.98           N  
    ATOM    557  N   LEU A  72     -12.568  28.449  35.312  1.00 32.46           N  
    ATOM    558  CA  LEU A  72     -11.997  27.711  36.437  1.00 31.22           C  
    ATOM    559  C   LEU A  72     -10.505  27.419  36.324  1.00 29.17           C  
    ATOM    560  O   LEU A  72      -9.908  26.849  37.244  1.00 30.87           O  
    ATOM    561  CB  LEU A  72     -12.751  26.399  36.641  1.00 34.89           C  
    ATOM    562  CG  LEU A  72     -14.195  26.439  37.135  1.00 38.38           C  
    ATOM    563  CD1 LEU A  72     -14.722  25.019  37.130  1.00 40.76           C  
    ATOM    564  CD2 LEU A  72     -14.294  27.078  38.528  1.00 29.47           C  
    ATOM    565  N   SER A  73      -9.914  27.771  35.185  1.00 26.96           N  
    ATOM    566  CA  SER A  73      -8.485  27.564  34.963  1.00 25.95           C  
    ATOM    567  C   SER A  73      -7.830  28.895  34.647  1.00 24.20           C  
    ATOM    568  O   SER A  73      -8.311  29.640  33.796  1.00 34.40           O  
    ATOM    569  CB  SER A  73      -8.236  26.564  33.829  1.00 23.98           C  
    ATOM    570  OG  SER A  73      -8.773  25.285  34.134  1.00 30.90           O  
    ATOM    571  N   THR A  74      -6.753  29.207  35.364  1.00 27.46           N  
    ATOM    572  CA  THR A  74      -6.021  30.460  35.168  1.00 31.87           C  
    ATOM    573  C   THR A  74      -5.183  30.423  33.891  1.00 31.35           C  
    ATOM    574  O   THR A  74      -4.931  29.344  33.339  1.00 31.04           O  
    ATOM    575  CB  THR A  74      -5.122  30.801  36.379  1.00 32.68           C  
    ATOM    576  OG1 THR A  74      -4.091  29.816  36.515  1.00 40.27           O  
    ATOM    577  CG2 THR A  74      -5.953  30.850  37.659  1.00 40.94           C  
    ATOM    578  N   LEU A  75      -4.755  31.608  33.450  1.00 31.82           N  
    ATOM    579  CA  LEU A  75      -3.963  31.810  32.232  1.00 27.52           C  
    ATOM    580  C   LEU A  75      -2.763  30.868  32.048  1.00 31.31           C  
    ATOM    581  O   LEU A  75      -1.923  30.725  32.944  1.00 38.19           O  
    ATOM    582  CB  LEU A  75      -3.541  33.283  32.133  1.00 22.20           C  
    ATOM    583  CG  LEU A  75      -3.011  33.860  30.815  1.00 17.59           C  
    ATOM    584  CD1 LEU A  75      -4.060  33.771  29.735  1.00  6.99           C  
    ATOM    585  CD2 LEU A  75      -2.602  35.305  31.029  1.00 19.13           C  
    ATOM    586  N   SER A  76      -2.729  30.219  30.879  1.00 28.06           N  
    ATOM    587  CA  SER A  76      -1.694  29.254  30.500  1.00 25.85           C  
    ATOM    588  C   SER A  76      -0.329  29.836  30.163  1.00 28.38           C  
    ATOM    589  O   SER A  76      -0.188  31.037  29.903  1.00 36.62           O  
    ATOM    590  CB  SER A  76      -2.167  28.403  29.310  1.00 15.99           C  
    ATOM    591  OG  SER A  76      -2.035  29.065  28.071  1.00 17.12           O  
    ATOM    592  N   ASP A  77       0.664  28.947  30.128  1.00 29.82           N  
    ATOM    593  CA  ASP A  77       2.036  29.300  29.773  1.00 31.71           C  
    ATOM    594  C   ASP A  77       2.063  29.578  28.262  1.00 31.47           C  
    ATOM    595  O   ASP A  77       1.206  29.066  27.522  1.00 25.48           O  
    ATOM    596  CB  ASP A  77       2.985  28.135  30.093  1.00 33.41           C  
    ATOM    597  CG  ASP A  77       3.304  28.020  31.570  1.00 40.11           C  
    ATOM    598  OD1 ASP A  77       3.539  29.059  32.235  1.00 40.77           O  
    ATOM    599  OD2 ASP A  77       3.345  26.872  32.059  1.00 46.88           O  
    ATOM    600  N   PRO A  78       2.986  30.443  27.791  1.00 36.45           N  
    ATOM    601  CA  PRO A  78       3.036  30.721  26.350  1.00 36.79           C  
    ATOM    602  C   PRO A  78       3.448  29.522  25.491  1.00 35.20           C  
    ATOM    603  O   PRO A  78       4.242  28.673  25.915  1.00 37.07           O  
    ATOM    604  CB  PRO A  78       4.082  31.832  26.258  1.00 40.67           C  
    ATOM    605  CG  PRO A  78       3.937  32.548  27.546  1.00 42.74           C  
    ATOM    606  CD  PRO A  78       3.840  31.404  28.520  1.00 41.72           C  
    ATOM    607  N   VAL A  79       2.808  29.418  24.332  1.00 27.95           N  
    ATOM    608  CA  VAL A  79       3.093  28.384  23.347  1.00 25.91           C  
    ATOM    609  C   VAL A  79       3.433  29.178  22.081  1.00 28.20           C  
    ATOM    610  O   VAL A  79       2.668  30.039  21.647  1.00 16.76           O  
    ATOM    611  CB  VAL A  79       1.900  27.394  23.202  1.00 27.76           C  
    ATOM    612  CG1 VAL A  79       1.817  26.770  21.814  1.00 29.54           C  
    ATOM    613  CG2 VAL A  79       2.071  26.281  24.214  1.00 31.10           C  
    ATOM    614  N   GLN A  80       4.630  28.940  21.558  1.00 29.74           N  
    ATOM    615  CA  GLN A  80       5.102  29.649  20.383  1.00 31.20           C  
    ATOM    616  C   GLN A  80       4.843  28.900  19.074  1.00 25.71           C  
    ATOM    617  O   GLN A  80       5.021  27.682  18.998  1.00 26.25           O  
    ATOM    618  CB  GLN A  80       6.593  29.968  20.551  1.00 36.30           C  
    ATOM    619  CG  GLN A  80       7.127  31.072  19.641  1.00 58.02           C  
    ATOM    620  CD  GLN A  80       8.638  31.246  19.737  1.00 74.21           C  
    ATOM    621  OE1 GLN A  80       9.383  30.282  19.946  1.00 82.37           O  
    ATOM    622  NE2 GLN A  80       9.098  32.479  19.564  1.00 78.56           N  
    ATOM    623  N   LEU A  81       4.362  29.651  18.081  1.00 21.32           N  
    ATOM    624  CA  LEU A  81       4.093  29.153  16.737  1.00 19.28           C  
    ATOM    625  C   LEU A  81       4.870  30.054  15.777  1.00 22.04           C  
    ATOM    626  O   LEU A  81       4.771  31.275  15.867  1.00 22.37           O  
    ATOM    627  CB  LEU A  81       2.594  29.206  16.381  1.00 22.16           C  
    ATOM    628  CG  LEU A  81       2.223  28.642  14.990  1.00 19.32           C  
    ATOM    629  CD1 LEU A  81       2.246  27.127  15.034  1.00 22.43           C  
    ATOM    630  CD2 LEU A  81       0.865  29.124  14.515  1.00 15.19           C  
    ATOM    631  N   GLU A  82       5.651  29.429  14.892  1.00 24.20           N  
    ATOM    632  CA  GLU A  82       6.462  30.105  13.877  1.00 16.62           C  
    ATOM    633  C   GLU A  82       5.963  29.663  12.505  1.00 21.65           C  
    ATOM    634  O   GLU A  82       5.721  28.471  12.292  1.00 20.14           O  
    ATOM    635  CB  GLU A  82       7.928  29.686  13.993  1.00 24.07           C  
    ATOM    636  CG  GLU A  82       8.719  30.270  15.150  1.00 37.41           C  
    ATOM    637  CD  GLU A  82      10.151  29.736  15.210  1.00 41.86           C  
    ATOM    638  OE1 GLU A  82      10.595  29.068  14.246  1.00 37.56           O  
    ATOM    639  OE2 GLU A  82      10.833  29.984  16.232  1.00 49.33           O  
    ATOM    640  N   VAL A  83       5.794  30.618  11.586  1.00 23.88           N  
    ATOM    641  CA  VAL A  83       5.329  30.325  10.220  1.00 17.17           C  
    ATOM    642  C   VAL A  83       6.363  30.774   9.182  1.00 16.06           C  
    ATOM    643  O   VAL A  83       6.689  31.958   9.082  1.00 15.32           O  
    ATOM    644  CB  VAL A  83       3.947  30.982   9.900  1.00 19.86           C  
    ATOM    645  CG1 VAL A  83       3.430  30.503   8.552  1.00 18.78           C  
    ATOM    646  CG2 VAL A  83       2.914  30.662  10.979  1.00 13.32           C  
    ATOM    647  N   HIS A  84       6.840  29.816   8.391  1.00 20.86           N  
    ATOM    648  CA  HIS A  84       7.848  30.069   7.359  1.00 21.82           C  
    ATOM    649  C   HIS A  84       7.316  29.979   5.937  1.00 20.93           C  
    ATOM    650  O   HIS A  84       6.199  29.510   5.697  1.00 13.45           O  
    ATOM    651  CB  HIS A  84       9.010  29.077   7.494  1.00 15.27           C  
    ATOM    652  CG  HIS A  84       9.798  29.232   8.755  1.00 23.78           C  
    ATOM    653  ND1 HIS A  84       9.309  28.864   9.989  1.00 24.60           N  
    ATOM    654  CD2 HIS A  84      11.049  29.700   8.968  1.00 29.89           C  
    ATOM    655  CE1 HIS A  84      10.228  29.097  10.909  1.00 23.84           C  
    ATOM    656  NE2 HIS A  84      11.293  29.604  10.316  1.00 29.73           N  
    ATOM    657  N   ILE A  85       8.116  30.496   5.005  1.00 23.78           N  
    ATOM    658  CA  ILE A  85       7.818  30.452   3.575  1.00 19.26           C  
    ATOM    659  C   ILE A  85       9.079  29.882   2.936  1.00 17.21           C  
    ATOM    660  O   ILE A  85      10.165  30.465   3.021  1.00 17.84           O  
    ATOM    661  CB  ILE A  85       7.472  31.851   2.967  1.00 17.57           C  
    ATOM    662  CG1 ILE A  85       6.182  32.385   3.584  1.00  9.04           C  
    ATOM    663  CG2 ILE A  85       7.248  31.743   1.451  1.00 12.34           C  
    ATOM    664  CD1 ILE A  85       6.057  33.849   3.528  1.00 12.62           C  
    ATOM    665  N   GLY A  86       8.910  28.709   2.345  1.00 14.10           N  
    ATOM    666  CA  GLY A  86       9.994  28.015   1.684  1.00  7.74           C  
    ATOM    667  C   GLY A  86       9.474  26.716   1.107  1.00 13.54           C  
    ATOM    668  O   GLY A  86       8.316  26.338   1.325  1.00 14.69           O  
    ATOM    669  N   TRP A  87      10.327  26.031   0.355  1.00 13.25           N  
    ATOM    670  CA  TRP A  87       9.968  24.753  -0.246  1.00 11.93           C  
    ATOM    671  C   TRP A  87      10.045  23.637   0.786  1.00 19.25           C  
    ATOM    672  O   TRP A  87       9.256  22.692   0.743  1.00 21.90           O  
    ATOM    673  CB  TRP A  87      10.909  24.426  -1.404  1.00 18.26           C  
    ATOM    674  CG  TRP A  87      10.671  25.218  -2.648  1.00 21.24           C  
    ATOM    675  CD1 TRP A  87      11.285  26.386  -3.008  1.00 11.74           C  
    ATOM    676  CD2 TRP A  87       9.775  24.885  -3.715  1.00 13.40           C  
    ATOM    677  NE1 TRP A  87      10.830  26.802  -4.240  1.00 28.69           N  
    ATOM    678  CE2 TRP A  87       9.904  25.903  -4.701  1.00 21.34           C  
    ATOM    679  CE3 TRP A  87       8.874  23.823  -3.940  1.00 16.54           C  
    ATOM    680  CZ2 TRP A  87       9.160  25.895  -5.906  1.00 28.08           C  
    ATOM    681  CZ3 TRP A  87       8.124  23.808  -5.149  1.00 25.27           C  
    ATOM    682  CH2 TRP A  87       8.279  24.844  -6.114  1.00 28.16           C  
    ATOM    683  N   LEU A  88      11.000  23.776   1.708  1.00 21.16           N  
    ATOM    684  CA  LEU A  88      11.258  22.809   2.775  1.00 23.85           C  
    ATOM    685  C   LEU A  88      11.358  23.413   4.174  1.00 23.79           C  
    ATOM    686  O   LEU A  88      11.722  24.579   4.337  1.00 23.10           O  
    ATOM    687  CB  LEU A  88      12.568  22.067   2.502  1.00 26.35           C  
    ATOM    688  CG  LEU A  88      12.607  20.762   1.719  1.00 34.63           C  
    ATOM    689  CD1 LEU A  88      14.031  20.232   1.729  1.00 29.43           C  
    ATOM    690  CD2 LEU A  88      11.673  19.756   2.350  1.00 32.91           C  
    ATOM    691  N   LEU A  89      11.111  22.571   5.176  1.00 17.27           N  
    ATOM    692  CA  LEU A  89      11.178  22.956   6.582  1.00 15.72           C  
    ATOM    693  C   LEU A  89      11.471  21.747   7.477  1.00 13.42           C  
    ATOM    694  O   LEU A  89      10.851  20.697   7.322  1.00 11.47           O  
    ATOM    695  CB  LEU A  89       9.863  23.612   7.020  1.00 11.83           C  
    ATOM    696  CG  LEU A  89       9.735  24.124   8.457  1.00  9.02           C  
    ATOM    697  CD1 LEU A  89      10.584  25.363   8.702  1.00 16.63           C  
    ATOM    698  CD2 LEU A  89       8.299  24.418   8.705  1.00  4.28           C  
    ATOM    699  N   LEU A  90      12.433  21.903   8.392  1.00 13.84           N  
    ATOM    700  CA  LEU A  90      12.779  20.850   9.353  1.00 20.21           C  
    ATOM    701  C   LEU A  90      11.999  21.118  10.638  1.00 21.41           C  
    ATOM    702  O   LEU A  90      12.196  22.139  11.311  1.00 28.32           O  
    ATOM    703  CB  LEU A  90      14.283  20.790   9.648  1.00 11.36           C  
    ATOM    704  CG  LEU A  90      14.675  19.633  10.571  1.00 14.01           C  
    ATOM    705  CD1 LEU A  90      14.830  18.343   9.817  1.00 14.83           C  
    ATOM    706  CD2 LEU A  90      15.931  19.972  11.256  1.00 16.43           C  
    ATOM    707  N   GLN A  91      11.101  20.191  10.948  1.00 16.43           N  
    ATOM    708  CA  GLN A  91      10.246  20.294  12.115  1.00 15.92           C  
    ATOM    709  C   GLN A  91      10.699  19.417  13.278  1.00 20.29           C  
    ATOM    710  O   GLN A  91      11.313  18.365  13.087  1.00 19.80           O  
    ATOM    711  CB  GLN A  91       8.799  19.962  11.732  1.00 17.09           C  
    ATOM    712  CG  GLN A  91       8.189  20.894  10.687  1.00 17.00           C  
    ATOM    713  CD  GLN A  91       6.786  20.491  10.255  1.00 20.12           C  
    ATOM    714  OE1 GLN A  91       6.541  19.345   9.879  1.00 24.48           O  
    ATOM    715  NE2 GLN A  91       5.863  21.443  10.288  1.00 18.18           N  
    ATOM    716  N   ALA A  92      10.381  19.888  14.480  1.00 14.46           N  
    ATOM    717  CA  ALA A  92      10.691  19.228  15.739  1.00 19.34           C  
    ATOM    718  C   ALA A  92       9.617  19.703  16.743  1.00 24.53           C  
    ATOM    719  O   ALA A  92       8.948  20.714  16.479  1.00 22.09           O  
    ATOM    720  CB  ALA A  92      12.091  19.648  16.207  1.00 18.68           C  
    ATOM    721  N   PRO A  93       9.354  18.932  17.836  1.00 22.88           N  
    ATOM    722  CA  PRO A  93       8.351  19.334  18.844  1.00 18.52           C  
    ATOM    723  C   PRO A  93       8.834  20.514  19.727  1.00 24.35           C  
    ATOM    724  O   PRO A  93       8.029  21.253  20.297  1.00 31.82           O  
    ATOM    725  CB  PRO A  93       8.186  18.064  19.685  1.00 18.51           C  
    ATOM    726  CG  PRO A  93       8.515  16.971  18.741  1.00 21.58           C  
    ATOM    727  CD  PRO A  93       9.728  17.518  18.048  1.00 19.45           C  
    ATOM    728  N   ARG A  94      10.157  20.623  19.871  1.00 25.18           N  
    ATOM    729  CA  ARG A  94      10.853  21.665  20.638  1.00 17.70           C  
    ATOM    730  C   ARG A  94      12.305  21.694  20.168  1.00 19.09           C  
    ATOM    731  O   ARG A  94      12.722  20.832  19.392  1.00 18.56           O  
    ATOM    732  CB  ARG A  94      10.777  21.411  22.152  1.00 24.06           C  
    ATOM    733  CG  ARG A  94      11.073  19.984  22.590  1.00 31.78           C  
    ATOM    734  CD  ARG A  94      10.044  19.529  23.610  1.00 39.97           C  
    ATOM    735  NE  ARG A  94      10.105  18.092  23.864  1.00 47.91           N  
    ATOM    736  CZ  ARG A  94      10.931  17.506  24.729  1.00 59.83           C  
    ATOM    737  NH1 ARG A  94      11.789  18.230  25.442  1.00 65.85           N  
    ATOM    738  NH2 ARG A  94      10.906  16.186  24.877  1.00 65.01           N  
    ATOM    739  N   TRP A  95      13.061  22.700  20.600  1.00 21.18           N  
    ATOM    740  CA  TRP A  95      14.463  22.824  20.206  1.00 18.95           C  
    ATOM    741  C   TRP A  95      15.416  22.475  21.343  1.00 16.99           C  
    ATOM    742  O   TRP A  95      16.631  22.485  21.168  1.00 19.71           O  
    ATOM    743  CB  TRP A  95      14.749  24.225  19.622  1.00 22.78           C  
    ATOM    744  CG  TRP A  95      14.420  25.409  20.523  1.00 32.67           C  
    ATOM    745  CD1 TRP A  95      15.296  26.101  21.311  1.00 27.22           C  
    ATOM    746  CD2 TRP A  95      13.137  26.030  20.710  1.00 40.37           C  
    ATOM    747  NE1 TRP A  95      14.646  27.108  21.981  1.00 38.65           N  
    ATOM    748  CE2 TRP A  95      13.322  27.096  21.634  1.00 42.95           C  
    ATOM    749  CE3 TRP A  95      11.845  25.797  20.185  1.00 50.91           C  
    ATOM    750  CZ2 TRP A  95      12.256  27.935  22.057  1.00 50.81           C  
    ATOM    751  CZ3 TRP A  95      10.771  26.634  20.604  1.00 54.75           C  
    ATOM    752  CH2 TRP A  95      10.996  27.693  21.533  1.00 56.44           C  
    ATOM    753  N   VAL A  96      14.852  22.230  22.523  1.00 20.30           N  
    ATOM    754  CA  VAL A  96      15.627  21.862  23.708  1.00 21.07           C  
    ATOM    755  C   VAL A  96      15.092  20.525  24.224  1.00 19.02           C  
    ATOM    756  O   VAL A  96      13.897  20.382  24.476  1.00 21.64           O  
    ATOM    757  CB  VAL A  96      15.544  22.945  24.836  1.00 19.41           C  
    ATOM    758  CG1 VAL A  96      16.364  22.520  26.059  1.00 26.82           C  
    ATOM    759  CG2 VAL A  96      16.075  24.266  24.343  1.00 23.73           C  
    ATOM    760  N   PHE A  97      15.986  19.545  24.330  1.00 15.60           N  
    ATOM    761  CA  PHE A  97      15.649  18.202  24.808  1.00 18.64           C  
    ATOM    762  C   PHE A  97      16.563  17.808  25.969  1.00 21.84           C  
    ATOM    763  O   PHE A  97      17.627  18.397  26.164  1.00 21.19           O  
    ATOM    764  CB  PHE A  97      15.821  17.151  23.693  1.00 10.47           C  
    ATOM    765  CG  PHE A  97      14.905  17.334  22.515  1.00  6.26           C  
    ATOM    766  CD1 PHE A  97      15.263  18.188  21.457  1.00  6.64           C  
    ATOM    767  CD2 PHE A  97      13.706  16.613  22.427  1.00  1.00           C  
    ATOM    768  CE1 PHE A  97      14.438  18.327  20.318  1.00 12.78           C  
    ATOM    769  CE2 PHE A  97      12.866  16.737  21.298  1.00  4.68           C  
    ATOM    770  CZ  PHE A  97      13.233  17.594  20.236  1.00 14.44           C  
    ATOM    771  N   LYS A  98      16.142  16.797  26.726  1.00 25.47           N  
    ATOM    772  CA  LYS A  98      16.913  16.278  27.855  1.00 31.47           C  
    ATOM    773  C   LYS A  98      17.694  15.065  27.316  1.00 26.70           C  
    ATOM    774  O   LYS A  98      17.500  14.680  26.164  1.00 35.81           O  
    ATOM    775  CB  LYS A  98      15.942  15.856  28.960  1.00 31.26           C  
    ATOM    776  CG  LYS A  98      16.518  15.794  30.361  1.00 45.55           C  
    ATOM    777  CD  LYS A  98      15.502  15.227  31.350  1.00 55.75           C  
    ATOM    778  CE  LYS A  98      15.323  13.717  31.185  1.00 59.71           C  
    ATOM    779  NZ  LYS A  98      14.305  13.165  32.117  1.00 61.33           N  
    ATOM    780  N   GLU A  99      18.617  14.508  28.098  1.00 26.58           N  
    ATOM    781  CA  GLU A  99      19.364  13.326  27.655  1.00 26.38           C  
    ATOM    782  C   GLU A  99      18.456  12.114  27.911  1.00 28.28           C  
    ATOM    783  O   GLU A  99      17.692  12.102  28.886  1.00 26.43           O  
    ATOM    784  CB  GLU A  99      20.706  13.199  28.381  1.00 29.13           C  
    ATOM    785  CG  GLU A  99      21.836  12.742  27.455  1.00 32.86           C  
    ATOM    786  CD  GLU A  99      23.216  12.709  28.105  1.00 37.64           C  
    ATOM    787  OE1 GLU A  99      23.457  13.439  29.092  1.00 48.09           O  
    ATOM    788  OE2 GLU A  99      24.075  11.941  27.624  1.00 30.37           O  
    ATOM    789  N   GLU A 100      18.504  11.144  26.990  1.00 28.29           N  
    ATOM    790  CA  GLU A 100      17.675   9.920  26.983  1.00 30.53           C  
    ATOM    791  C   GLU A 100      16.216  10.222  26.567  1.00 30.68           C  
    ATOM    792  O   GLU A 100      15.303   9.410  26.758  1.00 30.13           O  
    ATOM    793  CB  GLU A 100      17.736   9.140  28.312  1.00 39.09           C  
    ATOM    794  CG  GLU A 100      19.038   8.399  28.573  1.00 40.30           C  
    ATOM    795  CD  GLU A 100      19.046   7.734  29.930  1.00 51.47           C  
    ATOM    796  OE1 GLU A 100      18.228   6.811  30.150  1.00 55.42           O  
    ATOM    797  OE2 GLU A 100      19.863   8.144  30.782  1.00 59.50           O  
    ATOM    798  N   ASP A 101      16.035  11.396  25.965  1.00 27.66           N  
    ATOM    799  CA  ASP A 101      14.747  11.880  25.477  1.00 24.84           C  
    ATOM    800  C   ASP A 101      14.591  11.584  23.977  1.00 22.38           C  
    ATOM    801  O   ASP A 101      15.582  11.551  23.246  1.00 23.05           O  
    ATOM    802  CB  ASP A 101      14.656  13.395  25.684  1.00 30.12           C  
    ATOM    803  CG  ASP A 101      13.714  13.794  26.806  1.00 33.44           C  
    ATOM    804  OD1 ASP A 101      13.323  12.940  27.641  1.00 30.90           O  
    ATOM    805  OD2 ASP A 101      13.381  14.998  26.853  1.00 32.23           O  
    ATOM    806  N   PRO A 102      13.349  11.337  23.505  1.00 13.82           N  
    ATOM    807  CA  PRO A 102      13.107  11.055  22.088  1.00 16.79           C  
    ATOM    808  C   PRO A 102      13.083  12.314  21.213  1.00 26.46           C  
    ATOM    809  O   PRO A 102      12.504  13.337  21.601  1.00 28.42           O  
    ATOM    810  CB  PRO A 102      11.721  10.414  22.107  1.00 16.16           C  
    ATOM    811  CG  PRO A 102      11.577   9.915  23.471  1.00 13.33           C  
    ATOM    812  CD  PRO A 102      12.145  11.002  24.282  1.00 15.99           C  
    ATOM    813  N   ILE A 103      13.748  12.243  20.058  1.00 22.21           N  
    ATOM    814  CA  ILE A 103      13.762  13.349  19.101  1.00 21.10           C  
    ATOM    815  C   ILE A 103      13.084  12.869  17.810  1.00 22.42           C  
    ATOM    816  O   ILE A 103      13.534  11.913  17.170  1.00 26.89           O  
    ATOM    817  CB  ILE A 103      15.197  13.884  18.774  1.00 15.03           C  
    ATOM    818  CG1 ILE A 103      15.941  14.301  20.048  1.00 15.30           C  
    ATOM    819  CG2 ILE A 103      15.103  15.114  17.857  1.00 21.29           C  
    ATOM    820  CD1 ILE A 103      17.390  14.698  19.822  1.00 17.66           C  
    ATOM    821  N   HIS A 104      11.976  13.524  17.470  1.00 21.47           N  
    ATOM    822  CA  HIS A 104      11.211  13.226  16.266  1.00 15.24           C  
    ATOM    823  C   HIS A 104      11.368  14.410  15.305  1.00 19.68           C  
    ATOM    824  O   HIS A 104      10.803  15.490  15.524  1.00 20.53           O  
    ATOM    825  CB  HIS A 104       9.727  13.002  16.617  1.00  9.45           C  
    ATOM    826  CG  HIS A 104       8.884  12.462  15.496  1.00  5.26           C  
    ATOM    827  ND1 HIS A 104       7.581  12.060  15.691  1.00 26.18           N  
    ATOM    828  CD2 HIS A 104       9.147  12.244  14.183  1.00 15.07           C  
    ATOM    829  CE1 HIS A 104       7.080  11.614  14.552  1.00 22.95           C  
    ATOM    830  NE2 HIS A 104       8.011  11.716  13.621  1.00 24.07           N  
    ATOM    831  N   LEU A 105      12.142  14.181  14.245  1.00 15.85           N  
    ATOM    832  CA  LEU A 105      12.370  15.181  13.210  1.00 15.47           C  
    ATOM    833  C   LEU A 105      11.596  14.815  11.954  1.00 17.73           C  
    ATOM    834  O   LEU A 105      11.428  13.632  11.641  1.00 20.58           O  
    ATOM    835  CB  LEU A 105      13.858  15.325  12.885  1.00 12.01           C  
    ATOM    836  CG  LEU A 105      14.796  15.788  14.000  1.00 10.28           C  
    ATOM    837  CD1 LEU A 105      16.172  15.960  13.441  1.00 17.98           C  
    ATOM    838  CD2 LEU A 105      14.328  17.077  14.624  1.00  8.35           C  
    ATOM    839  N   ARG A 106      11.117  15.843  11.257  1.00 15.04           N  
    ATOM    840  CA  ARG A 106      10.345  15.683  10.027  1.00 20.19           C  
    ATOM    841  C   ARG A 106      10.779  16.713   8.979  1.00 19.79           C  
    ATOM    842  O   ARG A 106      10.887  17.908   9.269  1.00 15.15           O  
    ATOM    843  CB  ARG A 106       8.843  15.850  10.326  1.00 18.26           C  
    ATOM    844  CG  ARG A 106       7.887  15.604   9.156  1.00 28.70           C  
    ATOM    845  CD  ARG A 106       6.502  16.153   9.464  1.00 37.58           C  
    ATOM    846  NE  ARG A 106       5.443  15.608   8.616  1.00 55.36           N  
    ATOM    847  CZ  ARG A 106       4.646  14.593   8.959  1.00 67.46           C  
    ATOM    848  NH1 ARG A 106       4.774  13.989  10.142  1.00 60.49           N  
    ATOM    849  NH2 ARG A 106       3.708  14.178   8.114  1.00 72.32           N  
    ATOM    850  N   CYS A 107      11.034  16.232   7.765  1.00 16.48           N  
    ATOM    851  CA  CYS A 107      11.407  17.098   6.656  1.00 16.30           C  
    ATOM    852  C   CYS A 107      10.114  17.334   5.874  1.00 15.26           C  
    ATOM    853  O   CYS A 107       9.680  16.502   5.073  1.00 12.24           O  
    ATOM    854  CB  CYS A 107      12.476  16.441   5.788  1.00 23.45           C  
    ATOM    855  SG  CYS A 107      13.267  17.581   4.617  1.00 16.19           S  
    ATOM    856  N   HIS A 108       9.481  18.459   6.190  1.00 15.56           N  
    ATOM    857  CA  HIS A 108       8.207  18.874   5.619  1.00 16.26           C  
    ATOM    858  C   HIS A 108       8.328  19.802   4.408  1.00 18.47           C  
    ATOM    859  O   HIS A 108       8.945  20.868   4.483  1.00 18.28           O  
    ATOM    860  CB  HIS A 108       7.381  19.544   6.732  1.00 11.53           C  
    ATOM    861  CG  HIS A 108       6.004  19.958   6.318  1.00  3.21           C  
    ATOM    862  ND1 HIS A 108       5.104  19.087   5.742  1.00 16.89           N  
    ATOM    863  CD2 HIS A 108       5.378  21.155   6.388  1.00  7.02           C  
    ATOM    864  CE1 HIS A 108       3.983  19.732   5.474  1.00 14.30           C  
    ATOM    865  NE2 HIS A 108       4.124  20.988   5.856  1.00 13.63           N  
    ATOM    866  N   SER A 109       7.688  19.391   3.314  1.00 16.24           N  
    ATOM    867  CA  SER A 109       7.657  20.153   2.072  1.00 14.59           C  
    ATOM    868  C   SER A 109       6.331  20.921   1.966  1.00 19.87           C  
    ATOM    869  O   SER A 109       5.368  20.598   2.665  1.00 23.54           O  
    ATOM    870  CB  SER A 109       7.857  19.224   0.860  1.00 17.28           C  
    ATOM    871  OG  SER A 109       6.896  18.177   0.820  1.00 26.64           O  
    ATOM    872  N   TRP A 110       6.289  21.918   1.078  1.00 17.01           N  
    ATOM    873  CA  TRP A 110       5.116  22.767   0.830  1.00 17.76           C  
    ATOM    874  C   TRP A 110       3.857  21.953   0.502  1.00 18.62           C  
    ATOM    875  O   TRP A 110       3.966  20.818   0.041  1.00 27.07           O  
    ATOM    876  CB  TRP A 110       5.456  23.721  -0.321  1.00 19.59           C  
    ATOM    877  CG  TRP A 110       4.479  24.846  -0.568  1.00 21.19           C  
    ATOM    878  CD1 TRP A 110       3.503  24.886  -1.529  1.00 24.30           C  
    ATOM    879  CD2 TRP A 110       4.431  26.111   0.102  1.00 18.95           C  
    ATOM    880  NE1 TRP A 110       2.855  26.099  -1.504  1.00 20.80           N  
    ATOM    881  CE2 TRP A 110       3.397  26.874  -0.516  1.00 11.30           C  
    ATOM    882  CE3 TRP A 110       5.161  26.687   1.163  1.00 11.59           C  
    ATOM    883  CZ2 TRP A 110       3.071  28.186  -0.109  1.00  8.79           C  
    ATOM    884  CZ3 TRP A 110       4.839  28.006   1.577  1.00 14.61           C  
    ATOM    885  CH2 TRP A 110       3.798  28.735   0.932  1.00 13.71           C  
    ATOM    886  N   LYS A 111       2.678  22.526   0.762  1.00 21.05           N  
    ATOM    887  CA  LYS A 111       1.387  21.870   0.493  1.00 28.31           C  
    ATOM    888  C   LYS A 111       1.256  21.468  -0.986  1.00 33.63           C  
    ATOM    889  O   LYS A 111       1.582  22.257  -1.878  1.00 43.06           O  
    ATOM    890  CB  LYS A 111       0.219  22.789   0.902  1.00 29.13           C  
    ATOM    891  CG  LYS A 111      -1.171  22.129   0.869  1.00 44.50           C  
    ATOM    892  CD  LYS A 111      -2.285  23.093   1.255  1.00 55.39           C  
    ATOM    893  CE  LYS A 111      -3.674  22.485   1.052  1.00 62.52           C  
    ATOM    894  NZ  LYS A 111      -3.947  21.326   1.948  1.00 63.59           N  
    ATOM    895  N   ASN A 112       0.828  20.220  -1.210  1.00 36.91           N  
    ATOM    896  CA  ASN A 112       0.626  19.594  -2.536  1.00 43.05           C  
    ATOM    897  C   ASN A 112       1.897  19.239  -3.318  1.00 41.79           C  
    ATOM    898  O   ASN A 112       1.824  18.868  -4.495  1.00 45.27           O  
    ATOM    899  CB  ASN A 112      -0.360  20.386  -3.425  1.00 51.74           C  
    ATOM    900  CG  ASN A 112      -1.781  20.375  -2.883  1.00 63.91           C  
    ATOM    901  OD1 ASN A 112      -2.384  19.313  -2.696  1.00 71.95           O  
    ATOM    902  ND2 ASN A 112      -2.323  21.562  -2.625  1.00 67.33           N  
    ATOM    903  N   THR A 113       3.055  19.351  -2.661  1.00 33.23           N  
    ATOM    904  CA  THR A 113       4.329  19.010  -3.286  1.00 27.73           C  
    ATOM    905  C   THR A 113       4.879  17.724  -2.677  1.00 33.61           C  
    ATOM    906  O   THR A 113       4.638  17.426  -1.504  1.00 33.61           O  
    ATOM    907  CB  THR A 113       5.386  20.132  -3.160  1.00 32.15           C  
    ATOM    908  OG1 THR A 113       5.655  20.394  -1.780  1.00 32.56           O  
    ATOM    909  CG2 THR A 113       4.923  21.416  -3.859  1.00 32.01           C  
    ATOM    910  N   ALA A 114       5.585  16.956  -3.503  1.00 39.10           N  
    ATOM    911  CA  ALA A 114       6.184  15.692  -3.092  1.00 37.75           C  
    ATOM    912  C   ALA A 114       7.702  15.814  -3.018  1.00 37.03           C  
    ATOM    913  O   ALA A 114       8.331  16.499  -3.832  1.00 30.31           O  
    ATOM    914  CB  ALA A 114       5.785  14.578  -4.057  1.00 44.66           C  
    ATOM    915  N   LEU A 115       8.276  15.110  -2.051  1.00 36.65           N  
    ATOM    916  CA  LEU A 115       9.709  15.121  -1.819  1.00 36.59           C  
    ATOM    917  C   LEU A 115      10.324  13.727  -1.924  1.00 37.70           C  
    ATOM    918  O   LEU A 115       9.704  12.729  -1.540  1.00 41.78           O  
    ATOM    919  CB  LEU A 115       9.976  15.703  -0.431  1.00 36.44           C  
    ATOM    920  CG  LEU A 115      11.395  16.023   0.033  1.00 34.55           C  
    ATOM    921  CD1 LEU A 115      11.933  17.257  -0.655  1.00 38.02           C  
    ATOM    922  CD2 LEU A 115      11.367  16.222   1.530  1.00 41.04           C  
    ATOM    923  N   HIS A 116      11.524  13.674  -2.499  1.00 35.59           N  
    ATOM    924  CA  HIS A 116      12.294  12.436  -2.636  1.00 40.64           C  
    ATOM    925  C   HIS A 116      13.788  12.737  -2.481  1.00 32.60           C  
    ATOM    926  O   HIS A 116      14.201  13.898  -2.551  1.00 27.98           O  
    ATOM    927  CB  HIS A 116      11.974  11.662  -3.943  1.00 48.56           C  
    ATOM    928  CG  HIS A 116      12.162  12.451  -5.202  1.00 62.88           C  
    ATOM    929  ND1 HIS A 116      11.265  13.411  -5.617  1.00 71.01           N  
    ATOM    930  CD2 HIS A 116      13.141  12.419  -6.138  1.00 67.26           C  
    ATOM    931  CE1 HIS A 116      11.685  13.939  -6.754  1.00 75.02           C  
    ATOM    932  NE2 HIS A 116      12.821  13.355  -7.091  1.00 74.65           N  
    ATOM    933  N   LYS A 117      14.577  11.686  -2.239  1.00 33.24           N  
    ATOM    934  CA  LYS A 117      16.034  11.756  -2.036  1.00 33.93           C  
    ATOM    935  C   LYS A 117      16.411  12.698  -0.881  1.00 32.23           C  
    ATOM    936  O   LYS A 117      17.262  13.591  -1.005  1.00 34.13           O  
    ATOM    937  CB  LYS A 117      16.776  12.074  -3.348  1.00 32.26           C  
    ATOM    938  CG  LYS A 117      16.548  11.026  -4.432  1.00 30.85           C  
    ATOM    939  CD  LYS A 117      17.669  11.013  -5.432  1.00 48.66           C  
    ATOM    940  CE  LYS A 117      17.566   9.791  -6.318  1.00 56.52           C  
    ATOM    941  NZ  LYS A 117      18.745   9.658  -7.225  1.00 65.67           N  
    ATOM    942  N   VAL A 118      15.731  12.461   0.243  1.00 27.46           N  
    ATOM    943  CA  VAL A 118      15.870  13.226   1.486  1.00 24.10           C  
    ATOM    944  C   VAL A 118      17.093  12.849   2.309  1.00 23.01           C  
    ATOM    945  O   VAL A 118      17.325  11.672   2.583  1.00 17.14           O  
    ATOM    946  CB  VAL A 118      14.632  13.046   2.409  1.00 21.01           C  
    ATOM    947  CG1 VAL A 118      14.596  14.125   3.486  1.00 19.95           C  
    ATOM    948  CG2 VAL A 118      13.371  13.068   1.605  1.00 19.95           C  
    ATOM    949  N   THR A 119      17.809  13.872   2.769  1.00 22.07           N  
    ATOM    950  CA  THR A 119      18.996  13.694   3.595  1.00 24.86           C  
    ATOM    951  C   THR A 119      18.875  14.584   4.832  1.00 29.16           C  
    ATOM    952  O   THR A 119      18.601  15.780   4.716  1.00 31.76           O  
    ATOM    953  CB  THR A 119      20.285  14.072   2.824  1.00 20.61           C  
    ATOM    954  OG1 THR A 119      20.287  13.427   1.546  1.00 33.76           O  
    ATOM    955  CG2 THR A 119      21.513  13.629   3.589  1.00 23.56           C  
    ATOM    956  N   TYR A 120      19.036  13.980   6.008  1.00 21.78           N  
    ATOM    957  CA  TYR A 120      18.993  14.712   7.274  1.00 19.91           C  
    ATOM    958  C   TYR A 120      20.440  14.908   7.676  1.00 23.63           C  
    ATOM    959  O   TYR A 120      21.223  13.950   7.700  1.00 20.82           O  
    ATOM    960  CB  TYR A 120      18.260  13.921   8.357  1.00 11.25           C  
    ATOM    961  CG  TYR A 120      16.768  13.872   8.186  1.00  7.66           C  
    ATOM    962  CD1 TYR A 120      16.162  12.870   7.397  1.00  3.11           C  
    ATOM    963  CD2 TYR A 120      15.938  14.802   8.842  1.00  9.71           C  
    ATOM    964  CE1 TYR A 120      14.744  12.791   7.269  1.00 17.80           C  
    ATOM    965  CE2 TYR A 120      14.518  14.734   8.724  1.00 10.14           C  
    ATOM    966  CZ  TYR A 120      13.932  13.725   7.939  1.00 11.55           C  
    ATOM    967  OH  TYR A 120      12.561  13.640   7.842  1.00 15.46           O  
    ATOM    968  N   LEU A 121      20.810  16.160   7.929  1.00 23.40           N  
    ATOM    969  CA  LEU A 121      22.185  16.477   8.298  1.00 22.03           C  
    ATOM    970  C   LEU A 121      22.332  17.131   9.663  1.00 22.93           C  
    ATOM    971  O   LEU A 121      21.439  17.848  10.115  1.00 21.54           O  
    ATOM    972  CB  LEU A 121      22.848  17.372   7.233  1.00 21.01           C  
    ATOM    973  CG  LEU A 121      22.936  16.948   5.757  1.00 22.14           C  
    ATOM    974  CD1 LEU A 121      21.712  17.448   4.989  1.00 10.98           C  
    ATOM    975  CD2 LEU A 121      24.209  17.496   5.133  1.00 16.53           C  
    ATOM    976  N   GLN A 122      23.442  16.815  10.331  1.00 20.46           N  
    ATOM    977  CA  GLN A 122      23.794  17.387  11.630  1.00 19.79           C  
    ATOM    978  C   GLN A 122      25.226  17.890  11.482  1.00 24.32           C  
    ATOM    979  O   GLN A 122      26.129  17.127  11.104  1.00 25.30           O  
    ATOM    980  CB  GLN A 122      23.726  16.356  12.761  1.00 19.12           C  
    ATOM    981  CG  GLN A 122      24.056  16.931  14.147  1.00 17.73           C  
    ATOM    982  CD  GLN A 122      24.375  15.883  15.204  1.00 20.17           C  
    ATOM    983  OE1 GLN A 122      24.155  14.683  15.017  1.00 27.94           O  
    ATOM    984  NE2 GLN A 122      24.892  16.345  16.335  1.00 13.07           N  
    ATOM    985  N   ASN A 123      25.418  19.174  11.783  1.00 16.55           N  
    ATOM    986  CA  ASN A 123      26.713  19.855  11.694  1.00 24.17           C  
    ATOM    987  C   ASN A 123      27.420  19.663  10.344  1.00 24.79           C  
    ATOM    988  O   ASN A 123      28.631  19.400  10.277  1.00 26.19           O  
    ATOM    989  CB  ASN A 123      27.622  19.473  12.870  1.00 19.08           C  
    ATOM    990  CG  ASN A 123      27.040  19.871  14.213  1.00 26.36           C  
    ATOM    991  OD1 ASN A 123      26.268  20.828  14.321  1.00 24.77           O  
    ATOM    992  ND2 ASN A 123      27.409  19.130  15.249  1.00 36.40           N  
    ATOM    993  N   GLY A 124      26.627  19.793   9.278  1.00 22.97           N  
    ATOM    994  CA  GLY A 124      27.110  19.654   7.912  1.00 29.78           C  
    ATOM    995  C   GLY A 124      27.312  18.233   7.409  1.00 32.20           C  
    ATOM    996  O   GLY A 124      27.552  18.034   6.213  1.00 37.15           O  
    ATOM    997  N   LYS A 125      27.198  17.257   8.314  1.00 31.67           N  
    ATOM    998  CA  LYS A 125      27.394  15.842   7.994  1.00 39.52           C  
    ATOM    999  C   LYS A 125      26.102  15.031   7.963  1.00 37.88           C  
    ATOM   1000  O   LYS A 125      25.216  15.229   8.792  1.00 33.27           O  
    ATOM   1001  CB  LYS A 125      28.404  15.218   8.974  1.00 39.51           C  
    ATOM   1002  CG  LYS A 125      29.821  15.775   8.803  1.00 59.27           C  
    ATOM   1003  CD  LYS A 125      30.785  15.376   9.914  1.00 71.67           C  
    ATOM   1004  CE  LYS A 125      32.143  16.061   9.702  1.00 78.59           C  
    ATOM   1005  NZ  LYS A 125      33.154  15.758  10.758  1.00 79.61           N  
    ATOM   1006  N   ASP A 126      26.031  14.098   7.010  1.00 40.54           N  
    ATOM   1007  CA  ASP A 126      24.880  13.206   6.804  1.00 40.09           C  
    ATOM   1008  C   ASP A 126      24.615  12.282   7.995  1.00 36.12           C  
    ATOM   1009  O   ASP A 126      25.538  11.638   8.505  1.00 41.34           O  
    ATOM   1010  CB  ASP A 126      25.099  12.324   5.560  1.00 40.89           C  
    ATOM   1011  CG  ASP A 126      25.232  13.118   4.264  1.00 48.05           C  
    ATOM   1012  OD1 ASP A 126      25.620  14.308   4.294  1.00 59.50           O  
    ATOM   1013  OD2 ASP A 126      24.960  12.529   3.194  1.00 50.16           O  
    ATOM   1014  N   ARG A 127      23.363  12.253   8.447  1.00 28.14           N  
    ATOM   1015  CA  ARG A 127      22.953  11.390   9.556  1.00 28.77           C  
    ATOM   1016  C   ARG A 127      22.171  10.215   9.011  1.00 26.02           C  
    ATOM   1017  O   ARG A 127      22.431   9.064   9.374  1.00 26.49           O  
    ATOM   1018  CB  ARG A 127      22.090  12.150  10.561  1.00 29.81           C  
    ATOM   1019  CG  ARG A 127      22.876  13.003  11.503  1.00 30.67           C  
    ATOM   1020  CD  ARG A 127      23.539  12.187  12.595  1.00 38.86           C  
    ATOM   1021  NE  ARG A 127      22.818  12.298  13.861  1.00 33.58           N  
    ATOM   1022  CZ  ARG A 127      22.283  11.278  14.519  1.00 29.15           C  
    ATOM   1023  NH1 ARG A 127      22.371  10.038  14.044  1.00 40.96           N  
    ATOM   1024  NH2 ARG A 127      21.651  11.503  15.658  1.00 24.68           N  
    ATOM   1025  N   LYS A 128      21.202  10.523   8.150  1.00 25.23           N  
    ATOM   1026  CA  LYS A 128      20.368   9.508   7.527  1.00 27.77           C  
    ATOM   1027  C   LYS A 128      19.828   9.940   6.178  1.00 26.25           C  
    ATOM   1028  O   LYS A 128      19.438  11.093   5.982  1.00 19.99           O  
    ATOM   1029  CB  LYS A 128      19.207   9.098   8.441  1.00 28.16           C  
    ATOM   1030  CG  LYS A 128      18.870   7.616   8.343  1.00 35.94           C  
    ATOM   1031  CD  LYS A 128      18.159   7.134   9.587  1.00 52.96           C  
    ATOM   1032  CE  LYS A 128      18.049   5.625   9.635  1.00 59.94           C  
    ATOM   1033  NZ  LYS A 128      17.349   5.195  10.874  1.00 69.06           N  
    ATOM   1034  N   TYR A 129      19.862   8.988   5.250  1.00 27.76           N  
    ATOM   1035  CA  TYR A 129      19.363   9.159   3.898  1.00 23.85           C  
    ATOM   1036  C   TYR A 129      18.121   8.276   3.718  1.00 25.74           C  
    ATOM   1037  O   TYR A 129      18.018   7.178   4.277  1.00 37.11           O  
    ATOM   1038  CB  TYR A 129      20.447   8.798   2.861  1.00 20.52           C  
    ATOM   1039  CG  TYR A 129      19.981   8.867   1.413  1.00 27.32           C  
    ATOM   1040  CD1 TYR A 129      19.791  10.111   0.768  1.00 28.79           C  
    ATOM   1041  CD2 TYR A 129      19.664   7.687   0.696  1.00 31.74           C  
    ATOM   1042  CE1 TYR A 129      19.285  10.181  -0.566  1.00 26.24           C  
    ATOM   1043  CE2 TYR A 129      19.156   7.741  -0.628  1.00 28.33           C  
    ATOM   1044  CZ  TYR A 129      18.969   8.991  -1.251  1.00 35.11           C  
    ATOM   1045  OH  TYR A 129      18.463   9.045  -2.531  1.00 33.33           O  
    ATOM   1046  N   PHE A 130      17.179   8.788   2.938  1.00 25.76           N  
    ATOM   1047  CA  PHE A 130      15.940   8.098   2.615  1.00 27.56           C  
    ATOM   1048  C   PHE A 130      15.650   8.391   1.147  1.00 31.50           C  
    ATOM   1049  O   PHE A 130      15.813   9.527   0.689  1.00 27.12           O  
    ATOM   1050  CB  PHE A 130      14.778   8.612   3.477  1.00 30.96           C  
    ATOM   1051  CG  PHE A 130      14.799   8.132   4.906  1.00 38.31           C  
    ATOM   1052  CD1 PHE A 130      14.640   6.762   5.211  1.00 33.42           C  
    ATOM   1053  CD2 PHE A 130      14.950   9.058   5.961  1.00 31.73           C  
    ATOM   1054  CE1 PHE A 130      14.623   6.310   6.558  1.00 33.24           C  
    ATOM   1055  CE2 PHE A 130      14.937   8.630   7.315  1.00 43.33           C  
    ATOM   1056  CZ  PHE A 130      14.772   7.247   7.615  1.00 44.71           C  
    ATOM   1057  N   HIS A 131      15.252   7.358   0.406  1.00 34.72           N  
    ATOM   1058  CA  HIS A 131      14.920   7.500  -1.009  1.00 30.99           C  
    ATOM   1059  C   HIS A 131      13.596   8.236  -1.150  1.00 29.89           C  
    ATOM   1060  O   HIS A 131      13.451   9.106  -2.009  1.00 25.57           O  
    ATOM   1061  CB  HIS A 131      14.842   6.129  -1.686  1.00 33.40           C  
    ATOM   1062  CG  HIS A 131      16.178   5.552  -2.044  1.00 34.93           C  
    ATOM   1063  ND1 HIS A 131      16.725   4.474  -1.382  1.00 44.90           N  
    ATOM   1064  CD2 HIS A 131      17.068   5.893  -3.008  1.00 34.23           C  
    ATOM   1065  CE1 HIS A 131      17.893   4.174  -1.923  1.00 42.11           C  
    ATOM   1066  NE2 HIS A 131      18.124   5.020  -2.912  1.00 40.37           N  
    ATOM   1067  N   HIS A 132      12.654   7.891  -0.269  1.00 32.03           N  
    ATOM   1068  CA  HIS A 132      11.316   8.493  -0.223  1.00 38.73           C  
    ATOM   1069  C   HIS A 132      11.194   9.221   1.116  1.00 40.08           C  
    ATOM   1070  O   HIS A 132      11.854   8.836   2.082  1.00 49.33           O  
    ATOM   1071  CB  HIS A 132      10.237   7.408  -0.324  1.00 42.22           C  
    ATOM   1072  CG  HIS A 132      10.461   6.435  -1.442  1.00 58.96           C  
    ATOM   1073  ND1 HIS A 132      10.522   6.822  -2.765  1.00 55.09           N  
    ATOM   1074  CD2 HIS A 132      10.668   5.095  -1.432  1.00 62.90           C  
    ATOM   1075  CE1 HIS A 132      10.759   5.763  -3.520  1.00 60.69           C  
    ATOM   1076  NE2 HIS A 132      10.851   4.702  -2.736  1.00 62.57           N  
    ATOM   1077  N   ASN A 133      10.320  10.227   1.187  1.00 39.17           N  
    ATOM   1078  CA  ASN A 133      10.126  11.030   2.403  1.00 34.78           C  
    ATOM   1079  C   ASN A 133       9.590  10.317   3.655  1.00 36.33           C  
    ATOM   1080  O   ASN A 133       8.399  10.005   3.757  1.00 45.52           O  
    ATOM   1081  CB  ASN A 133       9.278  12.275   2.091  1.00 41.53           C  
    ATOM   1082  CG  ASN A 133       9.396  13.367   3.164  1.00 47.53           C  
    ATOM   1083  OD1 ASN A 133      10.459  13.566   3.765  1.00 46.20           O  
    ATOM   1084  ND2 ASN A 133       8.300  14.085   3.392  1.00 45.13           N  
    ATOM   1085  N   SER A 134      10.502  10.076   4.596  1.00 31.05           N  
    ATOM   1086  CA  SER A 134      10.205   9.445   5.881  1.00 30.45           C  
    ATOM   1087  C   SER A 134      10.739  10.364   6.971  1.00 26.25           C  
    ATOM   1088  O   SER A 134      11.572  11.236   6.698  1.00 28.01           O  
    ATOM   1089  CB  SER A 134      10.876   8.075   5.986  1.00 28.88           C  
    ATOM   1090  OG  SER A 134      10.294   7.152   5.083  1.00 48.96           O  
    ATOM   1091  N   ASP A 135      10.240  10.192   8.194  1.00 23.86           N  
    ATOM   1092  CA  ASP A 135      10.675  11.001   9.338  1.00 21.44           C  
    ATOM   1093  C   ASP A 135      11.994  10.454   9.893  1.00 23.38           C  
    ATOM   1094  O   ASP A 135      12.411   9.351   9.525  1.00 25.28           O  
    ATOM   1095  CB  ASP A 135       9.606  11.001  10.443  1.00 21.45           C  
    ATOM   1096  CG  ASP A 135       8.310  11.692  10.031  1.00 15.08           C  
    ATOM   1097  OD1 ASP A 135       8.252  12.304   8.946  1.00 26.06           O  
    ATOM   1098  OD2 ASP A 135       7.335  11.623  10.806  1.00 19.00           O  
    ATOM   1099  N   PHE A 136      12.665  11.245  10.731  1.00 21.52           N  
    ATOM   1100  CA  PHE A 136      13.930  10.837  11.346  1.00 17.82           C  
    ATOM   1101  C   PHE A 136      13.707  10.741  12.846  1.00 21.45           C  
    ATOM   1102  O   PHE A 136      13.433  11.735  13.524  1.00 16.75           O  
    ATOM   1103  CB  PHE A 136      15.049  11.835  11.000  1.00 15.79           C  
    ATOM   1104  CG  PHE A 136      16.428  11.480  11.539  1.00 16.25           C  
    ATOM   1105  CD1 PHE A 136      16.893  10.145  11.594  1.00 19.62           C  
    ATOM   1106  CD2 PHE A 136      17.271  12.505  12.002  1.00 17.46           C  
    ATOM   1107  CE1 PHE A 136      18.181   9.840  12.113  1.00 22.54           C  
    ATOM   1108  CE2 PHE A 136      18.557  12.219  12.518  1.00 23.57           C  
    ATOM   1109  CZ  PHE A 136      19.013  10.877  12.574  1.00 24.19           C  
    ATOM   1110  N   HIS A 137      13.865   9.520  13.344  1.00 22.39           N  
    ATOM   1111  CA  HIS A 137      13.675   9.219  14.751  1.00 22.06           C  
    ATOM   1112  C   HIS A 137      14.977   8.892  15.476  1.00 19.66           C  
    ATOM   1113  O   HIS A 137      15.816   8.145  14.968  1.00 26.75           O  
    ATOM   1114  CB  HIS A 137      12.715   8.031  14.903  1.00 12.39           C  
    ATOM   1115  CG  HIS A 137      11.453   8.157  14.105  1.00 10.85           C  
    ATOM   1116  ND1 HIS A 137      11.200   7.392  12.986  1.00 18.00           N  
    ATOM   1117  CD2 HIS A 137      10.375   8.959  14.261  1.00 17.59           C  
    ATOM   1118  CE1 HIS A 137      10.022   7.720  12.487  1.00 20.13           C  
    ATOM   1119  NE2 HIS A 137       9.500   8.670  13.242  1.00 17.41           N  
    ATOM   1120  N   ILE A 138      15.168   9.529  16.629  1.00 18.29           N  
    ATOM   1121  CA  ILE A 138      16.311   9.280  17.515  1.00 15.72           C  
    ATOM   1122  C   ILE A 138      15.573   8.952  18.827  1.00 22.77           C  
    ATOM   1123  O   ILE A 138      15.176   9.867  19.553  1.00 24.90           O  
    ATOM   1124  CB  ILE A 138      17.221  10.533  17.730  1.00 17.94           C  
    ATOM   1125  CG1 ILE A 138      17.744  11.078  16.394  1.00 20.93           C  
    ATOM   1126  CG2 ILE A 138      18.396  10.172  18.645  1.00 11.77           C  
    ATOM   1127  CD1 ILE A 138      18.382  12.463  16.472  1.00 16.12           C  
    ATOM   1128  N   PRO A 139      15.332   7.648  19.119  1.00 25.87           N  
    ATOM   1129  CA  PRO A 139      14.627   7.219  20.336  1.00 21.71           C  
    ATOM   1130  C   PRO A 139      15.247   7.699  21.642  1.00 20.92           C  
    ATOM   1131  O   PRO A 139      14.530   8.090  22.552  1.00 26.60           O  
    ATOM   1132  CB  PRO A 139      14.670   5.693  20.234  1.00 24.13           C  
    ATOM   1133  CG  PRO A 139      14.695   5.454  18.770  1.00 29.14           C  
    ATOM   1134  CD  PRO A 139      15.725   6.460  18.337  1.00 26.33           C  
    ATOM   1135  N   LYS A 140      16.577   7.660  21.720  1.00 27.72           N  
    ATOM   1136  CA  LYS A 140      17.313   8.095  22.904  1.00 26.47           C  
    ATOM   1137  C   LYS A 140      18.399   9.085  22.528  1.00 25.61           C  
    ATOM   1138  O   LYS A 140      19.394   8.727  21.882  1.00 31.60           O  
    ATOM   1139  CB  LYS A 140      17.924   6.902  23.656  1.00 33.22           C  
    ATOM   1140  CG  LYS A 140      17.131   6.475  24.886  1.00 37.40           C  
    ATOM   1141  CD  LYS A 140      17.893   5.482  25.748  1.00 41.03           C  
    ATOM   1142  CE  LYS A 140      17.063   5.077  26.963  1.00 47.34           C  
    ATOM   1143  NZ  LYS A 140      17.793   4.131  27.858  1.00 51.22           N  
    ATOM   1144  N   ALA A 141      18.198  10.333  22.942  1.00 24.64           N  
    ATOM   1145  CA  ALA A 141      19.139  11.416  22.669  1.00 27.12           C  
    ATOM   1146  C   ALA A 141      20.373  11.366  23.550  1.00 23.93           C  
    ATOM   1147  O   ALA A 141      20.292  11.050  24.740  1.00 22.40           O  
    ATOM   1148  CB  ALA A 141      18.457  12.765  22.813  1.00 17.28           C  
    ATOM   1149  N   THR A 142      21.515  11.622  22.919  1.00 29.30           N  
    ATOM   1150  CA  THR A 142      22.821  11.663  23.576  1.00 33.46           C  
    ATOM   1151  C   THR A 142      23.318  13.117  23.474  1.00 34.08           C  
    ATOM   1152  O   THR A 142      22.654  13.963  22.864  1.00 33.62           O  
    ATOM   1153  CB  THR A 142      23.844  10.699  22.888  1.00 30.14           C  
    ATOM   1154  OG1 THR A 142      24.006  11.059  21.510  1.00 41.62           O  
    ATOM   1155  CG2 THR A 142      23.372   9.259  22.963  1.00 28.62           C  
    ATOM   1156  N   LEU A 143      24.476  13.400  24.071  1.00 31.41           N  
    ATOM   1157  CA  LEU A 143      25.081  14.735  24.034  1.00 32.92           C  
    ATOM   1158  C   LEU A 143      25.732  14.976  22.662  1.00 32.97           C  
    ATOM   1159  O   LEU A 143      26.061  16.110  22.306  1.00 40.45           O  
    ATOM   1160  CB  LEU A 143      26.103  14.883  25.172  1.00 37.27           C  
    ATOM   1161  CG  LEU A 143      25.968  16.001  26.215  1.00 26.67           C  
    ATOM   1162  CD1 LEU A 143      24.579  16.056  26.842  1.00 35.35           C  
    ATOM   1163  CD2 LEU A 143      27.011  15.777  27.290  1.00 30.88           C  
    ATOM   1164  N   LYS A 144      25.878  13.893  21.893  1.00 29.40           N  
    ATOM   1165  CA  LYS A 144      26.433  13.919  20.542  1.00 29.37           C  
    ATOM   1166  C   LYS A 144      25.378  14.429  19.554  1.00 25.39           C  
    ATOM   1167  O   LYS A 144      25.709  14.833  18.439  1.00 28.03           O  
    ATOM   1168  CB  LYS A 144      26.866  12.516  20.110  1.00 37.40           C  
    ATOM   1169  CG  LYS A 144      28.211  12.042  20.631  1.00 53.22           C  
    ATOM   1170  CD  LYS A 144      28.561  10.699  19.986  1.00 67.01           C  
    ATOM   1171  CE  LYS A 144      29.970  10.235  20.319  1.00 71.72           C  
    ATOM   1172  NZ  LYS A 144      30.294   8.961  19.612  1.00 75.23           N  
    ATOM   1173  N   ASP A 145      24.117  14.409  19.987  1.00 22.57           N  
    ATOM   1174  CA  ASP A 145      22.981  14.851  19.184  1.00 22.01           C  
    ATOM   1175  C   ASP A 145      22.716  16.363  19.159  1.00 22.69           C  
    ATOM   1176  O   ASP A 145      21.785  16.814  18.495  1.00 22.53           O  
    ATOM   1177  CB  ASP A 145      21.713  14.101  19.604  1.00 19.03           C  
    ATOM   1178  CG  ASP A 145      21.774  12.629  19.276  1.00 22.83           C  
    ATOM   1179  OD1 ASP A 145      22.079  12.275  18.116  1.00 24.72           O  
    ATOM   1180  OD2 ASP A 145      21.523  11.818  20.184  1.00 29.27           O  
    ATOM   1181  N   SER A 146      23.517  17.140  19.884  1.00 17.92           N  
    ATOM   1182  CA  SER A 146      23.362  18.593  19.904  1.00 24.29           C  
    ATOM   1183  C   SER A 146      24.010  19.227  18.672  1.00 26.66           C  
    ATOM   1184  O   SER A 146      25.041  18.746  18.198  1.00 29.92           O  
    ATOM   1185  CB  SER A 146      23.954  19.181  21.186  1.00 28.07           C  
    ATOM   1186  OG  SER A 146      25.307  18.801  21.357  1.00 46.78           O  
    ATOM   1187  N   GLY A 147      23.362  20.249  18.113  1.00 22.56           N  
    ATOM   1188  CA  GLY A 147      23.918  20.913  16.947  1.00 18.96           C  
    ATOM   1189  C   GLY A 147      22.944  21.474  15.937  1.00 23.09           C  
    ATOM   1190  O   GLY A 147      21.735  21.523  16.178  1.00 28.37           O  
    ATOM   1191  N   SER A 148      23.495  21.886  14.795  1.00 23.29           N  
    ATOM   1192  CA  SER A 148      22.737  22.464  13.693  1.00 22.65           C  
    ATOM   1193  C   SER A 148      22.224  21.407  12.726  1.00 25.44           C  
    ATOM   1194  O   SER A 148      22.999  20.692  12.088  1.00 25.64           O  
    ATOM   1195  CB  SER A 148      23.591  23.476  12.936  1.00 26.04           C  
    ATOM   1196  OG  SER A 148      23.876  24.605  13.740  1.00 49.38           O  
    ATOM   1197  N   TYR A 149      20.904  21.318  12.632  1.00 16.84           N  
    ATOM   1198  CA  TYR A 149      20.261  20.365  11.752  1.00 14.87           C  
    ATOM   1199  C   TYR A 149      19.502  21.057  10.626  1.00 17.16           C  
    ATOM   1200  O   TYR A 149      19.021  22.183  10.777  1.00 20.50           O  
    ATOM   1201  CB  TYR A 149      19.255  19.510  12.525  1.00  8.77           C  
    ATOM   1202  CG  TYR A 149      19.797  18.497  13.502  1.00 20.95           C  
    ATOM   1203  CD1 TYR A 149      20.266  18.889  14.775  1.00 17.37           C  
    ATOM   1204  CD2 TYR A 149      19.738  17.117  13.206  1.00 20.15           C  
    ATOM   1205  CE1 TYR A 149      20.656  17.928  15.738  1.00 16.91           C  
    ATOM   1206  CE2 TYR A 149      20.121  16.146  14.160  1.00 16.66           C  
    ATOM   1207  CZ  TYR A 149      20.576  16.563  15.420  1.00 22.25           C  
    ATOM   1208  OH  TYR A 149      20.941  15.622  16.344  1.00 27.18           O  
    ATOM   1209  N   PHE A 150      19.419  20.363   9.496  1.00 14.96           N  
    ATOM   1210  CA  PHE A 150      18.671  20.791   8.317  1.00 12.89           C  
    ATOM   1211  C   PHE A 150      18.495  19.567   7.418  1.00 17.67           C  
    ATOM   1212  O   PHE A 150      19.192  18.553   7.588  1.00 19.26           O  
    ATOM   1213  CB  PHE A 150      19.295  22.020   7.585  1.00 10.95           C  
    ATOM   1214  CG  PHE A 150      20.522  21.726   6.747  1.00 14.77           C  
    ATOM   1215  CD1 PHE A 150      21.788  21.577   7.346  1.00 20.29           C  
    ATOM   1216  CD2 PHE A 150      20.418  21.626   5.340  1.00 13.35           C  
    ATOM   1217  CE1 PHE A 150      22.950  21.327   6.556  1.00 17.88           C  
    ATOM   1218  CE2 PHE A 150      21.560  21.378   4.536  1.00 18.59           C  
    ATOM   1219  CZ  PHE A 150      22.834  21.227   5.148  1.00 22.14           C  
    ATOM   1220  N   CYS A 151      17.514  19.634   6.524  1.00 18.04           N  
    ATOM   1221  CA  CYS A 151      17.264  18.544   5.597  1.00 19.59           C  
    ATOM   1222  C   CYS A 151      17.286  19.057   4.166  1.00 18.37           C  
    ATOM   1223  O   CYS A 151      16.932  20.204   3.898  1.00 20.61           O  
    ATOM   1224  CB  CYS A 151      15.975  17.787   5.943  1.00 12.61           C  
    ATOM   1225  SG  CYS A 151      14.430  18.737   5.799  1.00 19.27           S  
    ATOM   1226  N   ARG A 152      17.763  18.203   3.269  1.00 26.65           N  
    ATOM   1227  CA  ARG A 152      17.916  18.520   1.855  1.00 27.43           C  
    ATOM   1228  C   ARG A 152      17.242  17.433   1.034  1.00 24.63           C  
    ATOM   1229  O   ARG A 152      17.286  16.258   1.393  1.00 31.31           O  
    ATOM   1230  CB  ARG A 152      19.424  18.596   1.547  1.00 30.05           C  
    ATOM   1231  CG  ARG A 152      19.848  18.910   0.116  1.00 44.25           C  
    ATOM   1232  CD  ARG A 152      21.356  19.153   0.037  1.00 44.30           C  
    ATOM   1233  NE  ARG A 152      22.146  18.112   0.708  1.00 52.52           N  
    ATOM   1234  CZ  ARG A 152      22.643  17.020   0.122  1.00 61.11           C  
    ATOM   1235  NH1 ARG A 152      22.449  16.788  -1.175  1.00 63.24           N  
    ATOM   1236  NH2 ARG A 152      23.343  16.152   0.841  1.00 56.68           N  
    ATOM   1237  N   GLY A 153      16.623  17.837  -0.068  1.00 26.23           N  
    ATOM   1238  CA  GLY A 153      15.946  16.893  -0.935  1.00 28.10           C  
    ATOM   1239  C   GLY A 153      15.589  17.467  -2.287  1.00 27.65           C  
    ATOM   1240  O   GLY A 153      15.963  18.594  -2.623  1.00 24.21           O  
    ATOM   1241  N   LEU A 154      14.877  16.660  -3.069  1.00 26.12           N  
    ATOM   1242  CA  LEU A 154      14.442  17.026  -4.405  1.00 22.90           C  
    ATOM   1243  C   LEU A 154      12.925  17.152  -4.480  1.00 27.71           C  
    ATOM   1244  O   LEU A 154      12.214  16.231  -4.084  1.00 24.27           O  
    ATOM   1245  CB  LEU A 154      14.899  15.968  -5.424  1.00 29.93           C  
    ATOM   1246  CG  LEU A 154      16.311  15.835  -6.017  1.00 32.10           C  
    ATOM   1247  CD1 LEU A 154      16.635  17.048  -6.860  1.00 40.67           C  
    ATOM   1248  CD2 LEU A 154      17.370  15.611  -4.954  1.00 38.89           C  
    ATOM   1249  N   VAL A 155      12.442  18.313  -4.935  1.00 31.49           N  
    ATOM   1250  CA  VAL A 155      11.006  18.559  -5.140  1.00 37.64           C  
    ATOM   1251  C   VAL A 155      10.871  18.519  -6.668  1.00 44.40           C  
    ATOM   1252  O   VAL A 155      10.928  19.556  -7.344  1.00 48.20           O  
    ATOM   1253  CB  VAL A 155      10.534  19.938  -4.573  1.00 34.93           C  
    ATOM   1254  CG1 VAL A 155       9.056  20.167  -4.885  1.00 33.93           C  
    ATOM   1255  CG2 VAL A 155      10.715  19.979  -3.070  1.00 43.63           C  
    ATOM   1256  N   GLY A 156      10.749  17.298  -7.195  1.00 47.37           N  
    ATOM   1257  CA  GLY A 156      10.667  17.082  -8.631  1.00 43.59           C  
    ATOM   1258  C   GLY A 156      12.093  17.114  -9.156  1.00 46.27           C  
    ATOM   1259  O   GLY A 156      12.832  16.136  -9.038  1.00 45.71           O  
    ATOM   1260  N   SER A 157      12.495  18.288  -9.639  1.00 47.77           N  
    ATOM   1261  CA  SER A 157      13.837  18.526 -10.170  1.00 51.39           C  
    ATOM   1262  C   SER A 157      14.646  19.466  -9.278  1.00 51.86           C  
    ATOM   1263  O   SER A 157      15.879  19.393  -9.257  1.00 55.88           O  
    ATOM   1264  CB  SER A 157      13.749  19.129 -11.571  1.00 54.44           C  
    ATOM   1265  OG  SER A 157      12.749  20.133 -11.650  1.00 56.19           O  
    ATOM   1266  N   LYS A 158      13.938  20.338  -8.551  1.00 47.37           N  
    ATOM   1267  CA  LYS A 158      14.535  21.333  -7.651  1.00 40.81           C  
    ATOM   1268  C   LYS A 158      15.272  20.781  -6.440  1.00 33.03           C  
    ATOM   1269  O   LYS A 158      14.689  20.080  -5.612  1.00 23.90           O  
    ATOM   1270  CB  LYS A 158      13.478  22.334  -7.155  1.00 38.57           C  
    ATOM   1271  CG  LYS A 158      13.195  23.490  -8.094  1.00 45.74           C  
    ATOM   1272  CD  LYS A 158      12.981  24.797  -7.330  1.00 48.33           C  
    ATOM   1273  CE  LYS A 158      14.308  25.459  -6.940  1.00 50.67           C  
    ATOM   1274  NZ  LYS A 158      14.125  26.752  -6.220  1.00 48.06           N  
    ATOM   1275  N   ASN A 159      16.556  21.119  -6.349  1.00 30.20           N  
    ATOM   1276  CA  ASN A 159      17.392  20.709  -5.231  1.00 36.53           C  
    ATOM   1277  C   ASN A 159      17.283  21.833  -4.202  1.00 38.15           C  
    ATOM   1278  O   ASN A 159      17.831  22.926  -4.389  1.00 38.63           O  
    ATOM   1279  CB  ASN A 159      18.845  20.505  -5.673  1.00 38.34           C  
    ATOM   1280  CG  ASN A 159      19.501  19.323  -4.983  1.00 48.93           C  
    ATOM   1281  OD1 ASN A 159      19.542  19.241  -3.748  1.00 52.21           O  
    ATOM   1282  ND2 ASN A 159      20.000  18.386  -5.779  1.00 54.19           N  
    ATOM   1283  N   VAL A 160      16.478  21.581  -3.171  1.00 36.96           N  
    ATOM   1284  CA  VAL A 160      16.227  22.549  -2.105  1.00 31.73           C  
    ATOM   1285  C   VAL A 160      16.606  22.083  -0.701  1.00 27.67           C  
    ATOM   1286  O   VAL A 160      16.458  20.912  -0.355  1.00 29.90           O  
    ATOM   1287  CB  VAL A 160      14.743  23.038  -2.099  1.00 32.38           C  
    ATOM   1288  CG1 VAL A 160      14.485  23.973  -3.275  1.00 32.38           C  
    ATOM   1289  CG2 VAL A 160      13.777  21.862  -2.148  1.00 19.08           C  
    ATOM   1290  N   SER A 161      17.092  23.032   0.092  1.00 32.75           N  
    ATOM   1291  CA  SER A 161      17.504  22.808   1.474  1.00 28.77           C  
    ATOM   1292  C   SER A 161      16.627  23.641   2.406  1.00 26.58           C  
    ATOM   1293  O   SER A 161      16.078  24.674   2.004  1.00 34.31           O  
    ATOM   1294  CB  SER A 161      18.971  23.208   1.672  1.00 18.13           C  
    ATOM   1295  OG  SER A 161      19.839  22.422   0.876  1.00 41.73           O  
    ATOM   1296  N   SER A 162      16.505  23.187   3.651  1.00 22.04           N  
    ATOM   1297  CA  SER A 162      15.718  23.886   4.659  1.00 14.72           C  
    ATOM   1298  C   SER A 162      16.622  24.819   5.456  1.00 17.74           C  
    ATOM   1299  O   SER A 162      17.854  24.747   5.350  1.00 16.74           O  
    ATOM   1300  CB  SER A 162      15.078  22.876   5.615  1.00 19.02           C  
    ATOM   1301  OG  SER A 162      16.034  22.294   6.486  1.00 15.64           O  
    ATOM   1302  N   GLU A 163      16.007  25.679   6.267  1.00 14.08           N  
    ATOM   1303  CA  GLU A 163      16.760  26.573   7.135  1.00 18.10           C  
    ATOM   1304  C   GLU A 163      17.254  25.728   8.315  1.00 23.90           C  
    ATOM   1305  O   GLU A 163      16.797  24.596   8.521  1.00 19.46           O  
    ATOM   1306  CB  GLU A 163      15.892  27.733   7.622  1.00 26.70           C  
    ATOM   1307  CG  GLU A 163      15.594  28.769   6.556  1.00 37.32           C  
    ATOM   1308  CD  GLU A 163      14.557  29.769   7.007  1.00 49.56           C  
    ATOM   1309  OE1 GLU A 163      13.348  29.493   6.819  1.00 51.03           O  
    ATOM   1310  OE2 GLU A 163      14.949  30.825   7.551  1.00 51.00           O  
    ATOM   1311  N   THR A 164      18.217  26.267   9.052  1.00 29.78           N  
    ATOM   1312  CA  THR A 164      18.815  25.586  10.193  1.00 26.28           C  
    ATOM   1313  C   THR A 164      17.992  25.664  11.490  1.00 22.50           C  
    ATOM   1314  O   THR A 164      17.355  26.684  11.778  1.00 27.91           O  
    ATOM   1315  CB  THR A 164      20.274  26.111  10.382  1.00 27.93           C  
    ATOM   1316  OG1 THR A 164      21.093  25.574   9.340  1.00 38.79           O  
    ATOM   1317  CG2 THR A 164      20.880  25.729  11.722  1.00 47.13           C  
    ATOM   1318  N   VAL A 165      17.954  24.539  12.208  1.00 19.29           N  
    ATOM   1319  CA  VAL A 165      17.284  24.419  13.510  1.00 22.47           C  
    ATOM   1320  C   VAL A 165      18.410  24.005  14.461  1.00 24.18           C  
    ATOM   1321  O   VAL A 165      19.133  23.045  14.187  1.00 29.88           O  
    ATOM   1322  CB  VAL A 165      16.136  23.338  13.528  1.00 16.45           C  
    ATOM   1323  CG1 VAL A 165      15.506  23.237  14.914  1.00 14.57           C  
    ATOM   1324  CG2 VAL A 165      15.033  23.698  12.537  1.00 18.52           C  
    ATOM   1325  N   ASN A 166      18.585  24.768  15.542  1.00 27.04           N  
    ATOM   1326  CA  ASN A 166      19.623  24.491  16.531  1.00 22.14           C  
    ATOM   1327  C   ASN A 166      19.050  23.742  17.723  1.00 22.89           C  
    ATOM   1328  O   ASN A 166      18.235  24.273  18.482  1.00 28.38           O  
    ATOM   1329  CB  ASN A 166      20.306  25.783  16.978  1.00 26.03           C  
    ATOM   1330  CG  ASN A 166      20.937  26.539  15.824  1.00 38.06           C  
    ATOM   1331  OD1 ASN A 166      20.556  27.671  15.528  1.00 54.02           O  
    ATOM   1332  ND2 ASN A 166      21.893  25.909  15.154  1.00 31.61           N  
    ATOM   1333  N   ILE A 167      19.444  22.476  17.839  1.00 26.36           N  
    ATOM   1334  CA  ILE A 167      18.991  21.603  18.917  1.00 23.41           C  
    ATOM   1335  C   ILE A 167      20.002  21.550  20.075  1.00 22.93           C  
    ATOM   1336  O   ILE A 167      21.199  21.325  19.882  1.00 17.67           O  
    ATOM   1337  CB  ILE A 167      18.587  20.183  18.363  1.00 23.27           C  
    ATOM   1338  CG1 ILE A 167      17.330  20.317  17.489  1.00 17.76           C  
    ATOM   1339  CG2 ILE A 167      18.324  19.179  19.497  1.00 15.19           C  
    ATOM   1340  CD1 ILE A 167      16.918  19.059  16.763  1.00 18.86           C  
    ATOM   1341  N   THR A 168      19.483  21.799  21.273  1.00 19.63           N  
    ATOM   1342  CA  THR A 168      20.256  21.806  22.505  1.00 16.66           C  
    ATOM   1343  C   THR A 168      19.872  20.575  23.332  1.00 21.18           C  
    ATOM   1344  O   THR A 168      18.685  20.257  23.470  1.00 26.84           O  
    ATOM   1345  CB  THR A 168      19.954  23.118  23.296  1.00 19.13           C  
    ATOM   1346  OG1 THR A 168      20.374  24.247  22.519  1.00 23.07           O  
    ATOM   1347  CG2 THR A 168      20.648  23.147  24.656  1.00 15.66           C  
    ATOM   1348  N   ILE A 169      20.883  19.853  23.818  1.00 25.49           N  
    ATOM   1349  CA  ILE A 169      20.674  18.666  24.660  1.00 25.93           C  
    ATOM   1350  C   ILE A 169      21.214  18.979  26.067  1.00 23.26           C  
    ATOM   1351  O   ILE A 169      22.368  19.387  26.218  1.00 22.23           O  
    ATOM   1352  CB  ILE A 169      21.387  17.374  24.085  1.00 26.50           C  
    ATOM   1353  CG1 ILE A 169      20.953  17.081  22.639  1.00 18.48           C  
    ATOM   1354  CG2 ILE A 169      21.104  16.145  24.971  1.00 26.73           C  
    ATOM   1355  CD1 ILE A 169      19.475  16.804  22.407  1.00 22.35           C  
    ATOM   1356  N   THR A 170      20.346  18.857  27.073  1.00 25.01           N  
    ATOM   1357  CA  THR A 170      20.716  19.102  28.470  1.00 30.82           C  
    ATOM   1358  C   THR A 170      20.865  17.760  29.191  1.00 36.84           C  
    ATOM   1359  O   THR A 170      20.428  16.726  28.680  1.00 38.60           O  
    ATOM   1360  CB  THR A 170      19.665  19.977  29.211  1.00 24.07           C  
    ATOM   1361  OG1 THR A 170      18.425  19.271  29.324  1.00 38.57           O  
    ATOM   1362  CG2 THR A 170      19.419  21.269  28.467  1.00 30.76           C  
    ATOM   1363  N   GLN A 171      21.498  17.775  30.362  1.00 40.75           N  
    ATOM   1364  CA  GLN A 171      21.693  16.556  31.140  1.00 42.60           C  
    ATOM   1365  C   GLN A 171      20.543  16.315  32.101  1.00 44.01           C  
    ATOM   1366  O   GLN A 171      20.003  15.206  32.158  1.00 49.94           O  
    ATOM   1367  CB  GLN A 171      23.032  16.586  31.871  1.00 41.06           C  
    ATOM   1368  CG  GLN A 171      24.209  16.317  30.942  1.00 42.80           C  
    ATOM   1369  CD  GLN A 171      25.549  16.301  31.647  1.00 46.81           C  
    ATOM   1370  OE1 GLN A 171      25.876  17.207  32.415  1.00 57.96           O  
    ATOM   1371  NE2 GLN A 171      26.346  15.280  31.366  1.00 44.47           N  
    ATOM   1372  N   GLY A 172      20.142  17.371  32.806  1.00 40.89           N  
    ATOM   1373  CA  GLY A 172      19.044  17.274  33.753  1.00 47.41           C  
    ATOM   1374  C   GLY A 172      17.795  18.019  33.323  1.00 51.46           C  
    ATOM   1375  O   GLY A 172      16.756  17.881  34.003  1.00 58.08           O  
    ATOM   1376  OXT GLY A 172      17.846  18.747  32.309  1.00 51.83           O  
    TER    1377      GLY A 172                                                      
    HETATM 1378  O   HOH A2001      -5.014  21.998  32.622  1.00 34.52           O  
    HETATM 1379  O   HOH A2002       2.245  22.408  25.865  1.00 47.14           O  
    HETATM 1380  O   HOH A2003       2.409  21.472  22.890  1.00 48.26           O  
    HETATM 1381  O   HOH A2004      -4.451  19.531  12.824  1.00 50.15           O  
    HETATM 1382  O   HOH A2005       1.838  21.961  13.998  1.00 13.66           O  
    HETATM 1383  O   HOH A2006       1.829  14.585  15.655  1.00 50.45           O  
    HETATM 1384  O   HOH A2007       2.367  24.312   3.079  1.00 11.96           O  
    HETATM 1385  O   HOH A2008       1.024  36.141   4.489  1.00 32.84           O  
    HETATM 1386  O   HOH A2009      -4.778  24.446   3.666  1.00 50.35           O  
    HETATM 1387  O   HOH A2010      -2.356  30.207  -3.620  1.00 34.58           O  
    HETATM 1388  O   HOH A2011      -9.385  42.558  16.092  1.00 37.30           O  
    HETATM 1389  O   HOH A2012      -0.259  20.982   4.072  1.00 51.70           O  
    HETATM 1390  O   HOH A2013       2.611  21.006   9.109  1.00 13.98           O  
    HETATM 1391  O   HOH A2014       0.912  21.670  11.353  1.00 20.38           O  
    HETATM 1392  O   HOH A2015      -8.999  19.721  18.626  1.00 68.39           O  
    HETATM 1393  O   HOH A2016     -16.727  20.690  32.793  1.00 57.02           O  
    HETATM 1394  O   HOH A2017     -20.745  27.000  38.792  1.00 37.82           O  
    HETATM 1395  O   HOH A2018     -15.513  26.371  24.691  1.00 32.83           O  
    HETATM 1396  O   HOH A2019     -15.078  29.542  30.865  1.00 44.74           O  
    HETATM 1397  O   HOH A2020       0.252  16.309   9.227  1.00 41.08           O  
    HETATM 1398  O   HOH A2021      -9.990  20.817   3.821  1.00 49.82           O  
    HETATM 1399  O   HOH A2022      -0.419  37.766  27.386  1.00 54.16           O  
    HETATM 1400  O   HOH A2023      -7.892  40.755  18.155  1.00 46.94           O  
    HETATM 1401  O   HOH A2024     -11.761  39.790  15.088  1.00 48.83           O  
    HETATM 1402  O   HOH A2025     -11.346  34.386  23.127  1.00 10.20           O  
    HETATM 1403  O   HOH A2026     -14.109  29.473  17.376  1.00 43.27           O  
    HETATM 1404  O   HOH A2027       3.736   8.041  13.444  1.00 56.35           O  
    HETATM 1405  O   HOH A2028      -9.126  21.184  23.235  1.00 33.56           O  
    HETATM 1406  O   HOH A2029      -4.522  33.563  14.574  1.00 21.65           O  
    HETATM 1407  O   HOH A2030     -11.674  28.082  15.637  1.00 33.54           O  
    HETATM 1408  O   HOH A2031      -8.327  23.802   5.500  1.00 51.01           O  
    HETATM 1409  O   HOH A2032      13.157  26.356  10.881  1.00 33.55           O  
    HETATM 1410  O   HOH A2033      -5.113  31.695   9.120  1.00 22.31           O  
    HETATM 1411  O   HOH A2034       5.507  38.151  26.666  1.00 39.69           O  
    HETATM 1412  O   HOH A2035       0.870  32.559  34.070  1.00 43.70           O  
    HETATM 1413  O   HOH A2036       2.607  37.815  27.998  1.00 52.49           O  
    HETATM 1414  O   HOH A2037      12.213  31.559  19.628  1.00 29.69           O  
    HETATM 1415  O   HOH A2038      12.377  31.572   0.679  1.00 50.52           O  
    HETATM 1416  O   HOH A2039      13.092  26.561   5.159  1.00 18.72           O  
    HETATM 1417  O   HOH A2040       9.300  22.643  14.405  1.00 35.84           O  
    HETATM 1418  O   HOH A2041       6.401  18.979  23.245  1.00 41.82           O  
    HETATM 1419  O   HOH A2042       9.357  13.883  23.943  1.00 44.93           O  
    HETATM 1420  O   HOH A2043      15.561  28.801  24.966  1.00 87.04           O  
    HETATM 1421  O   HOH A2044      14.751  10.368  30.544  1.00 66.19           O  
    HETATM 1422  O   HOH A2045       2.102  14.985  11.341  1.00 51.57           O  
    HETATM 1423  O   HOH A2046      -6.623  20.356   2.521  1.00 50.09           O  
    HETATM 1424  O   HOH A2047      -1.649  24.600  -3.169  1.00 40.70           O  
    HETATM 1425  O   HOH A2048      -3.149  16.325  -2.674  1.00 42.69           O  
    HETATM 1426  O   HOH A2049       1.825  15.375  -1.224  1.00 36.36           O  
    HETATM 1427  O   HOH A2050      28.422  16.695  16.753  1.00 60.66           O  
    HETATM 1428  O   HOH A2051       9.519   9.301  -4.524  1.00 44.30           O  
    HETATM 1429  O   HOH A2052       5.900  12.242   6.948  1.00 37.46           O  
    HETATM 1430  O   HOH A2053      18.582   7.553  15.374  1.00 50.68           O  
    HETATM 1431  O   HOH A2054       6.598   8.667  12.058  1.00 26.23           O  
    HETATM 1432  O   HOH A2055      14.318   6.911  11.585  1.00 28.96           O  
    HETATM 1433  O   HOH A2056      21.727   9.020  18.950  1.00 28.27           O  
    HETATM 1434  O   HOH A2057      23.848  27.485  12.645  1.00 36.60           O  
    HETATM 1435  O   HOH A2058      24.101  21.378   9.626  1.00 32.76           O  
    HETATM 1436  O   HOH A2059      25.710  14.573  -1.360  1.00 32.48           O  
    HETATM 1437  O   HOH A2060      24.298  19.826   1.542  1.00 54.73           O  
    HETATM 1438  O   HOH A2061      19.942  15.837  -8.304  1.00 58.30           O  
    HETATM 1439  O   HOH A2062      13.980  24.390   8.699  1.00 24.13           O  
    HETATM 1440  O   HOH A2063      15.174  28.346  12.110  1.00 46.67           O  
    HETATM 1441  O   HOH A2064      23.857  24.793  16.730  1.00 42.96           O  
    HETATM 1442  O   HOH A2065      24.455  19.706  24.579  1.00 35.00           O  
    HETATM 1443  O   HOH A2066      23.153  19.975  31.160  1.00 31.01           O  
    HETATM 1444  O   HOH A2067      29.049  17.333  31.744  1.00 60.48           O  
    CONECT  203  532                                                                
    CONECT  532  203                                                                
    CONECT  855 1225                                                                
    CONECT 1225  855                                                                
    MASTER      270    0    0    2   20    0    0    6 1443    1    4   14          
    END                                                                             
    '''
    p = StructureAnalyser(121, structure, 'A','hello')

    print(p.get_structure_neighbours())
    print(p.get_superficiality())


# p=ProteinGatherer(uniprot='Q6ZN55').parse_uniprot().parse_pdb_blast()

# from protein.apriori_effect import WikiTable
# print(WikiTable(WikiTable.grantham).ndata)

def main():
    ## make everything!

    global_settings.error_tolerant = True

    ProteomeGatherer(skip=True, remake_pickles=True)

from protein.generate._proteome_gatherer2 import UniprotReader
import os, json
def mini_gene_data():
    genes = '''DOCK180
    DOCK2
    DOCK3
    DOCK4
    DOCK5
    DOCK6
    DOCK7
    DOCK8
    DOCK9
    DOCK10
    DOCK11
    '''.split()


    data = {}
    from pprint import PrettyPrinter
    pprint = PrettyPrinter().pprint
    namedex = json.load(open('data/human_prot_namedex.json'))
    for uni in set(namedex.values()):
        g = ProteinGatherer(uniprot=uni).parse_uniprot()
        data[g.gene_name] = {'name': g.gene_name, 'uniprot': g.uniprot, 'len': len(g), 'domains': {k: g.features[k] for k in ('active site','modified residue','topological domain','domain','region of interest','transmembrane region') if k in g.features}, 'disease': g.diseases}
        #print(g.gene_name,g.uniprot,len(g))
    json.dump(data,open('map.json','w'))

def make_pdb_dex():
    #I need to make a uniprot to pdb dex.
    from protein.generate._proteome_gatherer2 import UniprotReader
    master_file = os.path.join(ProteinGatherer.settings.temp_folder, 'uniprot_sprot.xml')
    UniprotReader.make_dictionary(uniprot_master_file=master_file, first_n_protein=0, chosen_attribute='uniprot')

def iterate_taxon(taxid):
    path = os.path.join(global_settings.pickle_folder,f'taxid{taxid}')
    for pf in os.listdir(path):
        protein = ProteinGatherer().load(file=os.path.join(path, pf))
        protein.get_offsets().parse_gNOMAD().compute_params()
        protein.dump()



if __name__ == '__main__':
    global_settings.verbose = True
    global_settings.init(data_folder='../test').retrieve_references(ask=False, refresh=False)
    #UniprotReader()

    #global_settings.init()

    #make_pdb_dex()
    iterate_taxon('9606')

    #p = ProteinGatherer(taxid='9606', uniprot='Q9BZ29').load().get_offsets()


    # fetch_binders is too slow. Pre-split the data like for gnomad.