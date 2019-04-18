from ._protein_gatherer import ProteinGatherer
from ._proteome_gatherer import ProteomeGatherer

### Some reference files have to be fetched manually...
refs = (
    'ftp://ftp.broadinstitute.org/pub/ExAC_release/release1/ExAC.r1.sites.vep.vcf.gz',
    'http://geneontology.org/gene-associations/goa_human.gaf.gz',
    'http://purl.obolibrary.org/obo/go.obo',
    'http://interactome.baderlab.org/data/Raul-Vidal(Nature_2005).psi',
    'http://slorth.biochem.sussex.ac.uk/download/h.sapiens_ssl_predictions.csv',
    'https://downloads.thebiogrid.org/Download/BioGRID/Release-Archive/BIOGRID-3.5.166/BIOGRID-ALL-3.5.166.mitab.zip',
    'https://stringdb-static.org/download/protein.links.v10.5/9606.protein.links.v10.5.txt.gz',
    'http://www.ensembl.org/biomart/martview/436b8b2b06f64bbee960592afda10817?VIRTUALSCHEMANAME=default&ATTRIBUTES=hsapiens_gene_ensembl.default.feature_page.ensembl_gene_id|hsapiens_gene_ensembl.default.feature_page.ensembl_transcript_id|hsapiens_gene_ensembl.default.feature_page.external_gene_name|hsapiens_gene_ensembl.default.feature_page.uniprotswissprot|hsapiens_gene_ensembl.default.feature_page.uniprot_gn|hsapiens_gene_ensembl.default.feature_page.ensembl_peptide_id&FILTERS=&VISIBLEPANEL=resultspanel',
    'ftp://ftp.nextprot.org/pub/current_release/mapping/nextprot_refseq.txt',
    'https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.exome_calling_intervals.sites.vcf.bgz'
    'ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/pdb_chain_uniprot.tsv.gz')


failed_blast_pdbs ='A0A0K0K1G8'
