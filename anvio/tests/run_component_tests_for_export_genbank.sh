#!/bin/bash
source 00.sh

# Setup #############################
SETUP_WITH_OUTPUT_DIR $1 $2 $3
#####################################

INFO "Setting up the anvi-export-genbank test directory"
mkdir $output_dir/export_genbank_test
cp $files/data/metagenomes/human_gut/IGD_SUBSET/CONTIGS.db $output_dir/export_genbank_test

cd $output_dir/export_genbank_test

INFO "Migrating the contigs database"
anvi-migrate CONTIGS.db --migrate-quickly

INFO "Running anvi-export-genbank"
anvi-export-genbank -c CONTIGS.db \
                    -o contigs.gbk \
                    --min-contig-length 1000

INFO "Validating the GenBank output"
python - <<'PY'
from Bio import SeqIO

records = list(SeqIO.parse('contigs.gbk', 'genbank'))
if not records:
    raise SystemExit("No GenBank records were found in contigs.gbk")

cds_features = [feature for record in records for feature in record.features if feature.type == 'CDS']
if not cds_features:
    raise SystemExit("No CDS features were found in contigs.gbk")

missing_locus_tags = [feature for feature in cds_features if 'locus_tag' not in feature.qualifiers]
if missing_locus_tags:
    raise SystemExit(f"{len(missing_locus_tags)} CDS features were missing locus_tag qualifiers")

missing_protein_ids = [feature for feature in cds_features if 'protein_id' not in feature.qualifiers]
if missing_protein_ids:
    raise SystemExit(f"{len(missing_protein_ids)} CDS features were missing protein_id qualifiers")

missing_translations = [feature for feature in cds_features if 'translation' not in feature.qualifiers]
if missing_translations:
    raise SystemExit(f"{len(missing_translations)} CDS features were missing translation qualifiers")

bad_locus_tags = [feature.qualifiers['locus_tag'][0] for feature in cds_features if not feature.qualifiers['locus_tag'][0].startswith('anvio_gene_')]
if bad_locus_tags:
    raise SystemExit(f"Some locus_tag qualifiers did not preserve the anvi'o gene ID prefix: {bad_locus_tags[:5]}")

for feature in cds_features:
    locus_tag = feature.qualifiers['locus_tag'][0]
    protein_id = feature.qualifiers['protein_id'][0]
    if locus_tag != protein_id:
        raise SystemExit(f"locus_tag and protein_id differ for a CDS feature: {locus_tag} != {protein_id}")

print(f"Validated {len(records)} GenBank records and {len(cds_features)} CDS features.")
PY

INFO "Checking that a missing gene caller fails clearly"
if anvi-export-genbank -c CONTIGS.db -o should-not-exist.gbk --gene-caller THIS_GENE_CALLER_DOES_NOT_EXIST
then
    echo "anvi-export-genbank unexpectedly succeeded with a missing gene caller"
    exit 1
else
    echo "anvi-export-genbank failed as expected with a missing gene caller"
fi
