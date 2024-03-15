import gffutils

db = gffutils.create_db(snakemake.input[0],
                        dbfn=snakemake.output.db,
                        force=True,
                        keep_order=True,
                        merge_strategy='merge',
                        sort_attribute_values=True,
                        disable_infer_genes=False,
                        disable_infer_transcripts=False)

with open(snakemake.output.bed, 'w') as outfileobj:
    for tx in db.features_of_type('transcript', order_by='start'):
        bed = [s.strip() for s in db.bed12(tx).split('\t')]
        bed[3] = tx.id
        outfileobj.write('{}\n'.format('\t'.join(bed)))
