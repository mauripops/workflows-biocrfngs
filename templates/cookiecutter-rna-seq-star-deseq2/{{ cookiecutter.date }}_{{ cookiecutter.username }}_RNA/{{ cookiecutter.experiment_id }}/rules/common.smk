def is_single_end(sample, unit):
    return pd.isnull(units.loc[(sample, unit), "fq2"])

def get_contrast(wildcards):
    return config["diffexp"]["contrasts"][wildcards.contrast]

def get_strandness(units):
    if "strandedness" in units.columns:
        return units["strandedness"].tolist()
    else:
        strand_list=["none"]
        return strand_list*units.shape[0]
