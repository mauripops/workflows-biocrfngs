# This file contains useful target rules to make calling snakemake easier

rule signal_tracks:
    """
    example usage: snakemake signal_tracks
    """
    input:
        expand("results/tracks/signal/{unit.sample}-{unit.unit}.bw", unit=units.itertuples())

rule interact_tracks:
    """
    example usage: snakemake interact_tracks
    """
    input:
        expand("results/tracks/interaction/{unit.sample}-{unit.unit}.bb", unit=units.itertuples())
