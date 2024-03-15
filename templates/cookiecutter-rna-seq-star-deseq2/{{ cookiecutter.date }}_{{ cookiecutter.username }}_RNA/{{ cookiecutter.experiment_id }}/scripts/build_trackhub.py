import trackhub
import os
import seaborn as sns
import common.helpers as helpers

# hack to get bigInteract working for now
trackhub.constants.track_typespecific_fields['bigInteract'] = trackhub.constants.track_typespecific_fields['bigBed']

# print("Snakemake parameters")
# print(snakemake.params)

config = snakemake.config
units = snakemake.params.units
samples = snakemake.params.samples

conditions = helpers.get_conditions(samples)
colors = [helpers.rgb2color(*color) for color in sns.color_palette("hls", len(conditions))]
# print(colors)

colormap = dict()
condition_mapping = dict()
mapping = dict()
for i, condition in enumerate(conditions):
    colormap[condition] = colors[i]
    condition_mapping[condition] = str(i)
    mapping[str(i)] = condition

# First we initialize the components of a track hub
hub, genome_file, genome, trackdb = trackhub.default_hub(
    hub_name=config["trackhub"]["name"],
    short_label=config["workflow"]["name"],
    long_label=config["workflow"]["name"],
    genome=config["ref"]["version"],
    email="biocrfngs@ust.hk")

subgroups = [
    # A contrived subgroup that will only tag it as "yes" if it's sample 1.
    trackhub.SubGroupDefinition(
        name='condition',
        label='Condition',
        mapping=mapping
    ),

    # While the two different views we create below are a good way of
    # turning on/off the signal or regions in bulk, this subgroup allows us
    # to sort the tracks by "name" and then by "kind". This is helpful for
    # ChIP-seq experiments where you want to have peaks under the
    # corresponding signal.
    trackhub.SubGroupDefinition(
        name='kind',
        label='Kind',
        mapping={
            'signal': 'signal',
            'interact': 'interact',
        }
    ),

]
#
# Create the composite track
composite = trackhub.CompositeTrack(
    name='composite',
    short_label='Signal and interact tracks',

    # The available options for dimensions are the `name` attributes of
    # each subgroup. Start with dimX and dimY (which become axes of the
    # checkbox matrix to select tracks), and then dimA, dimB, etc.
    dimensions='dimX=condition dimY=kind',

    # This enables a drop-down box under the checkbox matrix that lets us
    # select whatever dimA is (here, "kind").
    filterComposite='dimA',

    # The availalbe options here are the `name` attributes of each subgroup.
    sortOrder='condition=+ kind=-',
    tracktype='bed 3',
    visibility='full',
)

# Add those subgroups to the composite track
composite.add_subgroups(subgroups)

# Add the composite track to the trackDb
trackdb.add_tracks(composite)

# CompositeTracks compose different ViewTracks. We'll make one ViewTrack
# for signal, and one for bigBed regions.
signal_view = trackhub.ViewTrack(
    name='signalviewtrack',
    view='signal',
    visibility='full',
    tracktype='bigWig',
    short_label='Signal')

interact_view = trackhub.ViewTrack(
    name='interactviewtrack',
    view='interact',
    visibility='full',
    tracktype='bigInteract',
    short_label='interact_view')

# These need to be added to the composite.
composite.add_view(signal_view)
composite.add_view(interact_view)

for bigwig, unit in zip(snakemake.input.signals, units.itertuples()):
    print("Adding bigwig file {}".format(bigwig))
    sample_id = unit.sample
    condition = helpers.get_condition_for_sample(sample_id, samples)
    color = colormap[condition]
    name = trackhub.helpers.sanitize(os.path.basename(bigwig))
    short_label = sample_id.replace("_"," ")
    track_subgroup = {
        'kind': 'signal',
        'condition': condition_mapping[condition]
    }
    # track_subgroup[condition_mapping[condition]]
    track = trackhub.Track(
        name=name,
        shortLabel=short_label,
        longLabel=short_label,
        source=bigwig,
        visibility='full',
        tracktype='bigWig',
        color=color,
        autoScale='group',
        subgroups=track_subgroup,
        maxHeightPixels='100:32:28',
        windowingFunction='maximum',
        smoothingWindow=3,
        # darkerLabels='on'
    )

    signal_view.add_tracks(track)

for bigbed, unit in zip(snakemake.input.interacts, units.itertuples()):
    print("Adding bigbed file {}".format(bigbed))
    sample_id = unit.sample
    condition = helpers.get_condition_for_sample(sample_id, samples)
    color = colormap[condition]
    name = trackhub.helpers.sanitize(os.path.basename(bigbed))
    short_label = sample_id.replace("_"," ")
    track_subgroup = {
        'kind': 'interact',
        'condition': condition_mapping[condition]
    }
    track = trackhub.Track(
        name=name,
        shortLabel=short_label,
        longLabel=short_label,
        source=bigbed,
        subgroups=track_subgroup,
        visibility='full',
        maxHeightPixels="300:80:32",
        tracktype='bigInteract',
        # darkerLabels='on'
    )

    interact_view.add_tracks(track)


trackhub.upload.upload_hub(hub=hub, host='localhost', remote_dir=snakemake.output[0])
