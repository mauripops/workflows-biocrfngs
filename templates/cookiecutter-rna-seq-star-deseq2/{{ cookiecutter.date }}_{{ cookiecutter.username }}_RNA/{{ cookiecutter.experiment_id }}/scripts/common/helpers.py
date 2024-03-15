
# Get the condition for a given sample_id
def get_condition_for_sample(sample_id, samples):
    return samples.loc[(sample_id), 'condition']

def get_conditions(samples):
    return set([sample.condition for sample in samples.itertuples()])

def rgb2color(r,g,b):
    return '{},{},{}'.format(int(r*255),int(g*255),int(b*255))
