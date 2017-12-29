import numpy as np


def deseq2_normalization(list_of_profiles):
    """Perform DESeq2 like normalization position specific scores

    Parameters
    ----------
    list_of_profiles: array-like
        array of profiles across samples for one gene

    Returns
    -------
    normalized_profiles: array-like
        array of profiles across samples

    """

    # Transpose the profiles so that the rows represent the positions, the columns represent the samples
    list_of_profiles = np.array(list_of_profiles).T + 1
    geometric_mean_col = np.sqrt(list_of_profiles.prod(axis=1))
    profiles_ratio = list_of_profiles / geometric_mean_col[:, None]
    size_factors = np.median(profiles_ratio, axis=0)
    normalized_profiles = profiles_ratio / size_factors
    return normalized_profiles.T
