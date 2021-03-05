from Bio.SeqUtils import ProtParam, ProtParamData
from warnings import warn

# mod for DIWV
def mod(sequence):
    """
    This is a not implemented function. It is a fix for ProtParam.ProteinAnalysis().protein_scale and the DIWV scale.
    As the latter requires knowldge of the preceeding amino acid it will fail.
    >>> p = ProtParam.ProteinAnalysis(sequence)
    >>> p.protein_scale(ProtParamData.DIWV, window=9, edge=.4)
    hashtag epicfail.
    So this is the repalacement.
    :param sequence: sequence to score
    :type sequence: str
    :return: DIWV score.
    :rtype: list[int]
    """
    p = ProtParam.ProteinAnalysis(sequence)
    param_dict = ProtParamData.DIWV
    window = 9
    edge = 0.4
    weights = p._weight_list(window, edge)
    sum_of_weights = sum(weights) * 2 + 1
    scores = []
    for i in range(p.length - window):
        subsequence = p.sequence[i:i + window]
        score = 0.0
        for j in range(window // 2):
            try:
                front = param_dict[subsequence[j]][subsequence[j + 1]]
                back = param_dict[subsequence[window - j]][subsequence[window - j + 1]]
                score += weights[j] * front + weights[j] * back
            except KeyError:
                warn(f'warning: {subsequence[j]} or {subsequence[window - j - 1]} is not a standard amino acid.')
        middle = subsequence[window // 2]
        if middle in param_dict:
            score += param_dict[middle]
        else:
            warn(f'warning: {middle} is not a standard amino acid.')
        scores.append(score / sum_of_weights)
    return scores