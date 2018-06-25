import pytest
from riboraptor.infer_protocol import infer_protocol


def test_inferprotocol():
    protocol, ppnn, pnpn = infer_protocol(
        'tests/data/SRX2536403_subsampled.unique.bam',
        'tests/data/hg38_v24_refseq.bed12')
    assert protocol == 'forward'


if __name__ == '__main__':
    pytest.main([__file__])
