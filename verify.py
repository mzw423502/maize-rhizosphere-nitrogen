from pathlib import Path
import gzip
import hashlib
import sys

root = Path(__file__).resolve().parent
def verify(manifest, directory):
    entries = manifest.read_text(encoding='utf-8').splitlines()
    for entry in entries:
        expected, name = entry.split('  ', 1)
        with (directory / name).open('rb') as stream:
            if hashlib.file_digest(stream, 'sha256').hexdigest() != expected:
                raise ValueError(name)
    print(f'Verified {len(entries)} files from {manifest.name}.')

verify(root / 'SHA256SUMS-files.txt', root)
for experiment, expected in [('A', 3306), ('B', 2123), ('C', 93)]:
    path = root / '15_FINAL_AUDIT/PHASE5_D3_RANK_TRANSFER_20261010_V1/inputs' / f'{experiment}_COUNTS.tsv.gz'
    with gzip.open(path, 'rt') as stream:
        observed = sum(1 for line in stream) - 1
    if observed != expected:
        raise ValueError((experiment, observed, expected))
print('Empirical observations: 3,306 + 2,123 + 93 = 5,522.')
if len(sys.argv) > 1:
    verify(root / 'SHA256SUMS-assets.txt', Path(sys.argv[1]))
