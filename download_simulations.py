from pathlib import Path
from urllib.request import Request, urlopen
import argparse
import hashlib
import shutil

parser = argparse.ArgumentParser(description='Download and verify scenario archives from the v1.0.0 release.')
parser.add_argument('dest', type=Path)
parser.add_argument('scenarios', nargs='*', type=int, help='Optional scenario numbers, 0 through 25. Omit for all.')
args = parser.parse_args()
args.dest.mkdir(parents=True, exist_ok=True)
rows = (Path(__file__).resolve().parent / 'SHA256SUMS-assets.txt').read_text().splitlines()
expected = {line.split('  ', 1)[1]:line.split('  ', 1)[0] for line in rows}
def sha(path):
    with path.open('rb') as source:
        return hashlib.file_digest(source, 'sha256').hexdigest()
for cell in args.scenarios or range(26):
    if cell not in range(26):
        parser.error('scenario numbers must be between 0 and 25')
    name = f'simulation-c{cell:02d}.zip'
    path = args.dest / name
    if path.exists():
        if sha(path) != expected[name]:
            raise ValueError(f'Checksum mismatch in existing file: {path}')
        continue
    url = 'https://github.com/mzw423502/maize-rhizosphere-nitrogen/releases/download/v1.0.0/' + name
    part = path.with_suffix('.zip.partial')
    with urlopen(Request(url, headers={'User-Agent':'maize-rhizosphere-data/1.0'}), timeout=120) as response, part.open('wb') as out:
        shutil.copyfileobj(response, out)
    if sha(part) != expected[name]:
        raise ValueError(f'Checksum mismatch: {name}')
    part.replace(path)
    print(name, flush=True)
