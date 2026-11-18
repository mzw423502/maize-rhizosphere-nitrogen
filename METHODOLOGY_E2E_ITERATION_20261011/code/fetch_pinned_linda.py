from pathlib import Path
from urllib.request import urlopen
import json, hashlib
root = Path(__file__).resolve().parents[1]
manifest = json.loads((root / 'upstream/CURRENT_LINDA_MANIFEST.json').read_bytes())
for row in manifest['files']:
    with urlopen(row['url']) as response:
        data = response.read()
    if hashlib.sha256(data).hexdigest() != row['sha256']:
        raise ValueError(row['path'])
    (root / row['path']).write_bytes(data)
