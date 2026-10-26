"""Reproducible retrieval and alignment of a fixed public HCV HVR1 subset."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import re
import subprocess
import urllib.request
from .sequence import read_fasta, write_fasta

RECORD = 'https://zenodo.org/api/records/3266795'


def parse_frequency_fasta(path):
    """Retain frequency estimates from source headers; reject missing metadata."""
    records = read_fasta(path)
    headers = [line for line in path.read_text().splitlines() if line.startswith('>')]
    rows = []
    for header in headers:
        match = re.search(r'\bfreq=([0-9.eE+-]+)', header)
        if not match:
            raise ValueError('missing source haplotype frequency')
        frequency = float(match.group(1))
        if not math.isfinite(frequency) or frequency < 0:
            raise ValueError('invalid source frequency')
        name = header[1:].split()[0]
        rows.append((name, records[name], frequency))
    if not rows or sum(row[2] for row in rows) <= 0:
        raise ValueError('positive source frequency mass required')
    return rows


def fetch(folder):
    folder.mkdir(parents=True, exist_ok=True)
    metadata = json.load(urllib.request.urlopen(RECORD, timeout=60))
    files = sorted((f for f in metadata['files'] if re.fullmatch(
        r'\d+ chronic genotype 1 HCV infection\.fasta', f['key'])),
        key=lambda f: int(f['key'].split()[0]))[:8]
    if len(files) != 8:
        raise ValueError('expected eight eligible source specimens')
    provenance = {'record': 'https://zenodo.org/records/3266795',
                  'doi': '10.5281/zenodo.3266795', 'license': 'CC-BY-4.0',
                  'creator': 'Kamila Caraballo Cortes',
                  'selection': 'First eight chronic genotype 1 files sorted by numeric public sample identifier; chosen before running deconvolution.',
                  'files': []}
    (folder / 'raw').mkdir(exist_ok=True)
    for item in files:
        body = urllib.request.urlopen(item['links']['self'], timeout=60).read()
        if 'md5:' + hashlib.md5(body).hexdigest() != item['checksum']:
            raise ValueError('source checksum mismatch')
        source_id = item['key'].split()[0]
        local = f'raw/{source_id}.fasta'
        (folder / local).write_bytes(body)
        provenance['files'].append({'sample_id': source_id, 'original_filename': item['key'],
                                    'local_filename': local, 'url': item['links']['self'],
                                    'md5': item['checksum'], 'sha256': hashlib.sha256(body).hexdigest()})
    (folder / 'provenance.json').write_text(json.dumps(provenance, indent=2) + '\n')


def prepare(folder):
    """Remove source alignment gaps, then globally align all specimens together."""
    provenance = json.loads((folder / 'provenance.json').read_text())
    records, samples, mapping = {}, [], []
    for i, item in enumerate(provenance['files']):
        path = folder / item['local_filename']
        if hashlib.sha256(path.read_bytes()).hexdigest() != item['sha256']:
            raise ValueError('local source checksum mismatch')
        rows = parse_frequency_fasta(path)
        total = sum(f for _, _, f in rows)
        keys = []
        for j, (original, seq, frequency) in enumerate(rows):
            key = f's{i}_h{j}'
            records[key] = seq.replace('-', '')
            keys.append(key)
            mapping.append({'haplotype': key, 'source_sample': item['sample_id'],
                            'original_header_id': original, 'source_frequency_percent': frequency})
        samples.append({'id': i, 'source_sample': item['sample_id'], 'haplotypes': keys,
                        'frequencies': [f / total for _, _, f in rows]})
    write_fasta(folder / 'unaligned.fasta', records)
    version = subprocess.run(['mafft', '--version'], capture_output=True, text=True, check=True)
    with (folder / 'aligned.fasta').open('w') as output:
        subprocess.run(['mafft', '--auto', '--thread', '1', str(folder / 'unaligned.fasta')],
                             stdout=output, stderr=subprocess.PIPE, text=True, check=True)
    aligned = read_fasta(folder / 'aligned.fasta')
    if set(aligned) != set(records) or any(aligned[k].replace('-', '') != seq for k, seq in records.items()):
        raise ValueError('alignment changed sequence content or identifiers')
    manifest = {'synthetic': False, 'pooled_reads_simulated': True, 'titers_available': False,
                'alignment': {'program': version.stderr.strip() or version.stdout.strip(),
                              'arguments': ['--auto', '--thread', '1'],
                              'aligned_length': len(next(iter(aligned.values()))),
                              'sha256': hashlib.sha256((folder / 'aligned.fasta').read_bytes()).hexdigest()},
                'samples': samples, 'haplotype_mapping': mapping}
    (folder / 'samples.json').write_text(json.dumps(manifest, indent=2) + '\n')
    # Deliberately omit the full tool log: it may contain local filesystem paths.
    print(json.dumps({'samples': len(samples), 'haplotypes': len(records), 'alignment': manifest['alignment']}, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data', type=Path, default=Path('data/public_hcv'))
    parser.add_argument('--download', action='store_true', help='re-download pinned record subset')
    args = parser.parse_args()
    if args.download:
        fetch(args.data)
    prepare(args.data)


if __name__ == '__main__':
    main()
