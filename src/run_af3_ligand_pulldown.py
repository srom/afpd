"""
Script to run a pulldown using one protein (or protein complex) as bait against a library of ligands.

This script assumes that the MSA part of the pipeline has been run already (i.e. AF3 run with flag --norun_inference)
Hence, a JSON file in AF3 format must already be available.

AlphaFold 3 (command `alphafold`) is assumed to be available in $PATH.

Arguments:
- Path to AF3 compatible JSON spec with MSAs
- Path to CSV containing ligands in SMILES format, along with an ID (must contain headers)
- CSV ID column name
- CSV SMILES column name
- Output folder

Generate a new JSON spec per ligand and runs AF3.

NOTE: a way to specify MSA as external files is available in the input template version 2. Not released yet.
TODO: use template version 2 when released.
"""
import argparse
import copy
import csv
import json
import logging
from pathlib import Path
import random
import shutil
import subprocess
import sys
import tempfile


logger = logging.getLogger(__name__)


def main():
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(processName)-10s (%(levelname)s) %(message)s')

    parser = argparse.ArgumentParser(description='AF3 ligand pulldown')
    parser.add_argument(
        '-i', '--spec_path', 
        type=Path,
        required=True,
        help='Path to AlphaFold 3 JSON spec file with MSA included.',
    )
    parser.add_argument(
        '-l', '--ligands_path', 
        type=Path,
        required=True,
        help='Path to CSV containing ligands in SMILES format.',
    )
    parser.add_argument(
        '--id_col', 
        type=str,
        required=True,
        help='Name of ID column',
    )
    parser.add_argument(
        '--smiles_col', 
        type=str,
        required=True,
        help='Name of column containg SMILES data',
    )
    parser.add_argument(
        '-o', '--output_folder', 
        type=Path,
        required=True,
        help='Path to output folder where structures will be saved.',
    )
    parser.add_argument(
        '--n_models', 
        type=int,
        required=False,
        default=1,
        help='Path to output folder where structures will be saved.',
    )
    args = parser.parse_args()

    spec_path = args.spec_path
    ligands_path = args.ligands_path
    id_col = args.id_col
    smiles_col = args.smiles_col
    output_folder = args.output_folder
    n_models = n_models

    if not spec_path.is_file():
        logger.error(f'Spec path does not exist: {spec_path}')
        sys.exit(1)
    elif not ligands_path.is_file():
        logger.error(f'Ligands path does not exist: {ligands_path}')
        sys.exit(1)
    elif not output_folder.is_dir():
        logger.error(f'Output folder does not exist: {output_folder}')
        sys.exit(1)

    with spec_path.open() as f:
        spec = json.load(f)

    ligands_dict = parse_ligands_csv(ligands_path, id_col, smiles_col)

    tempdir = Path(tempfile.mkdtemp())
    try:
        save_json_specs(spec, ligands_dict, tempdir, n_models)
        returncode = run_af3(tempdir, output_folder)
    finally:
        shutil.rmtree(tempdir)

    if returncode == 0:
        logger.info('DONE')
    else:
        logger.error('DONE with errors')

    sys.exit(returncode)


def save_json_specs(base_spec, ligands_dict, specs_dir, n_models):
    for ligand_id, ligand_smiles in ligands_dict.items():
        spec = copy.deepcopy(base_spec)
        spec['name'] = base_spec['name'] + f'__{ligand_id}'
        spec['modelSeeds'] = gen_model_seeds(n_models)

        ligand_spec = {
            'ligand': {
                'id': 'Z',
                'smiles': ligand_smiles,
            }
        }
        spec['sequences'].append(ligand_spec)

        with (specs_dir / f'{ligand_id}.json').open('w') as f_out:
            json.dump(spec, f_out, indent=True)


def run_af3(specs_dir, output_folder):
    result = subprocess.run(
        [
            'alphafold',
            '--input_dir', specs_dir.resolve().as_posix(),
            '--output_dir', output_folder.resolve().as_posix(),
            '--norun_data_pipeline',
        ],
        stdout=sys.stdout, 
        stderr=sys.stderr,
    )
    return result.returncode


def parse_ligands_csv(ligands_path, id_col, smiles_col):
    ligands_dict = {}
    with open(ligands_path, mode='r', newline='') as f:
        csv_reader = csv.reader(f)
        
        header = None
        id_index = None
        smiles_index = None
        for row in csv_reader:
            if header is None:
                header = row
                if id_col not in header:
                    raise ValueError(f'No ID column in header with name: {id_col}')
                elif smiles_col not in header:
                    raise ValueError(f'No SMILES column in header with name: {smiles_col}')
                
                id_index = header.index(id_col)
                smiles_index = header.index(smiles_col)
            else:
                ligands_dict[row[id_index]] = row[smiles_index].strip()

    return ligands_dict


def gen_model_seeds(n):
    return [int(random.uniform(1, 100)) for _ in range(n)]


if __name__ == '__main__':
    main()
