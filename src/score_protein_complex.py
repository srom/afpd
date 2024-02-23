"""
Score protein complexes docked by AlphaFold multimer v3.

This script reports the following metrics:
- pLDDT  (1)
- pTM    (1)
- ipTM   (1)
- DockQ  (2)

Results are sorted by DockQ score (highest first) + ipTM to break even.

(1) Available directly from AlphaFold's output.
(2) For DockQ we report either:
    - pDockQ for protein complexes with 2 chains: [Bryant et al, March 2022](https://doi.org/10.1038/s41467-022-28865-w)
    - mpDockQ for protein complexes with >2 chains: [Bryant et al, October 2022](https://doi.org/10.1038/s41467-022-33729-4)

DockQ score implementations are adapted from [AlphaPulldown](https://github.com/KosinskiLab/AlphaPulldown).
"""
import argparse
import collections
import json
import logging
import math
from pathlib import Path
import re
import sys
from typing import List, Tuple

import numpy as np
import pandas as pd


logger = logging.getLogger(__name__)


def main():
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(processName)-10s (%(levelname)s) %(message)s')

    parser = argparse.ArgumentParser(description='Score protein complexes docked with AlphaFold multimer')
    parser.add_argument(
        '-i', '--af_folder', 
        type=Path,
        required=True,
        help='Path to AlphaFold output folder containing PDB structures and JSON scores.',
    )
    parser.add_argument(
        '-o', '--output_path', 
        type=Path,
        required=True,
        help='Path to output CSV file containing the protein complex scores.',
    )
    args = parser.parse_args()

    af_folder = args.af_folder
    output_path = args.output_path

    if not af_folder.is_dir():
        logger.error(f'AlphaFold predictions folder does not exist: {af_folder}')
        sys.exit(1)
    elif not output_path.parent.is_dir():
        logger.error(f'Output folder does not exist: {output_path.parent}')
        sys.exit(1)

    logger.info('Score protein complexes docked with AlphaFold multimer')
    logger.info(f'AlphaFold predictions folder : {af_folder.resolve().as_posix()}')
    logger.info(f'Output CSV path with scores  : {output_path.resolve().as_posix()}')

    protein_complex_files = load_protein_complex_files(af_folder)

    logger.info(f'Number of protein complexes found: {len(protein_complex_files):,}')

    scores_data = {
        'id'    : [],
        'plddt' : [],
        'ptm'   : [],
        'iptm'  : [],
        'dockq' : [],
    }
    for i, (complex_id, pdb_path, scores_path) in enumerate(protein_complex_files):
        if i == 0 or (i+1) % 100 == 0 or (i+1) == len(protein_complex_files):
            logger.info(f'Scoring protein complex {i+1:,} / {len(protein_complex_files):,}')

        plddt, plddt_avg, ptm, iptm = read_scores_from_json_file(scores_path)

        _, chain_coords, chain_CA_inds, chain_CB_inds = read_pdb(pdb_path)
        plddt_per_chain = read_plddt_per_chain(plddt, chain_CA_inds)
        complex_score, num_chains = score_complex(chain_coords, chain_CB_inds, plddt_per_chain)
        dockq_score = 0

        dockq_score = None
        if num_chains > 2:
            dockq_score = calculate_mpDockQ(complex_score)
        elif num_chains == 2:
            chain_coords, plddt_per_chain = read_pdb_pdockq(pdb_path)
            dockq_score = calc_pdockq(chain_coords,plddt_per_chain,t=8)

        scores_data['id'].append(complex_id)
        scores_data['plddt'].append(plddt_avg)
        scores_data['ptm'].append(ptm)
        scores_data['iptm'].append(iptm)
        scores_data['dockq'].append(np.round(dockq_score, 4))

    logger.info(f'Exporting sorted scores (best first) in CSV format to {output_path.resolve().as_posix()}')
    pd.DataFrame.from_dict(
        scores_data
    ).sort_values(
        ['dockq', 'iptm'], 
        ascending=False,
    ).to_csv(
        output_path,
        index=False,
    )

    logger.info('DONE')
    sys.exit(0)


def load_protein_complex_files(af_folder : Path) -> List[Tuple[str, Path, Path]]:
    """
    Iterate through AF predictions folder and find the PDB and JSON score files of the rank 1 model.
    """
    id_to_files = collections.defaultdict(dict)
    for p in af_folder.iterdir():
        if p.is_file():
            pdb_match = re.match(r'^(.+)_[^_]+_rank_001_.+.pdb$', p.name)
            scores_match = re.match(r'^(.+)_scores_rank_001_.+.json$', p.name)
            if pdb_match is not None:
                id_to_files[pdb_match[1]]['pdb'] = p
            elif scores_match is not None:
                id_to_files[scores_match[1]]['scores'] = p

    output = []
    for complex_id in id_to_files.keys():
        dct = id_to_files[complex_id]
        if 'pdb' not in dct:
            logger.warning(f'No PDB file found for complex {complex_id}. Skipping.')
        elif 'scores' not in dct:
            logger.warning(f'No JSON scores file found for complex {complex_id}. Skipping.')
        else:
            output.append((complex_id, dct['pdb'], dct['scores']))

    return output


def read_scores_from_json_file(json_scores : Path) -> Tuple[float, float, float]:
    with open(json_scores, 'r') as f:
        scores_dict = json.load(f)

    plddt = np.array(scores_dict.get('plddt', []), dtype=np.float32)
    plddt_avg = plddt.mean()
    ptm = scores_dict.get('ptm')
    iptm = scores_dict.get('iptm')

    return plddt, plddt_avg, ptm, iptm


###
# The functions below are either copied straight or lightly adapted from AlphaPulldown v1.0.4;
# https://github.com/KosinskiLab/AlphaPulldown/blob/1.0.4/alphapulldown/analysis_pipeline/calculate_mpdockq.py
###

def read_pdb(pdbfile):
    """
    Read a pdb file per chain.
    """
    pdb_chains = {}
    chain_coords = {}
    chain_CA_inds = {}
    chain_CB_inds = {}

    with open(pdbfile) as file:
        for line in file:
            if 'ATOM' in line:
                record = parse_atm_record(line)
                if record['chain'] in [*pdb_chains.keys()]:
                    pdb_chains[record['chain']].append(line)
                    chain_coords[record['chain']].append([record['x'],record['y'],record['z']])
                    coord_ind+=1
                    if record['atm_name'] == 'CA':
                        chain_CA_inds[record['chain']].append(coord_ind)
                    if record['atm_name'] == 'CB' or (record['atm_name'] == 'CA' and record['res_name'] == 'GLY'):
                        chain_CB_inds[record['chain']].append(coord_ind)
                else:
                    pdb_chains[record['chain']] = [line]
                    chain_coords[record['chain']]= [[record['x'],record['y'],record['z']]]
                    chain_CA_inds[record['chain']]= []
                    chain_CB_inds[record['chain']]= []
                    # Reset coord ind
                    coord_ind = 0

    return pdb_chains, chain_coords, chain_CA_inds, chain_CB_inds


def parse_atm_record(line):
    """
    Get the atm record.
    """
    record = collections.defaultdict()
    record['name'] = line[0:6].strip()
    record['atm_no'] = int(line[6:11].strip())
    record['atm_name'] = line[12:16].strip()
    record['atm_alt'] = line[17]
    record['res_name'] = line[17:20].strip()
    record['chain'] = line[21]
    record['res_no'] = int(line[22:26].strip())
    record['insert'] = line[26].strip()
    record['resid'] = line[22:29]
    record['x'] = float(line[30:38])
    record['y'] = float(line[38:46])
    record['z'] = float(line[46:54])
    record['occ'] = float(line[54:60])
    record['B'] = float(line[60:66])
    return record


def read_plddt_per_chain(plddt, chain_CA_inds):
    """
    Get the plDDT for each chain.
    """
    chain_names = chain_CA_inds.keys()
    chain_lengths = dict()
    for name in chain_names:
        curr_len = len(chain_CA_inds[name])
        chain_lengths[name] = curr_len
    
    plddt_per_chain = dict()
    curr_len = 0
    for k, v in chain_lengths.items():
        curr_plddt = plddt[curr_len:curr_len+v]
        plddt_per_chain[k] = curr_plddt
        curr_len += v

    return plddt_per_chain


def score_complex(path_coords, path_CB_inds, path_plddt):
    """
    Score all interfaces in the current complex.

    Modified from the score_complex() function in MoLPC repo: 
    https://gitlab.com/patrickbryant1/molpc/-/blob/main/src/complex_assembly/score_entire_complex.py#L106-154
    """

    chains = [*path_coords.keys()]
    chain_inds = np.arange(len(chains))
    complex_score = 0
    # Get interfaces per chain
    for i in chain_inds:
        chain_i = chains[i]
        chain_coords = np.array(path_coords[chain_i])
        chain_CB_inds = path_CB_inds[chain_i]
        l1 = len(chain_CB_inds)
        chain_CB_coords = chain_coords[chain_CB_inds]
        chain_plddt = path_plddt[chain_i]
 
        for int_i in np.setdiff1d(chain_inds, i):
            int_chain = chains[int_i]
            int_chain_CB_coords = np.array(path_coords[int_chain])[path_CB_inds[int_chain]]
            int_chain_plddt = path_plddt[int_chain]
            # Calc 2-norm
            mat = np.append(chain_CB_coords,int_chain_CB_coords,axis=0)
            a_min_b = mat[:,np.newaxis,:] -mat[np.newaxis,:,:]
            dists = np.sqrt(np.sum(a_min_b.T ** 2, axis=0)).T
            contact_dists = dists[:l1,l1:]
            contacts = np.argwhere(contact_dists <= 8)
            # The first axis contains the contacts from chain 1
            # The second the contacts from chain 2
            if contacts.shape[0] > 0:
                av_if_plDDT = np.concatenate((chain_plddt[contacts[:,0]], int_chain_plddt[contacts[:,1]])).mean()
                complex_score += np.log10(contacts.shape[0]+1)*av_if_plDDT

    return complex_score, len(chains)


def calculate_mpDockQ(complex_score):
    """
    A function that returns a complex's mpDockQ score after calculating complex_score.
    """
    L = 0.827
    x_0 = 261.398
    k = 0.036
    b = 0.221
    return L/(1+math.exp(-1*k*(complex_score-x_0))) + b


def read_pdb_pdockq(pdbfile):
    """
    Read a pdb file predicted with AF and rewritten to contain all chains.
    Adepted from FoldDock repo:
    https://gitlab.com/ElofssonLab/FoldDock/-/blob/main/src/pdockq.py#L34-59
    """
    chain_coords, chain_plddt = {}, {}
    with open(pdbfile, 'r') as file:
        for line in file:
            if not line.startswith('ATOM'):
                continue
            record = parse_atm_record(line)
            # Get CB - CA for GLY
            if record['atm_name']=='CB' or (record['atm_name']=='CA' and record['res_name']=='GLY'):
                if record['chain'] in [*chain_coords.keys()]:
                    chain_coords[record['chain']].append([record['x'],record['y'],record['z']])
                    chain_plddt[record['chain']].append(record['B'])
                else:
                    chain_coords[record['chain']] = [[record['x'],record['y'],record['z']]]
                    chain_plddt[record['chain']] = [record['B']]

    # Convert to arrays
    for chain in chain_coords:
        chain_coords[chain] = np.array(chain_coords[chain])
        chain_plddt[chain] = np.array(chain_plddt[chain])

    return chain_coords, chain_plddt


def calc_pdockq(chain_coords, chain_plddt, t):
    """
    Calculate the pDockQ scores
    pdockQ = L / (1 + np.exp(-k*(x-x0)))+b
    L= 0.724 x0= 152.611 k= 0.052 and b= 0.018

    Modified from the calc_pdockq() from FoldDock repo: 
    https://gitlab.com/ElofssonLab/FoldDock/-/blob/main/src/pdockq.py#L62
    """
    # Get coords and plddt per chain
    ch1, ch2 = [*chain_coords.keys()]
    coords1, coords2 = chain_coords[ch1], chain_coords[ch2]
    plddt1, plddt2 = chain_plddt[ch1], chain_plddt[ch2]

    # Calc 2-norm
    mat = np.append(coords1, coords2,axis=0)
    a_min_b = mat[:,np.newaxis,:] -mat[np.newaxis,:,:]
    dists = np.sqrt(np.sum(a_min_b.T ** 2, axis=0)).T
    l1 = len(coords1)
    contact_dists = dists[:l1,l1:] # upper triangular --> first dim = chain 1
    contacts = np.argwhere(contact_dists<=t)

    if contacts.shape[0] < 1:
        pdockq = 0
    else:
        # Get the average interface plDDT
        avg_if_plddt = np.average(np.concatenate([plddt1[np.unique(contacts[:,0])], plddt2[np.unique(contacts[:,1])]]))
        # Get the number of interface contacts
        n_if_contacts = contacts.shape[0]
        x = avg_if_plddt*np.log10(n_if_contacts)
        pdockq = 0.724 / (1 + np.exp(-0.052*(x-152.611)))+0.018

    return pdockq


if __name__ == '__main__':
    main()
