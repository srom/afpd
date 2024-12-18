"""
Score structures produced by AlphaFold 3.
"""
import argparse
import json
import logging
from pathlib import Path
import sys
from typing import Dict, List

import pandas as pd


logger = logging.getLogger(__name__)


def main():
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(processName)-10s (%(levelname)s) %(message)s')

    parser = argparse.ArgumentParser(description='Score protein complexes docked with AlphaFold multimer')
    parser.add_argument(
        '-i', '--af_folder', 
        type=Path,
        required=True,
        help='Path to AlphaFold 3 output folder containing scores.',
    )
    parser.add_argument(
        '-o', '--output_path', 
        type=Path,
        required=True,
        help='Path to output CSV file.',
    )
    args = parser.parse_args()

    af_folder = args.af_folder
    output_path = args.output_path

    if not af_folder.is_dir():
        logger.error(f'AlphaFold 3 predictions folder does not exist: {af_folder}')
        sys.exit(1)
    elif not output_path.parent.is_dir():
        logger.error(f'Output folder does not exist: {output_path.parent}')
        sys.exit(1)

    logger.info('Score structures docked with AlphaFold 3')
    logger.info(f'AlphaFold predictions folder : {af_folder.resolve().as_posix()}')
    logger.info(f'Output CSV path with scores  : {output_path.resolve().as_posix()}')

    scores_paths = load_scores_paths(af_folder)

    logger.info(f'Number of results found: {len(scores_paths):,}')

    score_names = [
        'fraction_disordered',
        'has_clash',
        'iptm',
        'ptm',
        'ranking_score',
    ]
    scores_data = {
        'id': [],
    }
    for n in score_names:
        scores_data[n] = []

    for i, scores_path in enumerate(scores_paths):
        if i == 0 or (i+1) % 100 == 0 or (i+1) == len(scores_paths):
            logger.info(f'Scoring structure {i+1:,} / {len(scores_paths):,}')

        structure_id = scores_path.name.replace('_summary_confidences.json', '')
        scores_dict = read_scores_from_json_file(scores_path, score_names)
        
        scores_data['id'].append(structure_id)
        for n in score_names:
            scores_data[n].append(scores_dict.get(n))

    logger.info(f'Exporting sorted scores (best first) in CSV format to {output_path.resolve().as_posix()}')
    out_df = pd.DataFrame.from_dict(
        scores_data
    )
    out_df['confidence'] = (0.8 * out_df['iptm'] + 0.2 * out_df['ptm']).round(4)
    out_df.sort_values(
        'confidence', 
        ascending=False,
    ).to_csv(
        output_path,
        index=False,
    )

    logger.info('DONE')
    sys.exit(0)


def load_scores_paths(af_folder) -> List[Path]:
    paths = []
    for f in af_folder.glob('**/*_summary_confidences.json'):
        if f.is_file():
            paths.append(f)
    return paths


def read_scores_from_json_file(json_scores : Path, score_names : List[str]) -> Dict[str, float]:
    with open(json_scores, 'r') as f:
        scores_dict = json.load(f)

    return {
        name: scores_dict.get(name)
        for name in score_names
    }


if __name__ == '__main__':
    main()
