"""
Prepare fasta file for AlphaFold pulldown.

Parameters:
- bait:           fasta file containing baait proteins
- target:         fasta file containing target proteins
- output_folder:  path to an existing folder where the results will be saved (one fasta per bait protein)
"""
import argparse
import logging
from pathlib import Path
import sys

from Bio.Seq import Seq
from Bio import SeqIO


logger = logging.getLogger(__name__)


def main():
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(processName)-10s (%(levelname)s) %(message)s')

    parser = argparse.ArgumentParser(description='Prepare fasta file for AlphaFold pulldown')
    parser.add_argument(
        '-b', '--bait', 
        type=Path,
        required=True,
        help='Path to fasta file containing bait proteins',
    )
    parser.add_argument(
        '-t', '--target', 
        type=Path,
        required=True,
        help='Path to fasta file containing target proteins',
    )
    parser.add_argument(
        '-o', '--output_folder', 
        type=Path,
        required=True,
        help='Path to existing output folder where the resuls will be saved',
    )
    args = parser.parse_args()

    bait_path = args.bait
    target_path = args.target
    output_folder = args.output_folder

    if not bait_path.is_file():
        logger.error(f'Bait fasta file does not exist: {bait_path}')
        sys.exit(1)
    elif not target_path.is_file():
        logger.error(f'Target fasta file does not exist: {target_path}')
        sys.exit(1)
    elif not output_folder.is_dir():
        logger.error(f'Output folder does not exist: {output_folder}')
        sys.exit(1)

    logger.info('Loading input fasta files')

    bait_dict = SeqIO.to_dict(SeqIO.parse(bait_path, 'fasta'))
    target_dict = SeqIO.to_dict(SeqIO.parse(target_path, 'fasta'))
    target_name = target_path.name.replace('.fasta', '').replace('.fa', '').replace('.faa', '')

    logger.info(f'Creating {len(bait_dict):,} pulldown fasta files')

    for i, (bait_id, bait_record) in enumerate(bait_dict.items()):
        output_path = output_folder / f'{bait_id}_{target_name}_pulldown.fasta'
        logger.info(f'Writing fasta file for bait {bait_id} to {output_path} ({i+1:,} / {len(bait_dict):,})')

        output_records = []
        bait_seq = str(bait_record.seq).upper()
        for target_id, target_record in target_dict.items():
            target_seq = str(target_record.seq).upper()

            seq_id = f'{bait_id}__{target_id}'
            seq = f'{bait_seq}:{target_seq}'

            record = SeqIO.SeqRecord(
                seq=Seq(seq),
                id=seq_id,
                name='',
                description='',
            )
            output_records.append(record)
        
        with output_path.open('w') as f_out:
            SeqIO.write(output_records, f_out, 'fasta')

    logger.info('DONE')
    sys.exit(0)


if __name__ == '__main__':
    main()
