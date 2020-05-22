import argparse
import sys
import os
import yaml
from pathlib import Path
from csv import DictReader
from operator import itemgetter, attrgetter
from itertools import groupby

from genetracks import Figure, Track, Multitrack, Coverage
import micall.core.plot_contigs as plot_contigs


def main(args):
    with open(args.blast_csv) as f:
        plot_blast(f, args.fasta)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'blast_csv',
        type=Path
    )
    parser.add_argument('fasta')
    args = parser.parse_args()
    return args


def get_contig_names(fasta):
    contig_nums = {}  # {contig_num: contig_name}
    with open(fasta) as f:
        for line in f:
            if line.startswith('>'):
                contig_name = line[1:-1]
                contig_nums[len(contig_nums) + 1] = contig_name
    return contig_nums


def get_landmark_tracks():
    tracks = []
    refname = 'HIV1-B-FR-K03455-seed'
    landmarks_path = (Path(os.path.realpath(__file__)).parent.parent / 'data' / 'landmark_references.yaml')
    landmark_groups = yaml.safe_load(landmarks_path.read_text())
    for reference_set in landmark_groups:
        if reference_set['coordinates'] != refname:
            continue
        prev_landmark = None
        for i, landmark in enumerate(sorted(reference_set['landmarks'],
                                            key=itemgetter('start'))):
            landmark.setdefault('frame', 0)
            if prev_landmark and 'end' not in prev_landmark:
                prev_landmark['end'] = landmark['start'] - 1
            prev_landmark = landmark
        for frame, frame_landmarks in groupby(reference_set['landmarks'],
                                              itemgetter('frame')):
            subtracks = []
            for landmark in frame_landmarks:
                subtracks.append(Track(landmark['start'],
                                       landmark['end'],
                                       label=landmark['name'],
                                       color=landmark['colour']))
            tracks.append(Multitrack(subtracks))
        break
    return tracks


def get_colors():
    colors = {
        'Intact': 'green',
	'5DEFECT': 'purple',
	'5DEFECT_GagNoATGGagFailed': 'purple',
	'5DEFECT_GagNoATGGagPassed': 'purple',
	'5DFECT_IntoGag': '',
	'Hypermut': 'red',
	'InternalInversion': 'purple',
	'LargeDeletion': 'orange',
	'NonHIV': 'purple',
	'PrematureStop_OR_AAtooLong_OR_AAtooShort': 'purple',
    }

    return colors


def plot_blast(blast_csv, fasta):
    padding = 2
    colors = get_colors()
    reader = DictReader(blast_csv)
    visited = set()
    figure = Figure(track_height=1)
    #figure.add(Track(1, 100, label='No contigs found', color='none'))
    contig_names = get_contig_names(fasta)
    with open('./contig_names', 'w') as f:
        yaml.dump(contig_names, f)
    
    min_start = 638
    max_end = 9604

    landmark_tracks = get_landmark_tracks()
    for track in landmark_tracks:
        figure.add(track)
    #for row in reader:
    #    print(contig_names[int(row['contig_num'])])
    multitracks = {
        'intact': {
            'data': [],
            'color': 'green'
        },
        'hypermut': {
            'data': [],
            'color': 'red'
        },
        'largedel': {
            'data': [],
            'color': 'orange'
        },
        'other_defect': {
            'data': [],
            'color': 'brown'
        }
    }
    for contig_num, group in groupby(reader, itemgetter('contig_num')):
        subtracks = []
        group = merge_group(group)
        for block in group:
            state = contig_names[int(contig_num)].split('::')[1]
            block_track = Track(
                block[0],
                block[1],
		color=colors[state],
                h=0.001
            )
            subtracks.append(block_track)
#            label_track = Track(
#                1,
#                9000,
#                label=contig_names[int(contig_num)],
#                color='none'
#            )
#            subtracks.append(label_track)
        multitrack = Multitrack(subtracks)
        if state == 'Intact':
            multitracks['intact']['data'].append(multitrack)
        elif state == 'Hypermut':
            multitracks['hypermut']['data'].append(multitrack)
        elif state == 'LargeDeletion':
            multitracks['largedel']['data'].append(multitrack)
        else:
            multitracks['other_defect']['data'].append(multitrack)
        #figure.add(multitrack)
    for _type in ('intact', 'hypermut', 'other_defect', 'largedel'):
        for track in multitracks[_type]['data']:
            figure.add(track, padding=padding)
    figure.show(w=970).saveSvg('./saved_figure.svg')


def merge_group(group):
    min_start = 638
    max_end = 9604
    target = 6623
    window = 50
    skip_next = False
    merged = []
    sorted_group = sorted(group, key=lambda x: int(x['ref_start']))
    for i, block in enumerate(sorted_group):
        if skip_next:
            skip_next = False
            continue
        if (
            int(block['ref_start']) < min_start
        ):
            continue
        try:
            next_block = sorted_group[i+1]
        except IndexError:
            pass
        if target-window <= int(block['ref_end']) <= target+window:
            new_block = (
                int(block['ref_start']),
                int(next_block['ref_end'])
            )
            skip_next = True
            merged.append(new_block)
            continue
        new_block = (
            int(block['ref_start']),
            min(int(block['ref_end']), max_end)
        )
        merged.append(new_block)
    return merged


if __name__ == '__main__':
    args = parse_args()
    main(args)
