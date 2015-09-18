#!/usr/bin/env python3

import simulate_data

def main():
    names = ['arabis_short', 'arabis_large', 'HIV', 'choristoneura', 'mimivirus']
    files = ['arabis_mosaic_virus.fasta', 'arabis_mosaic_large.fasta', 'HIV-1.fasta', 'choristoneura.fasta', 'mimivirus.fasta']
    coverages = [25, 50, 75, 100, 150, 200, 250]

    for i in range(0, len(names)):
        for c in coverages:
            simulate_data.main(['--snp', '-c', str(c), '-m', '0.05', '-s', '0.0', '-n', '3', '--no_ref', '--seed', '1337', '--tmp', '/media/storage/tmp/', '../data/sequences/' + files[i], '../data/sht_' + names[i] + '_' + str(3*c) + '_005.tar.gz'])
            simulate_data.main(['--snp', '-c', str(c), '-m', '0.01', '-s', '0.0', '-n', '3', '--no_ref', '--seed', '1337', '--tmp', '/media/storage/tmp/', '../data/sequences/' + files[i], '../data/sht_' + names[i] + '_' + str(3*c)  + '_001.tar.gz'])
            simulate_data.main(['--snp', '-c', str(c), '-m', '0.1', '-s', '0.0', '-n', '3', '--no_ref', '--seed', '1337', '--tmp', '/media/storage/tmp/', '../data/sequences/' + files[i], '../data/sht_' + names[i] + '_' + str(3*c)  + '_01.tar.gz'])

if __name__ == '__main__':
    main()
