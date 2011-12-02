# parses KO file into a "otu fasta-like" format
# each position in the output fasta sequences represents the
# presence/absence of a particular KO (sorted alphabetically)
import sys
from os.path import split, splitext

alphabet = {'present':'T', 'absent':'A'}

if __name__ == "__main__":
    fp_in = sys.argv[1]
    fp_out_fasta = splitext(split(fp_in)[-1])[0] + '.fna'
    fp_out_ko_list = splitext(split(fp_in)[-1])[0] + '_ko_list.txt'
    print 'Parsing input file:', fp_in
    print 'Writing output files:', fp_out_fasta, fp_out_ko_list

    ko_lines = open(fp_in,'U').readlines()
    genomes = dict()
    kos = set()
    
    # parse input
    for line in ko_lines:
        words = line.strip().split('\t')
        ko_id = words[0]
        genome_id = words[1]
        if not genomes.has_key(genome_id):
            genomes[genome_id] = {}
        genomes[genome_id][ko_id] = True
        kos.add(ko_id)
        
    # write output in "fasta-like" format
    # output fasta is <input filename>.fna
    f_out_fasta = open(fp_out_fasta,'w')
    kos = sorted(kos)
    for genome_id in sorted(genomes.keys()):
        f_out_fasta.write('>' + genome_id + '\n')
        genome_kos = genomes[genome_id].keys()
        genome_ko_string = ''
        for ko in kos:
            if ko in genome_kos:
                genome_ko_string += alphabet['present']
            else:
                genome_ko_string += alphabet['absent']
        f_out_fasta.write(genome_ko_string + '\n')
    f_out_fasta.close()

    # write list of KOs in same order as fasta-like output sequences
    f_out_ko_list = open(fp_out_ko_list,'w')
    f_out_ko_list.write('\n'.join(kos))
    f_out_ko_list.close()
