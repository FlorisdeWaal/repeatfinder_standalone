#imports from std
import argparse
import re
#imports from project
import scorer as sc
import output_repeatfinder as out

def populate_known_motifs(keylist):
    outlist = []
    for k in keylist:
        outlist = k.split(":")[-1].strip()
    return outlist

def check_known_classes(seq,known_classes):
    #known classes are checked for by comparing the repeat containing sequences to these regex snippets
    #known classes are stored as regex patterns in a fasta like txt file
    qualdict = {}
    for key in known_classes:
        pattern = known_classes[key]
        pattern_instances = []
        matches = re.findall(known_classes[key],seq)

        if len(matches) is not 0: 
            qualdict[key+": "+pattern] = matches

    return qualdict

def get_known_classes_from_file():
    known_classes = {}
    infile = open("known_motifs.txt", "r")

    current_key = ""

    for li in infile:
        line = re.sub("\\n", "", li)
        line = re.sub("\\r", "", line)
        if line.startswith(">"):
            current_key = line[1:]

        else:
            known_classes[current_key] = line

    infile.close()

    return known_classes

def index_seq(seq,k = 3):
    kmer_dict = {}

    for i in range(len(seq)-k):
        word = str(seq[i:i+k])
        if word in kmer_dict:
            kmer_dict[word].append(i)

        else:
            kmer_dict[word] = [i]

    return kmer_dict

def find_top_kmer(kmer_dict):
    
    top_kmer_word = ""
    top_kmer_hits = []
    n_hits = 0
    for key in kmer_dict:
        if len(kmer_dict[key])>n_hits:
            top_kmer_word = key
            top_kmer_hits = kmer_dict[key]
            n_hits = len(kmer_dict[key])

    return top_kmer_word,top_kmer_hits 

def create_table(seq,top_kmer_word,top_kmer_hits, cutoff= 2.5, max_len = 15):
    table = [top_kmer_word]*len(top_kmer_hits)
    
    cont = True
    while cont:
        cont = expand_forward(seq,top_kmer_hits,table, cutoff, max_len)

    return table

def expand_forward(seq, hits, pattern_table, cutoff, max_len):

    next_residue_list = []

    for i, hit in enumerate(hits):

        if hit+len(pattern_table[i]) < len(seq) and len(pattern_table[i]) < max_len:
            next_residue_list.append(seq[hit+len(pattern_table[i])])
        else:
            return False
    if check_cutoff(next_residue_list, cutoff):

        for i, char in enumerate(next_residue_list):

            pattern_table[i] += char

        return True

    else:
        return False


def check_cutoff(nrl, cutoff):

    max_res = most_common(nrl)
    order, matrix = sc.blosum62()
    score_list = []

    for res in nrl:
        current_score = sc.score(max_res, res, order, matrix)
        score_list.append(current_score)

    if len(score_list) > 0:
        avg_score = sum(score_list)/len(score_list)
    else:
        return False

    if avg_score > cutoff:
        return True
    else:
        return False


def most_common(lst):

    try:
        return max(set(lst), key=lst.count)
    except:
        "print error occurred in most common function in repeatfinder.py"


def make_pattern(pattern_table):
    #function creates a regex from an expanded kmer table
    pattern_string = ""
    current_residue_list = []

    for i in range(len(pattern_table[0])):
        for pattern in pattern_table:
            current_residue_list.append(pattern[i])

        unique_residues = "".join(set(current_residue_list))

        if len(unique_residues) > 1:

            pattern_string += "[" + unique_residues + "]"

        else:
            pattern_string += unique_residues
        current_residue_list = []


    return pattern_string


def parse_fasta(iterator):

    if not isinstance(iterator, basestring):
        iterator = "\n".join(iterator)

    sequence_blocks = [x.strip() for x in iterator.split(">") if x.strip()]
    name = None
    sequences = {}

    for sequence_block in sequence_blocks:

        sequence_block = sequence_block.splitlines()
        name = sequence_block[0].strip()
        sequence = [x.strip() for x in sequence_block[1:] if x.strip()]
        sequence = "".join(sequence)
        sequences[name] = re.sub("\n","",sequence)

    return sequences

if __name__ == "__main__":
    parser = argparse.ArgumentParser("\ninput(-i) and output(-o) are essential for functioning of the script.\nExample use: ")
    parser.add_argument('-i',dest='infile',help='path to the input file(fasta format)')
    parser.add_argument('-o',dest= 'output', help='path to the output file')
    parser.add_argument('--maxlen', dest='maxlen', help='set maximum pattern length')
    parser.add_argument('--cutoff', dest='cutoff', help='set score cutoff for residue comparison')
    parser.add_argument('--kmer-size', dest='kmer', help='set the size of the initial kmer')

    args = parser.parse_args()

    kmer_size = 3
    cutoff = 2.5
    max_len = 15

    if args.kmer:
        kmer_size = args.kmer

    if args.maxlen:
        max_len = args.maxlen
    
    if args.cutoff:
        cutoff = args.cutoff

    with open(args.infile,"r") as infile:
        seqs = parse_fasta(infile) 
    known_motifs = get_known_classes_from_file()

    full_output_html = "<div>"
    full_output_html += "Settings: kmer_size: %s, max_len: %s, cutoff: %s <br>"%(kmer_size, max_len, cutoff)
    for key in seqs:
        km_occurences = check_known_classes(seqs[key],known_motifs)
        kmer_dict = index_seq(seqs[key], kmer_size)
        top_kmer_word,top_kmer_hits = find_top_kmer(kmer_dict)
        table = create_table(seqs[key],top_kmer_word,top_kmer_hits, cutoff, max_len)

        output_html = "sequence " + key + "<br>"
        has_km = False
        has_repeat = False
        should_print = False


        if len(table[0]) > 3 and len(table) > 2:  
            has_repeat = True
            should_print = True

        if len(km_occurences)>0:
            has_km = True
            should_print = True

        if has_repeat and has_km:
            pattern = make_pattern(table)
            output_html += "Pattern %s has been found %s times <br>"%(pattern, len(top_kmer_hits))
            output_html += "Instances: %s<br>"%("|".join(table)) 
            
            for kmkey in km_occurences:
                output_html += "Known motif %s has been found %s times <br>"%(kmkey, len(km_occurences[kmkey]))

            kms = populate_known_motifs(km_occurences.keys())
            output_html += "Sequence:<br>"
            output_html += out.format_linebreaks(out.generate_html_seq(seqs[key],pattern,kms),64)
        elif has_repeat and not has_km:
            pattern = make_pattern(table)
            output_html += "Pattern %s has been found %s times <br>"%(pattern, len(top_kmer_hits))
            output_html += "Instances: %s<br>"%("|".join(table)) 
            output_html += "Sequence:<br>"
            output_html += out.format_linebreaks(out.generate_html_seq(seqs[key],pattern),64)

        elif has_km and not has_repeat:
            for kmkey in km_occurences:
                output_html += "Known motif %s has been found %s times <br>"%(kmkey, len(km_occurences[kmkey]))

            kms = populate_known_motifs(km_occurences.keys())
            output_html += "Sequence:<br>"
            output_html += out.format_linebreaks(out.generate_html_seq(seqs[key],kms),64)

        output_html += "<hr>"
        if should_print:
            full_output_html += output_html
            

    full_output_html += "</div>"

    with open(args.output,"w") as outfile:
        outfile.writelines(full_output_html)
