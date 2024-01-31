import pandas as pd
import requests
from Bio.Seq import Seq

def get_sequence(chr, start, end):
    server = "https://grch37.rest.ensembl.org"
    ext = f"/sequence/region/human/{chr}:{start}..{end}"
    headers={ "Content-Type" : "application/json"}
    
    r = requests.get(server+ext, headers=headers)

    if not r.ok:
        r.raise_for_status()
        return None

    decoded = r.json()
    return decoded['seq']

# translate sequence to AA
def translate_amino_acid_from_sequence(sequence):
    coding_dna = Seq(sequence)
    amino_acid = coding_dna.translate()
    return amino_acid

def is_inframe(whole_exon_skip, exon, start, exon_dict, transcript_id):
    exon_info = exon_dict[transcript_id]
    for e in exon_info['exons']:
        # find exon from dict
        if e['exon_number'] == exon:
            if whole_exon_skip:
                length = e['end'] - e['start'] + 1
                return length % 3 == 0
            else:
                length = e['end'] - start + 1
                return length % 3 == 0

def get_exon_by_number(exon_dict, exon_number, transcript_id):
    for exon in exon_dict[transcript_id]['exons']:
        if exon['exon_number'] == exon_number:
            return exon
    return None

def search_codon_position(sequence, codons):
    for i in range(0, len(sequence), 3):
        if sequence[i:i+3] in codons:
            return i // 3 + 1
    return None

def calculate_frameshift_end_codon(exon_dict, starting_exon, transcript_id, starting_sequence, chr):
    sequence = starting_sequence
    stop_codons = ['TAA', 'TAG', 'TGA']
    exon_number = starting_exon + 1

    while True:
        exon = get_exon_by_number(exon_dict, exon_number, transcript_id)
        if exon is None:
            break

        start = exon['start']
        end = exon['end']
        sequence += get_sequence(chr, start, end)

        end_codon_position = search_codon_position(sequence, stop_codons)
        if end_codon_position is not None:
            return end_codon_position - 1

        exon_number += 1

    return None

# calculate the start codon number of target exon, and how many nucleotides out of three are from the last codon of pre-exon
def get_exons_sequence_length(target_exon_number, transcript_id, exon_dict):
    exons = exon_dict[transcript_id]['exons']
    total_length = sum([e['end'] - e['start'] + 1 for e in exons if e['exon_number'] < target_exon_number])
    nt_number = total_length % 3
    target_exon_start_codon_number = total_length // 3 + 1
    return nt_number, target_exon_start_codon_number

# return genomic location of pre-exon, target-exon, post-exon
def get_exons_location(target_exon_number, transcript_id, exon_dict):
    exons = exon_dict[transcript_id]['exons']
    # pre-exon end genomic location
    pre_exon = next((exon for exon in exons if exon.get("exon_number") == target_exon_number - 1), None)
    # print(pre_exon)
    pre_exon_end = pre_exon.get("end") if pre_exon else None

    # target-exon start and end genomic location
    target_exon = next((exon for exon in exons if exon.get("exon_number") == target_exon_number), None)
    # print(target_exon)
    target_exon_start = target_exon.get("start") if target_exon else None
    target_exon_end = target_exon.get("end") if target_exon else None

    # post-exon start genomic location
    post_exon = next((exon for exon in exons if exon.get("exon_number") == target_exon_number + 1), None)
    # print(post_exon)
    post_exon_start = post_exon.get("start") if post_exon else None

    return pre_exon_end, target_exon_start, target_exon_end, post_exon_start


# get sequence for a single codon, it could have two parts from two exons (one from the end of last exon, the other is from the start of current exon)
def get_sequence_by_location(chr, frist_start, first_end, second_start, second_end):
    # print("frist_start:" + str(frist_start) )
    # print("first_end:" + str(first_end) )
    # print("second_start:" + str(second_start) )
    # print("second_end:" + str(second_end) )
    sequence1 = "" if frist_start is None or first_end is None or first_end < frist_start else get_sequence(chr, frist_start, first_end)
    sequence2 = "" if second_start is None or second_end is None or second_end < second_start else get_sequence(chr, second_start, second_end)
    # print("sequence1:" +sequence1)
    # print("sequence2:" + sequence2)
    return sequence1 + sequence2

# taking the first codon on a specific exon as the reference, calculate the codon number of giving genomic location on the exon
def calculate_codon_number_for_giving_location(start_genomic_location, start_codon_number, nt_number, target_genomic_location):
    # find the start location by concidering the nt number, create a sudo start for easy calculation
    sudo_start = start_genomic_location - nt_number
    if (target_genomic_location - sudo_start + 1) % 3 == 0:
        # when last AA is complete, we need to - 1
        target_condo_number = (target_genomic_location - sudo_start + 1) // 3 + start_codon_number - 1
    else:
        target_condo_number = (target_genomic_location - sudo_start + 1) // 3 + start_codon_number
     
    # print(target_genomic_location)
    # print(sudo_start)
    # print(start_codon_number)
    # print(target_genomic_location - sudo_start + 1)
    return target_condo_number

def correct_location(nt_number, is_start, frist_start, first_end, second_start, second_end):
    # nt_number == 0 means codon is complete on the exon
    if nt_number == 0: 
        if is_start == True:
            return None, None, second_start, second_end
        else:
            return first_end -2, first_end, None, None
    else:
        return frist_start, first_end, second_start, second_end

def get_exon_info(transcript_id, exon_dict):
    if transcript_id in exon_dict:
        return exon_dict[transcript_id]

    server = "https://grch37.rest.ensembl.org"
    ext = f"/lookup/id/{transcript_id}?expand=1&utr=1"
    r = requests.get(server+ext, headers={"Content-Type": "application/json"})

    if not r.ok:
        r.raise_for_status()
        return "Error occurred while fetching data."

    transcript = r.json()
    utrs = [feature for feature in transcript['UTR']]
    exons = [feature for feature in transcript['Exon']]


    five_prime_utr_end = max([utr['end'] for utr in utrs if utr['type'] == 'five_prime_utr'])
    three_prime_utr_start = min([utr['start'] for utr in utrs if utr['type'] == 'three_prime_utr'])
    # add exon number
    for i, exon in enumerate(exons):
        exon['exon_number'] = i + 1  # Exon numbers are 1-based
    
    # filters out exons that are entirely within the UTR regions
    exons = [exon for exon in exons if exon['end'] > five_prime_utr_end and exon['start'] < three_prime_utr_start]
    for exon in exons:
        # adjust exon start or end based on utr region
        if exon['start'] < five_prime_utr_end:
            exon['start'] = five_prime_utr_end + 1
        if exon['end'] > three_prime_utr_start:
            exon['end'] = three_prime_utr_start - 1

    exon_dict[transcript_id] = {"five_prime_utr_end": five_prime_utr_end, "three_prime_utr_start": three_prime_utr_start, "exons": exons}

    return exon_dict[transcript_id]

def calculate_protein_change(exon_dict, exon, transcript_id, start, end, whole_exon_skip, chr):

    is_inframe_mutation = is_inframe(whole_exon_skip, exon, start, exon_dict, transcript_id)

    if whole_exon_skip:
        if is_inframe_mutation:
            nt_number, target_exon_start_codon_number = get_exons_sequence_length(exon, transcript_id, exon_dict)
            pre_exon_end, target_exon_start, target_exon_end, post_exon_start = get_exons_location(exon, transcript_id, exon_dict)
            target_exon_end_codon_number = calculate_codon_number_for_giving_location(target_exon_start, target_exon_start_codon_number, nt_number, target_exon_end)
            start_codon_sequence = get_sequence_by_location(chr, *correct_location(nt_number, True, pre_exon_end - nt_number + 1, pre_exon_end, target_exon_start, target_exon_start + 2 - nt_number))
            end_codon_sequence = get_sequence_by_location(chr, *correct_location(nt_number, False, target_exon_end - nt_number + 1, target_exon_end, post_exon_start, post_exon_start + 2 - nt_number))
            start_codon_AA = translate_amino_acid_from_sequence(start_codon_sequence)
            end_codon_AA = translate_amino_acid_from_sequence(end_codon_sequence)
            printtext = "p." + start_codon_AA + str(target_exon_start_codon_number) + "_" + end_codon_AA + str(target_exon_end_codon_number) + "del"
            print(printtext)
        else:
            nt_number, target_exon_start_codon_number = get_exons_sequence_length(exon, transcript_id, exon_dict)
            pre_exon_end, target_exon_start, target_exon_end, post_exon_start = get_exons_location(exon, transcript_id, exon_dict)
            target_exon_end_codon_number = calculate_codon_number_for_giving_location(target_exon_start, target_exon_start_codon_number, nt_number, target_exon_end)
            # # if start of mutation is in the exon region
            # if start >= target_exon_start:
            #     # mutation is in the first codon of this exon and this first codon has complete 3 nt
            #     if nt_number == 0 and start <= target_exon_start + 2:
            #         mutation_codon_start = target_exon_start
            #         mutation_codon_end = target_exon_start + 2
            #         original_start_codon_sequence = get_sequence_by_location(chr, *correct_location(nt_number, True, pre_exon_end - nt_number + 1, pre_exon_end, mutation_codon_start, mutation_codon_end))

            #     # if mutation is in the first codon of this exon and this first codon does not have complete 3 nt
            #     elif nt_number != 0 and start <= target_exon_start + 2 - nt_number:
            #         mutation_codon_start = target_exon_start
            #         mutation_codon_end = target_exon_start + 2 - nt_number
            #         original_start_codon_sequence = get_sequence_by_location(chr, *correct_location(nt_number, True, pre_exon_end - nt_number + 1, pre_exon_end, mutation_codon_start, mutation_codon_end))

            #     # if mutation is not in the first codon
            #     else:
            #         mutation_codon_start_nt_number = (start - (target_exon_start - nt_number) + 1) % 3
            #         if mutation_codon_start_nt_number == 0:
            #             mutation_codon_start = start - 2
            #             mutation_codon_end = start
            #         elif mutation_codon_start_nt_number == 1:
            #             mutation_codon_start = start - 1
            #             mutation_codon_end = start + 1
            #         else:
            #             mutation_codon_start = start
            #             mutation_codon_end = start + 2
            #         original_start_codon_sequence = get_sequence_by_location(chr, *correct_location(nt_number, True, mutation_codon_start, mutation_codon_end, None, None))
            # else:
            #     # if start of mutation is in intron region
            #     original_start_codon_sequence = get_sequence_by_location(chr, *correct_location(nt_number, True, pre_exon_end - nt_number + 1, pre_exon_end, target_exon_start, target_exon_start + 2 - nt_number))
            print("nt_number:" + str(nt_number))
            start_codon_sequence_from_pre_exon_end = get_sequence_by_location(chr, *correct_location(nt_number, True, pre_exon_end - nt_number + 1, pre_exon_end, None, None))
            original_start_codon_sequence = get_sequence_by_location(chr, *correct_location(nt_number, True, pre_exon_end - nt_number + 1, pre_exon_end, target_exon_start, target_exon_start + 2 - nt_number))
            print("original_start_codon_sequence:" + original_start_codon_sequence)
            original_start_codon_AA = translate_amino_acid_from_sequence(original_start_codon_sequence)

            stop_codon_number = calculate_frameshift_end_codon(exon_dict, exon, transcript_id, start_codon_sequence_from_pre_exon_end, chr)
            print("stop_codon_number:" + str(stop_codon_number))

            new_codon_sequence = get_sequence_by_location(chr, *correct_location(nt_number, True, pre_exon_end - nt_number + 1, pre_exon_end, post_exon_start, post_exon_start + 2 - nt_number))
            print("new_codon_sequence:" + new_codon_sequence)
            new_codon_AA = translate_amino_acid_from_sequence(start_codon_sequence_from_pre_exon_end + new_codon_sequence)

            printtext = "p." + original_start_codon_AA + str(target_exon_start_codon_number) + new_codon_AA + "fs*" + str(stop_codon_number)
            print(printtext)

            

    pass

# Read the file
df = pd.read_csv('revue_bot_test.txt', sep='\t')



# Initialize the exon dictionary
exon_dict = {}

# Iterate over the rows of the dataframe
for index, row in df.iterrows():
    if row['strand'] == -1:
        # TODO: Filter out rows with strand equal to "-1"
        print(f"{row['gene']}_{row['transcript_id']} need manual review")
        # df = df[df['strand'] != -1]
    else:
        exon_info = get_exon_info(row['transcript_id'], exon_dict)
        calculate_protein_change(exon_dict, row['exon'], row['transcript_id'], row['start'], row['end'], row['whole_exon_skip'], row['chr'])
        # print(exon_info)
