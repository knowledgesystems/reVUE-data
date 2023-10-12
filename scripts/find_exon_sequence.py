import requests
import sys

def get_exon_sequence(transcript_id, exon_number):
    server = "https://grch37.rest.ensembl.org"
    ext = f"/overlap/id/{transcript_id}?feature=exon"

    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})

    if not r.ok:
        r.raise_for_status()
        sys.exit()

    decoded = r.json()
    exons = [exon for exon in decoded if exon['feature_type'] == 'exon' and  exon['Parent'] == transcript_id ]
    exons.sort(key=lambda exon: exon['rank'])

    try:
        selected_exon = exons[exon_number - 1]  # -1 because list indices start at 0
    except IndexError:
        print(f"Transcript {transcript_id} does not have an exon number {exon_number}")
        sys.exit()

    sequence_ext = f"/sequence/id/{selected_exon['id']}?content-type=text/x-fasta"

    r = requests.get(server+sequence_ext, headers={ "Content-Type" : "application/json"})

    if not r.ok:
        r.raise_for_status()
        sys.exit()

    sequence = r.json()['seq']

    return sequence

# Usage:
# python find_exon_sequence.py ENST00000278616 7
transcript_id = sys.argv[1]
exon_number = int(sys.argv[2])

print(get_exon_sequence(transcript_id, exon_number))
