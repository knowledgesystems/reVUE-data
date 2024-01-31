import requests
import sys
import json

def get_sequence(enst_id, query_type):
    if not enst_id:
        print("Please provide transcript id")
        return

    if not enst_id.startswith('ENST'):
        print("Transcript id is not valid")
        return
    server = "https://grch37.rest.ensembl.org"
    ext = f"/sequence/id/{enst_id}?content-type=text/x-fasta;type={query_type}"
    r = requests.get(server+ext, headers={"Content-Type": "application/json"})

    if not r.ok:
        r.raise_for_status()

    sequence_data = json.loads(r.text)
    return sequence_data['seq']

# python generate-vue-file.py ENST00000269305 protein
# python generate-vue-file.py ENST00000269305

if len(sys.argv) < 3:
    enst_id = sys.argv[1]
    query_type = "genomic"
else:
    enst_id = sys.argv[1]
    query_type = sys.argv[2]
sequence = get_sequence(enst_id, query_type)
print(sequence)