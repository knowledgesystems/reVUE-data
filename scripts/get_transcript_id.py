import requests
import json
import sys

def get_ENST_id(nm_id):
    server = "https://grch37.rest.ensembl.org"
    ext = f"/xrefs/symbol/homo_sapiens/{nm_id}?content-type=application/json"
    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})

    if not r.ok:
        r.raise_for_status()
        return

    decoded = r.json()
    for item in decoded:
        if 'ENST' in item['id']:
            return item['id']
        
    ids = json.loads(r.text)
    for id in ids:
        if id['type'] == "transcript":
            return id['id'].id

# python get_transcript_id.py NM_000051
nm_id = sys.argv[1]
transcript_id = get_ENST_id(nm_id)
print(transcript_id)
