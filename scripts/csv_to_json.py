import pandas as pd
import requests
import json
import sys

# Get input and output file paths from command line arguments
input_file = sys.argv[1]
output_file = sys.argv[2]

# Step 1: Read TSV file and transform to DataFrame
df = pd.read_csv(input_file, sep='\t')
df = df[['gene', 'genotype', 'chr', 'start', 'end', 'ref', 'var', 'transcript_id', 'variant_classification', 'protein_change']]

# Step 2: Generate hgvsg_id and call Genome Nexus API
def generate_hgvsg_id(row):
    # Todo: add other type of hgvsg
    if pd.isnull(row['var']) or row['var'] == '' or row['var'] == '-':
        return f"{row['chr']}:g.{row['start']}_{row['end']}del"
    else:
        return f"{row['chr']}:g.{row['start']}{row['ref']}>{row['var']}"

df['hgvsg_id'] = df.apply(generate_hgvsg_id, axis=1)

def get_annotation(hgvsg_id):
    url = f"https://www.genomenexus.org/annotation/{hgvsg_id}?fields=annotation_summary"
    response = requests.get(url)
    data = response.json()
    if 'annotation_summary' in data and 'transcriptConsequenceSummary' in data['annotation_summary']:
        return data['annotation_summary']['transcriptConsequenceSummary'].get('variantClassification'), data['annotation_summary']['transcriptConsequenceSummary'].get('hgvspShort'), data['annotation_summary']['genomicLocation'].get('referenceAllele')
    else:
        return None, None, None

df['vep_predicted_variant_classification'], df['vep_predicted_protein_change'], df['ref'] = zip(*df['hgvsg_id'].map(get_annotation))

# Step 3: Generate output JSON file
output = []
for gene, group in df.groupby('gene'):
    variants = group.to_dict('records')
    revised_protein_effects = [{
        "variant": variant['hgvsg_id'],
        "genomicLocation": f"{variant['chr']},{variant['start']},{variant['end']},{variant['ref']},{variant['var']}",
        "vepPredictedProteinEffect": variant['vep_predicted_protein_change'],
        "vepPredictedVariantClassification": variant['vep_predicted_variant_classification'],
        "transcriptId": variant['transcript_id'],
        "revisedProteinEffect": variant['protein_change'],
        "revisedVariantClassification": variant['variant_classification'],
        "pubmedId": 31843900,
        "referenceText": "Casadei S et al., 2019",
        "confirmed": False
    } for variant in variants]
    output.append({
        "hugoGeneSymbol": gene,
        "transcriptId": variants[0]['transcript_id'],
        "genomicLocationDescription": "",
        "defaultEffect": "",
        "comment": "",
        "context": "",
        "revisedProteinEffects": revised_protein_effects
    })

with open(output_file, 'w') as f:
    json.dump(output, f)
