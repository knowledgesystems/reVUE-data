import pandas as pd
import json
from pathlib import Path

def normalize_value(val):
    if pd.isna(val) or str(val).strip() == '':
        return None
    return str(val).strip()

def build_variant_dict(row):
    pubmed_ids = row['pubmedId'].split(';') if row['pubmedId'] else []
    reference_texts = row['referenceText'].split(';') if row['referenceText'] else []

    references = [
        {
            "pubmedId": pid,
            "referenceText": rtxt
        }
        for pid, rtxt in zip(pubmed_ids, reference_texts)
    ]

    variant = {
        "variant": row['variant'],
        "genomicLocation": row['genomicLocation'],
        "transcriptId": row['transcriptId'],
        "vepPredictedProteinEffect": row['vepPredictedProteinEffect'],
        "vepPredictedVariantClassification": row['vepPredictedVariantClassification'],
        "revisedProteinEffect": row['revisedProteinEffect'],
        "revisedVariantClassification": row['revisedVariantClassification'],
        "revisedStandardVariantClassification": row['revisedStandardVariantClassification'],
        "hgvsc": row['hgvsc'],
        "confirmed": str(row['confirmed']).lower() == 'true' if isinstance(row['confirmed'], str) else bool(row['confirmed']),
        "references": references
    }

    return variant


input_path = Path('../VUEs.txt')
output_path = Path('../generated/VUEs.json')

df = pd.read_csv(input_path, sep='\t', dtype=str).fillna('')
df = df.applymap(normalize_value)

json_output = []

for hugo, group in df.groupby('hugoGeneSymbol', sort=False):
    first = group.iloc[0]
    gene_entry = {
        "hugoGeneSymbol": hugo,
        "transcriptId": first['transcriptId'],
        "genomicLocationDescription": first['genomicLocationDescription'],
        "defaultEffect": first['defaultEffect'],
        "comment": first['comment'],
        "context": first['context'],
        "revisedProteinEffects": [build_variant_dict(row) for _, row in group.iterrows()]
    }
    json_output.append(gene_entry)

output_path.parent.mkdir(parents=True, exist_ok=True)
with open(output_path, 'w', encoding='utf-8') as f:
    json.dump(json_output, f, indent=4, ensure_ascii=False)
