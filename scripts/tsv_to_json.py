import pandas as pd
import json

# Read the file as a dataframe
df = pd.read_csv('../VUEs.txt', delimiter='\t', encoding='utf-8')

# Set columns to string type for splitting
df['pubmedId'] = df['pubmedId'].astype(str)
df['referenceText'] = df['referenceText'].astype(str)

# Group data by hugoGeneSymbol while preserving the order
grouped_data = df.groupby('hugoGeneSymbol', sort=False)

# Initialize the final JSON array
json_array = []

# Iterate over each group
for name, group in grouped_data:
    # Initialize the dictionary for each group
    group_dict = {
        "hugoGeneSymbol": name,
        "transcriptId": group['transcriptId'].iloc[0],
        "genomicLocationDescription": group['genomicLocationDescription'].iloc[0],
        "defaultEffect": group['defaultEffect'].iloc[0],
        "comment": group['comment'].iloc[0] if not pd.isnull(group['comment'].iloc[0]) else "",
        "context": group['context'].iloc[0] if not pd.isnull(group['context'].iloc[0]) else "",
        "revisedProteinEffects": []
    }

    # Iterate over each row in the group
    for index, row in group.iterrows():
        # Create a dictionary for each variant
        variant_dict = {
            "variant": row['variant'],
            "genomicLocation": row['genomicLocation'],
            "transcriptId": row['transcriptId'],
            "vepPredictedProteinEffect": row['vepPredictedProteinEffect'],
            "vepPredictedVariantClassification": row['vepPredictedVariantClassification'],
            "revisedProteinEffect": row['revisedProteinEffect'],
            "revisedVariantClassification": row['revisedVariantClassification'],
            "revisedStandardVariantClassification": row['revisedStandardVariantClassification'],
            "hgvsc": row['hgvsc'],
            "confirmed": bool(row['confirmed']),
            "references": []
        }
        if not pd.isnull(row['mutationOrigin']):
            variant_dict['mutationOrigin'] = row['mutationOrigin']
        if not pd.isnull(row['variantNote']):
            variant_dict['variantNote'] = row['variantNote']
        if not pd.isnull(row['otherVariation']):
            variant_dict['otherVariation'] = row['otherVariation']

        # Parse pubmedId and referenceText
        pubmed_ids = row['pubmedId'].split(';')
        reference_texts = row['referenceText'].split(';')

        # Create pairs of pubmedId and referenceText
        pairs = []
        for i in range(min(len(pubmed_ids), len(reference_texts))):
            pairs.append({
                "pubmedId": pubmed_ids[i],
                "referenceText": reference_texts[i]
            })

        # Add pairs to the reference list
        variant_dict['references'] = pairs

        # Add the variant dictionary to the revisedProteinEffects list
        group_dict['revisedProteinEffects'].append(variant_dict)

    # Add the group dictionary to the final JSON array
    json_array.append(group_dict)

# Write the JSON array to a file
with open('../generated/VUEs.json', 'w', encoding='utf-8') as f:
    json.dump(json_array, f, indent=4, ensure_ascii=False)
