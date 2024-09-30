import pandas as pd
import json

# Read the file as a dataframe
df = pd.read_csv('../VUEs_table.txt', delimiter='\t', encoding='utf-8')

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
            "references": [],
            "counts": {}
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

        # Check and add counts if not empty
        if not pd.isnull(row['mskimpact_germlineVariantsCount']) and not pd.isnull(row['mskimpact_somaticVariantsCount']) and not pd.isnull(row['mskimpact_unknownVariantsCount']) and not pd.isnull(row['mskimpact_totalPatientCount']) and not pd.isnull(row['mskimpact_genePatientCount']):
            variant_dict['counts']['mskimpact'] = {
                "germlineVariantsCount": int(row['mskimpact_germlineVariantsCount']),
                "somaticVariantsCount": int(row['mskimpact_somaticVariantsCount']),
                "unknownVariantsCount": int(row['mskimpact_unknownVariantsCount']),
                "totalPatientCount": int(row['mskimpact_totalPatientCount']),
                "genePatientCount": int(row['mskimpact_genePatientCount'])
            }

        if not pd.isnull(row['tcga_germlineVariantsCount']) and not pd.isnull(row['tcga_somaticVariantsCount']) and not pd.isnull(row['tcga_unknownVariantsCount']) and not pd.isnull(row['tcga_totalPatientCount']) and not pd.isnull(row['tcga_genePatientCount']):
            variant_dict['counts']['tcga'] = {
                "germlineVariantsCount": int(row['tcga_germlineVariantsCount']),
                "somaticVariantsCount": int(row['tcga_somaticVariantsCount']),
                "unknownVariantsCount": int(row['tcga_unknownVariantsCount']),
                "totalPatientCount": int(row['tcga_totalPatientCount']),
                "genePatientCount": int(row['tcga_genePatientCount'])
            }

        if not pd.isnull(row['genie_germlineVariantsCount']) and not pd.isnull(row['genie_somaticVariantsCount']) and not pd.isnull(row['genie_unknownVariantsCount']) and not pd.isnull(row['genie_totalPatientCount']) and not pd.isnull(row['genie_genePatientCount']):
            variant_dict['counts']['genie'] = {
                "germlineVariantsCount": int(row['genie_germlineVariantsCount']),
                "somaticVariantsCount": int(row['genie_somaticVariantsCount']),
                "unknownVariantsCount": int(row['genie_unknownVariantsCount']),
                "totalPatientCount": int(row['genie_totalPatientCount']),
                "genePatientCount": int(row['genie_genePatientCount'])
            }

        if not pd.isnull(row['mskimpact_nonsignedout_germlineVariantsCount']) and not pd.isnull(row['mskimpact_nonsignedout_somaticVariantsCount']) and not pd.isnull(row['mskimpact_nonsignedout_unknownVariantsCount']) and not pd.isnull(row['mskimpact_nonsignedout_totalPatientCount']) and not pd.isnull(row['mskimpact_nonsignedout_genePatientCount']):
            variant_dict['counts']['mskimpact_nonsignedout'] = {
                "germlineVariantsCount": int(row['mskimpact_nonsignedout_germlineVariantsCount']),
                "somaticVariantsCount": int(row['mskimpact_nonsignedout_somaticVariantsCount']),
                "unknownVariantsCount": int(row['mskimpact_nonsignedout_unknownVariantsCount']),
                "totalPatientCount": int(row['mskimpact_nonsignedout_totalPatientCount']),
                "genePatientCount": int(row['mskimpact_nonsignedout_genePatientCount'])
            }

        # Add the variant dictionary to the revisedProteinEffects list
        group_dict['revisedProteinEffects'].append(variant_dict)

    # Add the group dictionary to the final JSON array
    json_array.append(group_dict)

# Write the JSON array to a file
with open('../VUEs.json', 'w', encoding='utf-8') as f:
    json.dump(json_array, f, indent=4, ensure_ascii=False)
