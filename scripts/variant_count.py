import json
import os
import pandas as pd
from typing import List, Dict, Set, Any, Tuple
import requests

# Download files first
# There will be internal MSK-IMAPCT, GENIE v15 public cohort, and TCGA Pan-Cancer Atlas (32 cohorts)
# Each one contains xxx_mutations.txt and xxx_clinical_samples.txt
# TCGA Pan-Cancer data can be downloaded by running download_files.py. Files are stored in /files/tcga
# MSK-IMAPCT data can be downloaded from https://github.mskcc.org/cdsi/msk-impact.git in msk_solid_heme folder. Files are stored in /files/mskimpact
# GENIE data can be downloaded from https://www.synapse.org/#!Synapse:syn53210170. Files are stored in /files/genie

# Define gene panels
target_gene_panels_by_cohorts = {
    "mskimpact": ["IMPACT341", "IMPACT410", "IMPACT468", "IMPACT505"],
    "mskimpact_nonsignedout": ["IMPACT341", "IMPACT410", "IMPACT468", "IMPACT505"],
    # "genie": ["MSK-IMPACT341", "MSK-IMPACT410", "MSK-IMPACT468", "MSK-IMPACT505"]
}
vues = pd.read_json('../generated/VUEs.json')

vue_df = pd.json_normalize(
    data=vues.to_dict(orient='records'),
    record_path='revisedProteinEffects'
)

vue_df = vue_df.rename(columns={'genomicLocation': 'vue'})
vue_df = vue_df.set_index('vue')[[]]

# fetch therapeutic level from OncoKB
def get_therapeutic_level(genomic_location):
    url = f"https://www.oncokb.org/api/v1/annotate/mutations/byGenomicChange?genomicLocation={genomic_location}&referenceGenome=GRCh37"
    headers = {"Authorization": f"Bearer {os.environ.get('ONCOKB_TOKEN')}"}
    response = requests.get(url, headers=headers)
    if response.status_code == 200:
        return response.json().get('highestSensitiveLevel', None)
    return None

# Add the therapeutic levels
vue_df['therapeuticLevel'] = vue_df.index.to_series().apply(get_therapeutic_level)



def read_clinical_df(path, cohort_name, target_gene_panels_by_cohorts):
    """Read and normalize a clinical file. Filter by GENE_PANEL if gene panels are provided for this cohort."""
    df = pd.read_csv(path, sep="\t", low_memory=False)

    gene_panel_column = "SEQ_ASSAY_ID" if cohort_name == "genie" else "GENE_PANEL"
    if cohort_name in target_gene_panels_by_cohorts and gene_panel_column in df.columns:
        target_panels = target_gene_panels_by_cohorts[cohort_name]
        df = df[df[gene_panel_column].isin(target_panels)]

    df = df[['SAMPLE_ID', 'PATIENT_ID', 'CANCER_TYPE']].copy()
    df[['SAMPLE_ID', 'PATIENT_ID']] = df[['SAMPLE_ID', 'PATIENT_ID']].apply(
        lambda col: col.str.replace("GENIE-MSK-", "", regex=False)
    )
    return df

def read_mutation_df(path):
    df = pd.read_csv(path, sep="\t", low_memory=False)
    df = df[[
        'Chromosome', 'Start_Position', 'End_Position',
        'Reference_Allele', 'Tumor_Seq_Allele2',
        'Tumor_Sample_Barcode', 'Mutation_Status'
    ]].copy()

    df['Tumor_Sample_Barcode'] = df['Tumor_Sample_Barcode'].str.replace("GENIE-MSK-", "", regex=False)
    df['genomicLocation'] = df[['Chromosome', 'Start_Position', 'End_Position',
                                'Reference_Allele', 'Tumor_Seq_Allele2']].astype(str).agg(','.join, axis=1)

    return df

def build_variant_and_status_maps(mutation_df, vue_df):
    matched = mutation_df[mutation_df['genomicLocation'].isin(vue_df.index)]
    variant_to_sample = matched.groupby('genomicLocation')['Tumor_Sample_Barcode'].apply(list).to_dict()
    sample_to_mutation_status = matched.set_index('Tumor_Sample_Barcode')['Mutation_Status'].to_dict()
    return matched, variant_to_sample, sample_to_mutation_status

def build_sample_to_patient_map(clinical_df, valid_sample_ids):
    filtered = clinical_df[clinical_df['SAMPLE_ID'].isin(valid_sample_ids)]
    return filtered.set_index('SAMPLE_ID')['PATIENT_ID'].to_dict()

def process_cohort(clinical_df, mutation_df, vue_df):
    matched_df, variant_to_sample, sample_to_mutation_status = build_variant_and_status_maps(mutation_df, vue_df)
    sample_to_patient = build_sample_to_patient_map(clinical_df, matched_df['Tumor_Sample_Barcode'].unique())
    sample_to_cancer_type = clinical_df.set_index('SAMPLE_ID')['CANCER_TYPE'].to_dict()

    return variant_to_sample, sample_to_mutation_status, sample_to_patient, sample_to_cancer_type

def find_tcga_file_pairs(folder_path):
    files = os.listdir(folder_path)
    mutation_files = [f for f in files if f.endswith('_mutations.txt')]
    pairs = []
    for mut_file in mutation_files:
        prefix = mut_file.replace('_mutations.txt', '')
        clinical_file = f"{prefix}_clinical_samples.txt"
        if clinical_file in files:
            pairs.append((
                os.path.join(folder_path, clinical_file),
                os.path.join(folder_path, mut_file)
            ))
    return pairs

def read_all_tcga(folder_path):
    pairs = find_tcga_file_pairs(folder_path)

    all_clinical = [
        read_clinical_df(clinical_path, "tcga", target_gene_panels_by_cohorts)
        for clinical_path, _ in pairs
    ]
    all_mutation = [
        read_mutation_df(mutation_path)
        for _, mutation_path in pairs
    ]

    clinical_df = pd.concat(all_clinical, ignore_index=True)
    mutation_df = pd.concat(all_mutation, ignore_index=True)

    return clinical_df, mutation_df

def categorize_mutation_status(mutation_status_series):
    # gemrline, somatic, unknown
    return mutation_status_series.str.lower().where(mutation_status_series.str.lower().isin(['germline', 'somatic']), 'unknown')

def count_unique_patients_per_vue(vue, variant_to_sample, sample_to_mutation_status, sample_to_patient, sample_to_cancer_type):
    sample_ids = variant_to_sample.get(vue, [])
    df = pd.DataFrame({'sample_id': sample_ids})

    # mutation status category: gemrline, somatic, unknown
    raw_status = df['sample_id'].map(sample_to_mutation_status).fillna('unknown').str.lower()
    df['status'] = raw_status.where(raw_status.isin(['germline', 'somatic']), 'unknown')

    df['patient_id'] = df['sample_id'].map(sample_to_patient)
    df['cancer_type'] = df['sample_id'].map(sample_to_cancer_type)

    df = df.dropna(subset=['patient_id'])

    counts = {
        'germlineVariantsCount': 0,
        'somaticVariantsCount': 0,
        'unknownVariantsCount': 0,
        'germlineVariantsCountByCancerType': {},
        'somaticVariantsCountByCancerType': {},
        'unknownVariantsCountByCancerType': {}
    }

    for status in ['germline', 'somatic', 'unknown']:
        group = df[df['status'] == status].drop_duplicates('patient_id')
        counts[f'{status}VariantsCount'] = len(group)

        counts_by_cancer_type = group.groupby('cancer_type')['patient_id'].nunique()
        counts[f'{status}VariantsCountByCancerType'] = counts_by_cancer_type.to_dict()

    return counts

def add_vue_counts(vue_df,
                       variant_to_sample,
                       sample_to_mutation_status,
                       sample_to_cancer_type,
                       sample_to_patient,
                       clinical_df):
    vue_df = vue_df.copy()
    total_patients = clinical_df['PATIENT_ID'].nunique()

    vue_df['count'] = vue_df.index.to_series().apply(
        lambda vue: {
            **count_unique_patients_per_vue(
                vue,
                variant_to_sample,
                sample_to_mutation_status,
                sample_to_patient,
                sample_to_cancer_type
            ),
            'totalPatientCount': total_patients
        }
    )

    return vue_df

def merge_sample_maps(sample_maps_by_cohort):
    return dict(chain.from_iterable(d.items() for d in sample_maps_by_cohort.values()))

def update_vue_counts_json(vues_json, vue_df):
    for vue in vues_json:
        for effect in vue.get("revisedProteinEffects", []):
            vue_key = effect.get("genomicLocation")
            if vue_key in vue_df.index:
                counts_by_cohort = {}
                for col in vue_df.columns:
                    if col.startswith("count_"):
                        cohort = col.replace("count_", "")
                        count_data = vue_df.at[vue_key, col]
                        if isinstance(count_data, dict):
                            # Directly store the count structure as-is
                            counts_by_cohort[cohort] = count_data
                if "counts" in effect:
                    del effect["counts"]
                effect["counts"] = counts_by_cohort
    return vues_json


# Process each cohort
mskimpact_clinical_df = read_clinical_df("./files/mskimpact/mskimpact_data_clinical_sample.txt", "mskimpact", target_gene_panels_by_cohorts)
mskimpact_mutation_df = read_mutation_df("./files/mskimpact/mskimpact_data_mutations_extended.txt")
mskimpact_variant_to_sample, mskimpact_sample_to_mutation_status, mskimpact_sample_to_patient, mskimpact_sample_to_cancer_type = process_cohort(
    mskimpact_clinical_df,
    mskimpact_mutation_df,
    vue_df=vue_df
)

mskimpact_nonsignedout_clinical_df = read_clinical_df("./files/mskimpact_nonsignedout/data_clinical_sample.txt", "mskimpact_nonsignedout", target_gene_panels_by_cohorts)
mskimpact_nonsignedout_mutation_df = read_mutation_df("./files/mskimpact_nonsignedout/data_nonsignedout_mutations.txt")
mskimpact_nonsignedout_variant_to_sample, mskimpact_nonsignedout_sample_to_mutation_status, mskimpact_nonsignedout_sample_to_patient, mskimpact_nonsignedout_sample_to_cancer_type = process_cohort(
    mskimpact_nonsignedout_clinical_df,
    mskimpact_nonsignedout_mutation_df,
    vue_df=vue_df
)

genie_clinical_df = read_clinical_df("./files/genie/genie_data_clinical_sample.txt", "genie", target_gene_panels_by_cohorts)
genie_mutation_df = read_mutation_df("./files/genie/genie_data_mutations_extended.txt")
genie_variant_to_sample, genie_sample_to_mutation_status, genie_sample_to_patient, genie_sample_to_cancer_type = process_cohort(
    genie_clinical_df,
    genie_mutation_df,
    vue_df=vue_df
)

tcga_clinical_df, tcga_mutation_df = read_all_tcga("./files/tcga")
tcga_variant_to_sample, tcga_sample_to_mutation_status, tcga_sample_to_patient, tcga_sample_to_cancer_type = process_cohort(
    tcga_clinical_df,
    tcga_mutation_df,
    vue_df=vue_df
)

# Add counts
vue_df['count_mskimpact'] = add_vue_counts(
    vue_df=vue_df,
    variant_to_sample=mskimpact_variant_to_sample,
    sample_to_mutation_status=mskimpact_sample_to_mutation_status,
    sample_to_cancer_type=mskimpact_sample_to_cancer_type,
    sample_to_patient=mskimpact_sample_to_patient,
    clinical_df=mskimpact_clinical_df
    )['count']
# mskimpact_nonsignedout
vue_df['count_mskimpact_nonsignedout'] = add_vue_counts(
    vue_df=vue_df,
    variant_to_sample=mskimpact_nonsignedout_variant_to_sample,
    sample_to_mutation_status=mskimpact_nonsignedout_sample_to_mutation_status,
    sample_to_cancer_type=mskimpact_nonsignedout_sample_to_cancer_type,
    sample_to_patient=mskimpact_nonsignedout_sample_to_patient,
    clinical_df=mskimpact_nonsignedout_clinical_df
)['count']

# genie
vue_df['count_genie'] = add_vue_counts(
    vue_df=vue_df,
    variant_to_sample=genie_variant_to_sample,
    sample_to_mutation_status=genie_sample_to_mutation_status,
    sample_to_cancer_type=genie_sample_to_cancer_type,
    sample_to_patient=genie_sample_to_patient,
    clinical_df=genie_clinical_df
)['count']

# tcga
vue_df['count_tcga'] = add_vue_counts(
    vue_df=vue_df,
    variant_to_sample=tcga_variant_to_sample,
    sample_to_mutation_status=tcga_sample_to_mutation_status,
    sample_to_cancer_type=tcga_sample_to_cancer_type,
    sample_to_patient=tcga_sample_to_patient,
    clinical_df=tcga_clinical_df
)['count']

# Total counts
all_variant_to_sample = {
    "mskimpact": mskimpact_variant_to_sample,
    "mskimpact_nonsignedout": mskimpact_nonsignedout_variant_to_sample,
    "genie": genie_variant_to_sample,
    "tcga": tcga_variant_to_sample
}

all_sample_to_mutation_status = {
    "mskimpact": mskimpact_sample_to_mutation_status,
    "mskimpact_nonsignedout": mskimpact_nonsignedout_sample_to_mutation_status,
    "genie": genie_sample_to_mutation_status,
    "tcga": tcga_sample_to_mutation_status
}

all_sample_to_patient = {
    "mskimpact": mskimpact_sample_to_patient,
    "mskimpact_nonsignedout": mskimpact_nonsignedout_sample_to_patient,
    "genie": genie_sample_to_patient,
    "tcga": tcga_sample_to_patient
}

all_sample_to_cancer_type = {
    "mskimpact": mskimpact_sample_to_cancer_type,
    "mskimpact_nonsignedout": mskimpact_nonsignedout_sample_to_cancer_type,
    "genie": genie_sample_to_cancer_type,
    "tcga": tcga_sample_to_cancer_type
}

total_variant_to_sample = {}
for cohort_map in all_variant_to_sample.values():
    for vue, samples in cohort_map.items():
        total_variant_to_sample.setdefault(vue, []).extend(samples)

total_sample_to_mutation_status = merge_sample_maps(all_sample_to_mutation_status)
total_sample_to_patient = merge_sample_maps(all_sample_to_patient)
total_sample_to_cancer_type = merge_sample_maps(all_sample_to_cancer_type)
total_clinical_df = pd.concat([
    mskimpact_clinical_df[['PATIENT_ID']],
    mskimpact_nonsignedout_clinical_df[['PATIENT_ID']],
    genie_clinical_df[['PATIENT_ID']],
    tcga_clinical_df[['PATIENT_ID']]
], ignore_index=True).drop_duplicates()

vue_df['count_total'] = add_vue_counts(
    vue_df=vue_df,
    variant_to_sample=total_variant_to_sample,
    sample_to_mutation_status=total_sample_to_mutation_status,
    sample_to_cancer_type=total_sample_to_cancer_type,
    sample_to_patient=total_sample_to_patient,
    clinical_df=total_clinical_df
)['count']


with open("../generated/VUEs.json", "r", encoding="utf-8") as f:
    original_vues = json.load(f)

with open("../generated/VUEs.json", "w", encoding="utf-8") as f:
    json.dump(update_vue_counts_json(original_vues, vue_df), f, indent=2, ensure_ascii=False)
