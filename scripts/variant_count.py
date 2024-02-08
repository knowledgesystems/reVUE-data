import json
import os
import pandas as pd
from typing import List


# Download files first
# There will be internal MSK-IMAPCT, GENIE v15 public cohort, and TCGA Pan-Cancer Atlas (32 cohorts)
# Each one contains xxx_mutations.txt and xxx_clinical_samples.txt
# TCGA Pan-Cancer data can be downloaded by running download_files.py. Files are stored in /files/tcga
# MSK-IMAPCT data can be downloaded from https://github.mskcc.org/cdsi/msk-impact.git in msk_solid_heme folder. Files are stored in /files/mskimpact
# GENIE data can be downloaded from https://www.synapse.org/#!Synapse:syn53210170. Files are stored in /files/genie

class StudyInfo:
    def __init__(self, panel_name: str, panel_list: List[str]):
        self.panel_name = panel_name
        self.panel_list = panel_list

def find_file_pairs(folder_path: str):
    # Get all files in the directory
    files = os.listdir(folder_path)
    mutation_files = [f for f in files if f.endswith('_mutations.txt')]
    clinical_files = [f for f in files if f.endswith('_clinical_samples.txt')]

    # Find pairs of files that have the same prefix
    pairs = []
    for m_file in mutation_files:
        prefix = m_file.replace('_mutations.txt', '')
        for c_file in clinical_files:
            if c_file.startswith(prefix):
                pairs.append((m_file, c_file))

    return pairs

with open('../VUEs.json', 'r') as f:
    data = json.load(f)

def updateCounts(mutation_file: str, clinical_file: str, study_id: str, study_info: StudyInfo, multiple_study: bool, first_study: bool):
    jsonSet = {vue['genomicLocation'] for vueSet in data for vue in vueSet['revisedProteinEffects']}

    counts =  {'germlineVariantsCount': {}, 'somaticVariantsCount': {}, 'unknownVariantsCount': {}, 'totalPatientCount': 0, 'patientCountByGene': {}}
    mutationPatientsId = set()
    checkedPatientId = set()
    # Read clinical sample data
    clinical_sample_df = pd.read_csv(clinical_file, sep='\t')

    # Filter data by gene panel
    if study_info is not None:
        filtered_clinical_sample_df = clinical_sample_df[clinical_sample_df[study_info.panel_name].isin(study_info.panel_list)]
    else:
        filtered_clinical_sample_df = clinical_sample_df

    # Count unique PATIENT_ID
    counts['totalPatientCount'] = len(filtered_clinical_sample_df['PATIENT_ID'].unique())
    sample_to_patient_map = dict(zip(filtered_clinical_sample_df['SAMPLE_ID'], filtered_clinical_sample_df['PATIENT_ID']))
    mutations_df = pd.read_csv(mutation_file, sep='\t', usecols=['Hugo_Symbol', 'Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2', 'Mutation_Status', 'Tumor_Sample_Barcode'])

    filtered_mutations_df = mutations_df[mutations_df['Tumor_Sample_Barcode'].isin(sample_to_patient_map.keys())]

    # Count sample number by gene
    # Group the DataFrame by 'Hugo_Symbol' and count unique 'Tumor_Sample_Barcode' values
    for gene, group in filtered_mutations_df.groupby('Hugo_Symbol'):
        for index, row in group.iterrows():
            patientId = sample_to_patient_map.get(row['Tumor_Sample_Barcode'])
            if patientId and patientId not in checkedPatientId:
                counts['patientCountByGene'][gene] = counts['patientCountByGene'].get(gene, 0) + 1
                checkedPatientId.add(patientId)

    for index, row in filtered_mutations_df.iterrows():
        genomicLocationString = f"{row['Chromosome']},{row['Start_Position']},{row['End_Position']},{row['Reference_Allele']},{row['Tumor_Seq_Allele2']}"
        if genomicLocationString in jsonSet:
            if row['Tumor_Sample_Barcode'] in sample_to_patient_map and sample_to_patient_map[row['Tumor_Sample_Barcode']] not in mutationPatientsId:
                mutationPatientsId.add(sample_to_patient_map[row['Tumor_Sample_Barcode']])
                if row['Mutation_Status'].lower() == 'germline':
                    counts['germlineVariantsCount'][genomicLocationString] = counts['germlineVariantsCount'].get(genomicLocationString, 0) + 1
                elif row['Mutation_Status'].lower() == 'somatic':
                    counts['somaticVariantsCount'][genomicLocationString] = counts['somaticVariantsCount'].get(genomicLocationString, 0) + 1
                else:
                    counts['unknownVariantsCount'][genomicLocationString] = counts['unknownVariantsCount'].get(genomicLocationString, 0) + 1

    # Add the count numbers to the corresponding JSON objects
    for vueSet in data:
        for vue in vueSet['revisedProteinEffects']:
            if 'counts' in vue:
                if multiple_study and study_id in vue['counts'] and not first_study:
                    vue['counts'][study_id] = { 
                            "germlineVariantsCount": counts["germlineVariantsCount"].get(vue['genomicLocation'], 0) + vue['counts'][study_id]['germlineVariantsCount'],
                            "somaticVariantsCount": counts["somaticVariantsCount"].get(vue['genomicLocation'], 0) + vue['counts'][study_id]['somaticVariantsCount'],
                            "unknownVariantsCount": counts["unknownVariantsCount"].get(vue['genomicLocation'], 0) + vue['counts'][study_id]['unknownVariantsCount'],
                            "totalPatientCount": counts.get("totalPatientCount") + vue['counts'][study_id]['totalPatientCount'],
                            "genePatientCount": counts['patientCountByGene'].get(vueSet['hugoGeneSymbol'], 0) + vue['counts'][study_id]['genePatientCount']
                        }
                else:
                    vue['counts'][study_id] = { 
                            "germlineVariantsCount": counts["germlineVariantsCount"].get(vue['genomicLocation'], 0),
                            "somaticVariantsCount": counts["somaticVariantsCount"].get(vue['genomicLocation'], 0),
                            "unknownVariantsCount": counts["unknownVariantsCount"].get(vue['genomicLocation'], 0),
                            "totalPatientCount": counts.get("totalPatientCount"),
                            "genePatientCount": counts['patientCountByGene'].get(vueSet['hugoGeneSymbol'], 0)
                        }
            else:
                vue['counts']= { 
                    study_id: {
                        "germlineVariantsCount": counts["germlineVariantsCount"].get(vue['genomicLocation'], 0),
                        "somaticVariantsCount": counts["somaticVariantsCount"].get(vue['genomicLocation'], 0),
                        "unknownVariantsCount": counts["unknownVariantsCount"].get(vue['genomicLocation'], 0),
                        "totalPatientCount": counts.get("totalPatientCount"),
                        "genePatientCount": counts['patientCountByGene'].get(vueSet['hugoGeneSymbol'], 0)
                }}

# Update mskimpact counts
study_info_mskimpact = StudyInfo('GENE_PANEL', ["IMPACT341", "IMPACT410", "IMPACT468", "IMPACT505"])
updateCounts('./files/mskimpact/mskimpact_data_mutations_extended.txt', './files/mskimpact/mskimpact_data_clinical_sample.txt', 'mskimpact', study_info_mskimpact, False, False)
# Update genie counts
study_info_genie = StudyInfo('SEQ_ASSAY_ID', ["MSK-IMPACT341","MSK-IMPACT410", "MSK-IMPACT468", "MSK-IMPACT505"])
updateCounts('./files/genie/genie_data_mutations_extended.txt', './files/genie/genie_data_clinical_sample.txt', 'genie', study_info_genie, False, False)
# Update tcga counts
tcga_folder = './files/tcga'
tcga_file_pairs = find_file_pairs(tcga_folder)
for index, pair in enumerate(tcga_file_pairs):
    mutation_file_name, clinical_file_name = pair
    mutation_file_path = tcga_folder + '/' + mutation_file_name
    clinical_file_path = tcga_folder + '/' + clinical_file_name
    updateCounts(mutation_file_path, clinical_file_path, 'tcga', None, True, index == 0)

# Write the updated JSON data back to the file
with open('../VUEs.json', 'w') as f:
    json.dump(data, f, indent=4)
