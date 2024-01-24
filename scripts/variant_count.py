import pandas as pd
import json
import glob


### read from multiple files
# Get a list of all TXT files in the folder
files = glob.glob('./files/*.txt')

### read from single file
# Read the file into a DataFrame, 
# df = pd.read_csv('internal_msk_impact.txt',  sep="\t", usecols=['Hugo_Symbol', 'Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2', 'Mutation_Status'])

# Create maps for variants count
germlineVariantsCount = {}
somaticVariantsCount = {}
unknownVariantsCount = {}

### start
# Read the JSON file
with open('../VUEs.json', 'r') as f:
    data = json.load(f)

# Create a set from the JSON data
jsonSet = {vue['genomicLocation'] for vueSet in data for vue in vueSet['revisedProteinEffects']}

# Iterate over each file
for file in files:
    # Read each TXT file into a DataFrame and append it to a list
    df = pd.read_csv(file, sep="\t", usecols=['Hugo_Symbol', 'Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2', 'Mutation_Status'])
    
    print("Counting from " + file)
    # Iterate over the DataFrame rows
    for index, row in df.iterrows():
        # Generate the genomicLocationString
        genomicLocationString = f"{row['Chromosome']},{row['Start_Position']},{row['End_Position']},{row['Reference_Allele']},{row['Tumor_Seq_Allele2']}"
        
        # Check if the genomicLocationString is in the jsonSet
        if genomicLocationString in jsonSet:
            # If it is, increment the count in the variantsCount maps
            if row['Mutation_Status'].lower() == 'germline':
                germlineVariantsCount[genomicLocationString] = germlineVariantsCount.get(genomicLocationString, 0) + 1
            elif row['Mutation_Status'].lower() == 'somatic':
                somaticVariantsCount[genomicLocationString] = somaticVariantsCount.get(genomicLocationString, 0) + 1
            else:
                unknownVariantsCount[genomicLocationString] = unknownVariantsCount.get(genomicLocationString, 0) + 1

# Iterate over the DataFrame rows
for index, row in df.iterrows():
    # Generate the genomicLocationString
    genomicLocationString = f"{row['Chromosome']},{row['Start_Position']},{row['End_Position']},{row['Reference_Allele']},{row['Tumor_Seq_Allele2']}"
    
    # Check if the genomicLocationString is in the jsonSet
    if genomicLocationString in jsonSet:
        # If it is, increment the count in the variantsCount maps
        if row['Mutation_Status'].lower() == 'germline':
            germlineVariantsCount[genomicLocationString] = germlineVariantsCount.get(genomicLocationString, 0) + 1
        elif row['Mutation_Status'].lower() == 'somatic':
            somaticVariantsCount[genomicLocationString] = somaticVariantsCount.get(genomicLocationString, 0) + 1
        else:
            unknownVariantsCount[genomicLocationString] = unknownVariantsCount.get(genomicLocationString, 0) + 1

# Add the count number to the corresponding JSON objects
for item in data:
    for vue in item['revisedProteinEffects']:
        vue['germlineVariantsCount'] = germlineVariantsCount.get(vue['genomicLocation'], 0)
        vue['somaticVariantsCount'] = somaticVariantsCount.get(vue['genomicLocation'], 0)
        vue['unknownMutationStatusVariantsCount'] = unknownVariantsCount.get(vue['genomicLocation'], 0)

# Write the updated JSON data back to the file
with open('../VUEs.json', 'w') as f:
    json.dump(data, f, indent=4)
