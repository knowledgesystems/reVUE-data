import pandas as pd
import json
import glob

# Get a list of all TXT files in the mskimpact and tcga folders
mskimpact_files = glob.glob('./files/mskimpact/*.txt')
tcga_files = glob.glob('./files/tcga/*.txt')

# Create maps for variants count
mskimpactCounts = {'germlineVariantsCount': {}, 'somaticVariantsCount': {}, 'unknownVariantsCount': {}, 'totalSampleCount': 0, 'sampleCountByGene': {}}
tcgaCounts = {'germlineVariantsCount': {}, 'somaticVariantsCount': {}, 'unknownVariantsCount': {}, 'totalSampleCount': 0, 'sampleCountByGene': {}}

# Read the JSON file
with open('../VUEs.json', 'r') as f:
    data = json.load(f)

# Create a set from the JSON data
jsonSet = {vue['genomicLocation'] for vueSet in data for vue in vueSet['revisedProteinEffects']}

# Define a function to update the counts
def updateCounts(files, counts):
    # Iterate over each file
    for file in files:
        # Read each TXT file into a DataFrame
        df = pd.read_csv(file, sep="\t", usecols=['Hugo_Symbol', 'Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2', 'Mutation_Status', 'Tumor_Sample_Barcode'])
        
        print("Counting from " + file)
        # Update the total sample count
        counts['totalSampleCount'] += df['Tumor_Sample_Barcode'].nunique()

        # Group the DataFrame by 'Hugo_Symbol' and count unique 'Tumor_Sample_Barcode' values
        sampleCountByGene = df.groupby('Hugo_Symbol')['Tumor_Sample_Barcode'].nunique()
        
        # Update the sampleCountByGene map
        for gene, count in sampleCountByGene.iteritems():
            counts['sampleCountByGene'][gene] = counts['sampleCountByGene'].get(gene, 0) + count

        # Iterate over the DataFrame rows
        for index, row in df.iterrows():
            # Generate the genomicLocationString
            genomicLocationString = f"{row['Chromosome']},{row['Start_Position']},{row['End_Position']},{row['Reference_Allele']},{row['Tumor_Seq_Allele2']}"
            
            # Check if the genomicLocationString is in the jsonSet
            if genomicLocationString in jsonSet:
                # If it is, increment the count in the corresponding variantsCount map
                if row['Mutation_Status'].lower() == 'germline':
                    counts['germlineVariantsCount'][genomicLocationString] = counts['germlineVariantsCount'].get(genomicLocationString, 0) + 1
                elif row['Mutation_Status'].lower() == 'somatic':
                    counts['somaticVariantsCount'][genomicLocationString] = counts['somaticVariantsCount'].get(genomicLocationString, 0) + 1
                else:
                    counts['unknownVariantsCount'][genomicLocationString] = counts['unknownVariantsCount'].get(genomicLocationString, 0) + 1

# Update the counts for the mskimpact and tcga files
updateCounts(mskimpact_files, mskimpactCounts)
updateCounts(tcga_files, tcgaCounts)

# Add the count numbers to the corresponding JSON objects
for item in data:
    for vue in item['revisedProteinEffects']:
        vue['counts'] = {
            'mskimpact': {
                "germlineVariantsCount": mskimpactCounts["germlineVariantsCount"].get(vue['genomicLocation'], 0),
                "somaticVariantsCount": mskimpactCounts["somaticVariantsCount"].get(vue['genomicLocation'], 0),
                "unknownVariantsCount": mskimpactCounts["unknownVariantsCount"].get(vue['genomicLocation'], 0),
                "totalSampleCount": mskimpactCounts.get("totalSampleCount"),
                "geneSampleCount": mskimpactCounts['sampleCountByGene'].get(item['hugoGeneSymbol'], 0)
            },
            'tcga': {
                "germlineVariantsCount": tcgaCounts["germlineVariantsCount"].get(vue['genomicLocation'], 0),
                "somaticVariantsCount": tcgaCounts["somaticVariantsCount"].get(vue['genomicLocation'], 0),
                "unknownVariantsCount": tcgaCounts["unknownVariantsCount"].get(vue['genomicLocation'], 0),
                "totalSampleCount": tcgaCounts.get("totalSampleCount"),
                "geneSampleCount": tcgaCounts['sampleCountByGene'].get(item['hugoGeneSymbol'], 0)
            }
        }

# Write the updated JSON data back to the file
with open('../VUEs.json', 'w') as f:
    json.dump(data, f, indent=4)
