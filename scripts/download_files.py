import requests
import os
import requests
import json

# This is the script to download TCGA pancan data files
url = "https://api.github.com/repos/cBioPortal/datahub/git/trees/master?recursive=1"

# Target path to save the files
target_path = "./files/tcga"

# Send a GET request to the GitHub API
response = requests.get(url)
data = response.json()

# Filter the data for tree type (which represents directories)
directories = [item for item in data['tree'] if item['type'] == 'tree']

# Filter directories that contain "tcga_pan_can_atlas_2018" in their path
tcga_dirs = [dir for dir in directories if dir['path'].endswith("tcga_pan_can_atlas_2018")]

# Print the count of such directories
print(f"Found {len(tcga_dirs)} directories.")

# Go to each directory "data_mutations.txt" and "data_clinical_sample.txt"
for directory in tcga_dirs:
    # Construct the raw URL of the file
    print(directory['path'][7:])
    mutation_file_url = f"https://media.githubusercontent.com/media/cBioPortal/datahub/master/{directory['path']}/data_mutations.txt"
    clinical_file_url = f"https://media.githubusercontent.com/media/cBioPortal/datahub/master/{directory['path']}/data_clinical_sample.txt"
    
    mutation_file_response = requests.get(mutation_file_url)
    clinical_file_response = requests.get(clinical_file_url)

    if mutation_file_response.status_code == 200:
        mutation_filename = os.path.join(target_path, directory['path'].split('/')[-1] + "_mutations.txt")
        with open(mutation_filename, 'wb') as f:
            f.write(mutation_file_response.content)
        print(f"Downloaded file to {mutation_filename}")
    else:
        print(f"No 'data_mutations.txt' found in {directory['path']}")
    
    if clinical_file_response.status_code == 200:
        clinical_filename = os.path.join(target_path, directory['path'].split('/')[-1] + "_clinical_samples.txt")
        with open(clinical_filename, 'wb') as f:
            lines = clinical_file_response.content.splitlines()
            for line in lines[4:]:
                f.write(line + b'\n')

        print(f"Downloaded file to {clinical_filename}")
    else:
        print(f"No 'data_mutations.txt' found in {directory['path']}")
