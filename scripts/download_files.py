import requests
import os
import requests
import json

# This is the script to download TCGA pancan data files
url = "https://api.github.com/repos/cBioPortal/datahub/git/trees/master?recursive=1"

# Target path to save the files
target_path = "./files/"

# Send a GET request to the GitHub API
response = requests.get(url)
data = response.json()

# Filter the data for tree type (which represents directories)
directories = [item for item in data['tree'] if item['type'] == 'tree']

# Filter directories that contain "tcga_pan_can_atlas_2018" in their path
tcga_dirs = [dir for dir in directories if dir['path'].endswith("tcga_pan_can_atlas_2018")]

# Print the count of such directories
print(f"Found {len(tcga_dirs)} directories.")

# Go to each directory and download "data_mutations.txt"
for directory in tcga_dirs:
    # Construct the raw URL of the file
    print(directory['path'][7:])
    file_url = f"https://media.githubusercontent.com/media/cBioPortal/datahub/master/{directory['path']}/data_mutations.txt"
    
    # Send a GET request to the raw URL
    file_response = requests.get(file_url)

    # If the file exists, the status code of the response will be 200
    if file_response.status_code == 200:
        # Save the file with the directory's name
        filename = os.path.join(target_path, directory['path'].split('/')[-1] + ".txt")
        with open(filename, 'wb') as f:
            f.write(file_response.content)
        print(f"Downloaded file to {filename}")
    else:
        print(f"No 'data_mutations.txt' found in {directory['path']}")
