# Instuction
The `VUEs.json` file located in the `generated` folder is generated from the `VUEs.txt` file in the root directory. 

### How to generate VUEs.json
1. Update the `VUEs.txt` file in the root directory.
2. Run the `../scripts/tsv_to_json.py` script to generate the `VUEs.json` file.
```
python ../scripts/tsv_to_json.py
```
3. Run the `../scripts/variant_count.py` script to add therapeutic level and variant count to the `VUEs.json` file. Please note that the therapeutic level requires an OncoKB license. Make sure to export your OncoKB token to your local environment.
```
export ONCOKB_TOKEN=<put_your_token_here>
```
Then run the script:
```
python ../scripts/variant_count.py
```

