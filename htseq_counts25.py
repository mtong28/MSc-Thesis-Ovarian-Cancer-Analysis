import pandas as pd

# Define the relative path to the file containing all sample names
# (Replaced local absolute path for security and reproducibility)
path = "./sample_list.txt" 

# Read sample names
file_name = pd.read_csv(path, header=None)
file_name = file_name.values.tolist()

total = []
for r in range(len(file_name)):
    name = ''.join(file_name[r])
    # Read individual count file for each sample
    each_file = pd.read_csv(name + '.csv', sep=",", header=None, names=[name])
    total.append(each_file)
    
    # Concatenate all dataframes horizontally by columns
    df = pd.concat(total, axis=1)
    df.index.name = "gene_id"

# Output the consolidated raw count matrix
df.to_csv('htseq_all_count.csv', sep=',', index=True)
