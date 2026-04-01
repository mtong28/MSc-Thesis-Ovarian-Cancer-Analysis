# ==============================================================================
# Commands Log: Upstream Quantification & Cluster Environment Setup
# Institution: Queen Mary University of London / UCL Cluster
# Author: [Your Name]
# ==============================================================================

# 1. View merged count matrix in terminal with aligned columns
cat htseq_all_count.csv | column -t -s, | less -S

# 2. Run custom python script to merge multi-sample HTSeq counts
python3 /path/to/project_dir/python_scripts/htseq_counts25.py


# 3. Environment Setup on HPC Cluster (R and Python)
# Active R-4.0.2 in cluster
export PATH=/share/apps/R-4.0.2/bin:$PATH
R

# Active Python-3.8.3
export PATH=/share/apps/python-3.8.3/bin:$PATH
export LD_LIBRARY_PATH=/share/apps/python-3.8.3/lib:$LD_LIBRARY_PATH


# 4. Manual local installation of HTSeq (to circumvent root permission constraints)
wget --no-check-certificate https://pypi.python.org/packages/source/H/HTSeq/HTSeq-0.6.1p1.tar.gz
tar -zxvf HTSeq-0.6.1p1.tar.gz
cd HTSeq-0.6.1p1/
python setup.py install --user
chmod +x scripts/htseq-count
./scripts/htseq-count


# 5. Data cleaning in R (filtering out HTSeq metadata flags)
subset(counts, ! grepl("^__",rownames(counts)))


# 6. Secure Data Transfer from Cluster (via SSH Tunnel / Jump Host)
# Step 1: Establish SSH tunnel via jump host
ssh -L 2222:cluster_node.university.edu:22 username@jump_host.university.edu

# Step 2: SCP transfer file via forwarded port to local machine
scp -P 2222 username@localhost:/path/to/project_dir/htseq_counts_results/htseq_all_count.csv ./local_folder
