import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np


def filter_genes(input_file, output_file, gene_list):
    # Read the input TSV file into a DataFrame
    df = pd.read_csv(input_file, sep='\t')

    # Filter the rows based on the 'Gene.symbol' column
    filtered_df = df[df['Gene.symbol'].isin(gene_list)]

    # Write the filtered DataFrame to the output TSV file
    filtered_df.to_csv(output_file, sep='\t', index=False)

def read_dir_filter(dir_path, choice):
    for file_name in os.listdir(dir_path):
        if file_name.endswith('logFCNorm.tsv'):
            file_path = os.path.join(dir_path, file_name)
            if choice == 1:
                genes = ["TMED9", "NFKB1", "SRR", "GAPDH", "PDE4B"]
                # create name for output file
                name = os.path.splitext(file_name)[0].split("_")[0] + "_n_control_genes.tsv"
            if choice == 3:
                genes = ["SFRP4", "SLC2A2", "IGF1R", "TMED9"]
                name  = os.path.splitext(file_name)[0].split("_")[0] + "_significant_genes.tsv"
            if choice == 4:
                # genes = list
                name = os.path.splitext(file_name)[0].split("_")[0] + "_50_top_sig_genes.tsv"
            # else:
            #     genes = ["TMED1", "TMED2", "TMED3", "TMED4", "TMED5", "TMED6", "TMED7", "TMED9", "TMED10"]
            #     name  = os.path.splitext(file_name)[0].split("_")[0] + "_n_TMED_genes.tsv"

            output_file = os.path.join(dir_path, name)
            filter_genes(file_path, output_file, genes)
    

def min_max_normalize(input_dir, identifier="_logFCNorm"):
    # Get a list of all TSV files in the input directory
    tsv_files = [file for file in os.listdir(input_dir) if file.endswith('output_table.tsv')]

    for file in tsv_files:
        # Read the TSV file into a DataFrame
        file_path = os.path.join(input_dir, file)
        df = pd.read_csv(file_path, sep='\t')

        # Step 1: Find the current minimum and maximum values in the 'logFC' column
        current_min = df['logFC'].min()
        current_max = df['logFC'].max()

        # Step 2: Define the new range you want the values to be normalized to
        new_min = -4
        new_max = 4

        # Step 3: Perform Min-Max normalization
        df['normalized_logFC'] = ((df['logFC'] - current_min) / (current_max - current_min)) * (new_max - new_min) + new_min

        # Construct the output file name
        output_file_name = os.path.splitext(file)[0] + identifier + ".tsv"
        output_file_path = os.path.join(input_dir, output_file_name)

        # Write the DataFrame with normalized values to the output file
        df.to_csv(output_file_path, sep='\t', index=False)

def extract_gene_name(row):
    gene_symbol = row['Gene.symbol']
    # if gene_symbol.startswith('TMED'):
    #     return gene_symbol[4:]  # Extract the number part from TMED#
    # else:
    return gene_symbol   # Take the first 3 letter for other gene symbols


def plot_scatter_colored_by_gene(data):
    data['name'] = data.apply(extract_gene_name, axis=1)
    
    # Set up the figure and axis
    plt.figure(figsize=(12, 10))
    ax = plt.gca()

    # Define the color map for genes
    cmap = plt.get_cmap('tab10')
    gene_colors = {gene: cmap(i % 10) for i, gene in enumerate(data['name'].unique())}

    # Plot the scatter points
    for gene, color in gene_colors.items():
        gene_data = data[data['name'] == gene]
        plt.scatter(gene_data['normalized_logFC'], gene_data['P.Value'], color=color, label=gene, alpha=0.7)
        # Add gene names as labels above the dots
        for i in range(len(gene_data)):
            plt.text(gene_data.iloc[i]['normalized_logFC'], gene_data.iloc[i]['P.Value'],
                     gene_data.iloc[i]['name'], fontsize=8, ha='center', va='bottom')

    # Add labels and title
    plt.xlabel('LogFC')
    plt.ylabel('P-Value')
    plt.title('Scatter Plot: P-Value vs. LogFC (Control Genes)')

    # Add legend for genes outside the plot
    plt.legend(title='Gene', loc='center left', bbox_to_anchor=(1, 0.5))

    # Show grid and set plot limits
    plt.grid(True)
    plt.xlim(data['normalized_logFC'].min() - 1, data['normalized_logFC'].max() + 1)
    plt.ylim(data['P.Value'].min() - 0.1, data['P.Value'].max() + 0.1)

    # Adjust layout to make space for the legend
    plt.tight_layout()

    # Show the plot
    plt.show()

def plot_scatter_with_labels(data, choice):
    data['name'] = data.apply(extract_gene_name, axis=1)
    # Set up the figure and axis
    plt.figure(figsize=(10, 8))
    ax = plt.gca()

    # Define the color map for studies
    cmap = plt.get_cmap('tab20')
    study_colors = {study: cmap(i) for i, study in enumerate(data['study'].unique())}

    # Plot the scatter points
    for study, color in study_colors.items():
        study_data = data[data['study'] == study]
        plt.scatter(study_data['normalized_logFC'], study_data['P.Value'], color=color, label=study, alpha=0.7)
        # Add gene names as labels above the dots
        for i in range(len(study_data)):
            plt.text(study_data.iloc[i]['normalized_logFC'], study_data.iloc[i]['P.Value'],
                     study_data.iloc[i]['name'], fontsize=8, ha='center', va='bottom')

    # Add labels and title
    plt.xlabel('LogFC')
    plt.ylabel('P-Value')
    if choice == 1:
        plt.title('Scatter Plot: P-Value vs. LogFC (TMED Genes)')
    elif choice == 3:
        plt.title('Scatter Plot: P-Value vs. LogFC (Significant Genes)')
    else:
        plt.title('Scatter Plot: P-Value vs. LogFC (Control Genes)')
    plt.legend(title='Study', loc='upper right')

    # Show grid and set plot limits
    plt.grid(True)
    plt.xlim(data['normalized_logFC'].min() - 1, data['normalized_logFC'].max() + 1)
    plt.ylim(data['P.Value'].min() - 0.1, data['P.Value'].max() + 0.1)

    # Show the plot
    plt.show()

def plot_scatter_colored_by_gene_labled_tissue(data):
    data['name'] = data.apply(extract_gene_name, axis=1)
    data['tissue'] = data['study'].map(study_to_tissue)  # map study to tissue
    
    # Set up the figure and axis
    plt.figure(figsize=(12, 10))
    ax = plt.gca()

    # Define the color map for genes
    cmap = plt.get_cmap('tab10')
    gene_colors = {gene: cmap(i % 10) for i, gene in enumerate(data['name'].unique())}

    # Plot the scatter points
    for gene, color in gene_colors.items():
        gene_data = data[data['name'] == gene]
        plt.scatter(gene_data['normalized_logFC'], gene_data['P.Value'], color=color, label=gene, alpha=0.7)
        # Add tissue names as labels above the dots
        for i in range(len(gene_data)):
            plt.text(gene_data.iloc[i]['normalized_logFC'], gene_data.iloc[i]['P.Value'],
                     gene_data.iloc[i]['tissue'], fontsize=8, ha='center', va='bottom')

    # Add labels and title
    plt.xlabel('LogFC')
    plt.ylabel('P-Value')
    plt.title('Scatter Plot: P-Value vs. LogFC (Control Genes)')

    # Add legend for genes outside the plot
    plt.legend(title='Gene', loc='center left', bbox_to_anchor=(1, 0.5))

    # Show grid and set plot limits
    plt.grid(True)
    plt.xlim(data['normalized_logFC'].min() - 1, data['normalized_logFC'].max() + 1)
    plt.ylim(data['P.Value'].min() - 0.1, data['P.Value'].max() + 0.1)

    # Adjust layout to make space for the legend
    plt.tight_layout()

    # Show the plot
    plt.show()



def read_and_combine_tsv_files_from_directory(dir_path):
    data_structure = []  # List to store data from all TSV files

    # Loop through all files in the directory
    for file_name in os.listdir(dir_path):
        # Check for both TMED and control genes TSV files
        if file_name.endswith("n_TMED_genes.tsv") or file_name.endswith("n_control_genes.tsv"):
            print("Reading file: {}".format(file_name))
            file_path = os.path.join(dir_path, file_name)

            # Read the TSV file into a Pandas DataFrame
            df = pd.read_csv(file_path, delimiter="\t")

            # Add a new column 'study' containing the study ID (file name without extension)
            df['study'] = os.path.splitext(file_name)[0].split("_")[0]

            no_dup_df = df.drop_duplicates(subset=['Gene.symbol'], keep='first')

            # Append the DataFrame to the data_structure list
            data_structure.append(no_dup_df)

    # Combine all DataFrames into a single DataFrame
    combined_df = pd.concat(data_structure, ignore_index=True)

    return combined_df


def read_tsv_files_from_directory(dir_path, choice):
    data_structure = []  # List to store data from all TSV files

    # Loop through all files in the directory
    for file_name in os.listdir(dir_path):
        if choice == 1:
            if file_name.endswith("n_TMED_genes.tsv"):
                print("Reading file: {}".format(file_name))
                file_path = os.path.join(dir_path, file_name)

                # Read the TSV file into a Pandas DataFrame
                df = pd.read_csv(file_path, delimiter="\t")

                # Add a new column 'study' containing the study ID (file name without extension)
                df['study'] = os.path.splitext(file_name)[0].split("_")[0]

                no_dup_df = df.drop_duplicates(subset=['Gene.symbol'], keep='first')

                # Append the DataFrame to the data_structure list
                data_structure.append(no_dup_df)
        if choice == 3:
            if file_name.endswith("_significant_genes.tsv"):
                print("Reading file: {}".format(file_name))
                file_path = os.path.join(dir_path, file_name)

                # Read the TSV file into a Pandas DataFrame
                df = pd.read_csv(file_path, delimiter="\t")

                # Add a new column 'study' containing the study ID (file name without extension)
                df['study'] = os.path.splitext(file_name)[0].split("_")[0]

                no_dup_df = df.drop_duplicates(subset=['Gene.symbol'], keep='first')

                # Append the DataFrame to the data_structure list
                data_structure.append(no_dup_df)

        else:
            if file_name.endswith("n_control_genes.tsv"):
                print("Reading file: {}".format(file_name))
                file_path = os.path.join(dir_path, file_name)

                # Read the TSV file into a Pandas DataFrame
                df = pd.read_csv(file_path, delimiter="\t")

                # Add a new column 'study' containing the study ID (file name without extension)
                df['study'] = os.path.splitext(file_name)[0].split("_")[0]

                no_dup_df = df.drop_duplicates(subset=['Gene.symbol'], keep='first')

                # Append the DataFrame to the data_structure list
                data_structure.append(no_dup_df)

    # Combine all DataFrames into a single DataFrame
    combined_df = pd.concat(data_structure, ignore_index=True)

    return combined_df
def read_results(path):
    genes = []
    with open(path, 'r') as f:
        lines = f.readlines()
        for line in lines:
            l = line.strip()
            gene = l.split('.')[1]
            g = gene.strip()
            genes.append(g)
    return genes[:50]

def select_top_genes(data, top_n=100):
    # Create a dictionary to store top genes for each study
    study_genes_dict = {}
    
    # Group the data by the 'study' column
    grouped_data = data.groupby('study')
    
    # Iterate over each study group and select the top genes
    for study, group in grouped_data:
        # Sort the group by 'normalized_logFC' in descending order and select the top n genes
        top_genes = group.nlargest(top_n, 'B')
        study_genes_dict[study] = set(top_genes['Gene.symbol'])
    
    # Find the intersection of genes across all studies
    common_genes = set.intersection(*study_genes_dict.values())
    
    # Filter the data to include only the common genes
    common_genes_data = data[data['Gene.symbol'].isin(common_genes)]
    
    return common_genes_data


def plot_heatmap_sort_country(data):
    # Group the data and form the heatmap structure
    heatmap_df = data.groupby(['Gene.symbol', 'study'])['normalized_logFC'].mean().unstack()
    
    # Sort rows based on the 'TMED9' values
    if 'TMED9' in heatmap_df.index:
        tmed9_sorted_order = heatmap_df.loc['TMED9'].sort_values(ascending=False).index
        heatmap_df = heatmap_df[tmed9_sorted_order]

    plt.figure(figsize=(12, 14))  # Adjust the figure size as per your preference
    
    # Create a grid to place two subplots (1 row, 1 column)
    grid = plt.GridSpec(2, 1, height_ratios=[1, 15], hspace=0.5)
    
    # Top axis for the title
    ax_top = plt.subplot(grid[0])
    ax_top.text(0.5, 0.5, "Heatmap of TMED vs. Studies (logFC) Country", ha='center', va='center', fontsize=12, transform=ax_top.transAxes)
    ax_top.axis('off')  # Hide all axis markings
    
    # Main heatmap axis
    ax_main = plt.subplot(grid[1])
    sns.heatmap(heatmap_df, ax=ax_main, cmap='coolwarm', annot=False, fmt=".2f", linewidths=0.5, cbar=False)
    plt.yticks(rotation=45)
    
    # Set study labels on x-axis ticks
    ax_main.set_xticks([x + 0.5 for x in range(len(heatmap_df.columns))])  # Adjust the position to center of each column
    ax_main.set_xticklabels(heatmap_df.columns, rotation=45, ha='center')

    # Use twinx() to create a second x-axis for tissue labels
    ax_countries = ax_main.twiny()
    ax_countries.set_xlim(ax_main.get_xlim())
    ax_countries.set_xticks(ax_main.get_xticks())
    ax_countries.set_xticklabels([study_to_country.get(study, study) for study in heatmap_df.columns], rotation=45, ha='center')
    ax_countries.set_xlabel("Country")
    
    ax_main.set_xlabel("Study")
    ax_main.set_ylabel("Gene Symbol")
    
    plt.show()


def plot_heatmap_sort_tissue(data):
    # Group the data and form the heatmap structure
    heatmap_df = data.groupby(['Gene.symbol', 'study'])['normalized_logFC'].mean().unstack()
    
    # Sort rows based on the 'TMED9' values
    if 'TMED9' in heatmap_df.index:
        tmed9_sorted_order = heatmap_df.loc['TMED9'].sort_values(ascending=False).index
        heatmap_df = heatmap_df[tmed9_sorted_order]

    plt.figure(figsize=(12, 14))  # Adjust the figure size as per your preference
    
    # Create a grid to place two subplots (1 row, 1 column)
    grid = plt.GridSpec(2, 1, height_ratios=[1, 15], hspace=0.5)
    
    # Top axis for the title
    ax_top = plt.subplot(grid[0])
    ax_top.text(0.5, 0.5, "Heatmap of TMED vs. Studies (logFC) Tissue", ha='center', va='center', fontsize=12, transform=ax_top.transAxes)
    ax_top.axis('off')  # Hide all axis markings
    
    # Main heatmap axis
    ax_main = plt.subplot(grid[1])
    sns.heatmap(heatmap_df, ax=ax_main, cmap='coolwarm', annot=False, fmt=".2f", linewidths=0.5, cbar=False)
    plt.yticks(rotation=45)
    
    # Set study labels on x-axis ticks
    ax_main.set_xticks([x + 0.5 for x in range(len(heatmap_df.columns))])  # Adjust the position to center of each column
    ax_main.set_xticklabels(heatmap_df.columns, rotation=45, ha='center')

    # Use twinx() to create a second x-axis for tissue labels
    ax_tissues = ax_main.twiny()
    ax_tissues.set_xlim(ax_main.get_xlim())
    ax_tissues.set_xticks(ax_main.get_xticks())
    ax_tissues.set_xticklabels([study_to_tissue.get(study, study) for study in heatmap_df.columns], rotation=45, ha='center')
    ax_tissues.set_xlabel("Tissue")
    
    ax_main.set_xlabel("Study")
    ax_main.set_ylabel("Gene Symbol")
    
    plt.show()



def plot_heatmap_sort_tissue_and_country_order(data):
    # Group the data and form the heatmap structure
    heatmap_df = data.groupby(['Gene.symbol', 'study'])['normalized_logFC'].mean().unstack()

    # Sort rows based on the 'TMED9' values
    if 'TMED9' in heatmap_df.index:
        tmed9_sorted_order = heatmap_df.loc['TMED9'].sort_values(ascending=False).index
        heatmap_df = heatmap_df[tmed9_sorted_order]

    # Create a mapping from study to "tissue-country"
    study_to_country_tissue = {study: f"{study_to_country[study]} - {study_to_tissue[study]}" for study in heatmap_df.columns}
    
    # Separate genes into control genes, TMED9, and other TMEDs
    control_genes = heatmap_df[heatmap_df.index.str.startswith('TMED') == False]
    tmed9 = heatmap_df[heatmap_df.index == 'TMED9']
    other_tmeds = heatmap_df[(heatmap_df.index.str.startswith('TMED')) & (heatmap_df.index != 'TMED9')]
    
    # Concatenate in the desired order
    ordered_heatmap_df = pd.concat([tmed9, control_genes, other_tmeds])

    plt.figure(figsize=(12, 14))  # Adjust the figure size as per your preference
    
    # Create a grid to place two subplots (1 row, 1 column)
    grid = plt.GridSpec(2, 1, height_ratios=[1, 15], hspace=0.5)
    
    # Top axis for the title
    ax_top = plt.subplot(grid[0])
    ax_top.text(0.5, 0.5, "Heatmap of TMED vs. Studies (logFC)", ha='center', va='center', fontsize=12, transform=ax_top.transAxes)
    ax_top.axis('off')  # Hide all axis markings
    
    # Main heatmap axis
    ax_main = plt.subplot(grid[1])
    sns.heatmap(ordered_heatmap_df, ax=ax_main, cmap='coolwarm', annot=False, fmt=".2f", linewidths=0.5, cbar=False)
    plt.yticks(rotation=45)

    # Set study labels on x-axis ticks
    ax_main.set_xticks([x + 0.5 for x in range(len(ordered_heatmap_df.columns))])  # Adjust the position to center of each column
    ax_main.set_xticklabels(ordered_heatmap_df.columns, rotation=45, ha='center')

    # Use twinx() to create a second x-axis for tissue labels
    ax_countries = ax_main.twiny()
    ax_countries.set_xlim(ax_main.get_xlim())
    ax_countries.set_xticks(ax_main.get_xticks())
    ax_countries.set_xticklabels([study_to_country_tissue.get(study, study) for study in heatmap_df.columns], rotation=45, ha='left')
    ax_countries.set_xlabel("Coutnry - Tissue")
    
    ax_main.set_xlabel("Study")
    ax_main.set_ylabel("Gene Symbol")
    
    plt.show()

def plot_heatmap_sort_tissue_order(data):
    # Group the data and form the heatmap structure
    heatmap_df = data.groupby(['Gene.symbol', 'study'])['normalized_logFC'].mean().unstack()

    # Sort rows based on the 'TMED9' values
    if 'TMED9' in heatmap_df.index:
        tmed9_sorted_order = heatmap_df.loc['TMED9'].sort_values(ascending=False).index
        heatmap_df = heatmap_df[tmed9_sorted_order]
    
    # Separate genes into control genes, TMED9, and other TMEDs
    control_genes = heatmap_df[heatmap_df.index.str.startswith('TMED') == False]
    tmed9 = heatmap_df[heatmap_df.index == 'TMED9']
    other_tmeds = heatmap_df[(heatmap_df.index.str.startswith('TMED')) & (heatmap_df.index != 'TMED9')]
    
    # Concatenate in the desired order
    ordered_heatmap_df = pd.concat([tmed9, control_genes, other_tmeds])

    plt.figure(figsize=(12, 14))  # Adjust the figure size as per your preference
    
    # Create a grid to place two subplots (1 row, 1 column)
    grid = plt.GridSpec(2, 1, height_ratios=[1, 15], hspace=0.5)
    
    # Top axis for the title
    ax_top = plt.subplot(grid[0])
    ax_top.text(0.5, 0.5, "Heatmap of TMED vs. Studies (logFC) Tissue", ha='center', va='center', fontsize=12, transform=ax_top.transAxes)
    ax_top.axis('off')  # Hide all axis markings
    
    # Main heatmap axis
    ax_main = plt.subplot(grid[1])
    sns.heatmap(ordered_heatmap_df, ax=ax_main, cmap='coolwarm', annot=False, fmt=".2f", linewidths=0.5, cbar=False)
    plt.yticks(rotation=45)

    # Set study labels on x-axis ticks
    ax_main.set_xticks([x + 0.5 for x in range(len(ordered_heatmap_df.columns))])  # Adjust the position to center of each column
    ax_main.set_xticklabels(ordered_heatmap_df.columns, rotation=45, ha='center')

    # Use twinx() to create a second x-axis for tissue labels
    ax_tissues = ax_main.twiny()
    ax_tissues.set_xlim(ax_main.get_xlim())
    ax_tissues.set_xticks(ax_main.get_xticks())
    ax_tissues.set_xticklabels([study_to_tissue.get(study, study) for study in ordered_heatmap_df.columns], rotation=45, ha='center')
    ax_tissues.set_xlabel("Tissue")
    
    ax_main.set_xlabel("Study")
    ax_main.set_ylabel("Gene Symbol")
    
    plt.show()

# Call the function with your data
# plot_heatmap_sort(your_data)

def plot_heatmap_sort(data):
    # Group the data and form the heatmap structure
    heatmap_df = data.groupby(['Gene.symbol', 'study'])['normalized_logFC'].mean().unstack()
    
    # Sort rows based on the 'TMED9' values
    if 'TMED9' in heatmap_df.index:
        tmed9_sorted_order = heatmap_df.loc['TMED9'].sort_values(ascending=False).index
        heatmap_df = heatmap_df[tmed9_sorted_order]

    plt.figure(figsize=(12, 14))  # Adjust the figure size as per your preference
    sns.heatmap(heatmap_df, cmap='coolwarm', annot=False, fmt=".2f", linewidths=0.5, cbar=False)
    plt.yticks(rotation=45)
    plt.title("Heatmap of TMED vs. Studies (logFC)")
    plt.xlabel("Study")
    plt.ylabel("Gene Symbol")
    plt.show()


def plot_heatmap(data):
    heatmap_df = data.groupby(['Gene.symbol', 'study'])['normalized_logFC'].mean().unstack()
    plt.figure(figsize=(12, 14))  # Adjust the figure size as per your preference
    sns.heatmap(heatmap_df, cmap='coolwarm', annot=False, fmt=".2f", linewidths=0.5, cbar=False)
    plt.yticks(rotation=45)
    plt.title("Heatmap of TMED vs. Studies (logFC)")
    plt.xlabel("Study")
    plt.ylabel("Gene Symbol")
    plt.show()

def parse_config(file_name):
    study_to_tissue = {}
    study_to_country = {}

    with open(file_name, 'r') as f:
        for line in f.readlines():
            print(line)
            fields = line.strip().split(',')
            study_id = fields[0]
            country = fields[1]
            tissue = fields[2]

            study_to_tissue[study_id] = tissue
            study_to_country[study_id] = country
            if study_id == "GSE15773":
                break

    return study_to_tissue, study_to_country

# Call the function to create the mapping
study_to_tissue, study_to_country = parse_config("Configs.txt")



print("hello")
input_directory = "BeforeNormDataSetsOutputs"
# gene_list = read_results('top100_results.txt')
# read_dir_filter(input_directory, choice=5)
# # Example usage:
# directory_path = "BeforeNormDataSetsOutputs"
combined_data = read_and_combine_tsv_files_from_directory(input_directory)
print(combined_data.head(15))
# plot_heatmap_sort_tissue_order(combined_data)
plot_heatmap_sort_tissue_and_country_order(combined_data)

# top = select_top_genes(combined_data, top_n=100)
# print(top.head(150))
# # Call the function to generate the scatter plot
# plot_scatter_colored_by_gene(combined_data)
# plot_scatter_with_labels(combined_data, 5)
# plot_heatmap(combined_data)
