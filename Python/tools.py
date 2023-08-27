# from rpy2.robjects.packages import importr
# import rpy2.robjects as ro
import csv




# function that recieves 0/1 string and switches every 0 to 1 and every 1 to 0
def switch(string):
    result = ""
    countZ = 0
    countO = 0

    for char in string:
        if char == "0":
            countZ += 1
            result += "1"
        elif char == "1":
            countO += 1
            result += "0"
        else:
            result += char
    print("Number of 0's: " + str(countZ))
    print("Number of 1's: " + str(countO))
    return result

class Study:
    def __init__(self, path, country, tissue, n):
        self.path = path
        self.country = country
        self.tissue = tissue
        self.n = n
        temp = path.split('_')[0]
        self.name = temp.split('/')[-1]


def read_top_config_file(filename):
    paths = []

    with open(filename, 'r') as file:
        for line in file:
            paths.append()
            line = line.strip()



def read_config_file(filename):
    studies = []

    with open(filename, 'r') as file:
        for line in file:
            line = line.strip()
            if line:
                elements = line.split(',')

                if len(elements) == 4:
                    path, country, tissue, n = elements
                    study = Study(path, country, tissue, n)
                    studies.append(study)
                else:
                    print(f"Ignoring invalid line: {line}")

    return studies

def print_studies(studies):
    for study in studies:
        print(f"Name: {study.name}")
        print(f"Country: {study.country}")
        print(f"Tissue: {study.tissue}")
        print(f"N: {study.n}")
        print(f"Group String: {study.group_string}")
        print(f"File Format: {study.format}")
        print()

def count_tsv_rows(file_path):
    with open(file_path, 'r') as tsv_file:
        row_count = sum(1 for row in tsv_file)
        print(f"The file '{file_path}' has {row_count} rows.")

def check_gene_logfc(gene, tsv_path):
    with open(tsv_path, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        next(reader)  # Skip the header row
        for row in reader:
            logfc = float(row[5])
            gene_name = row[6].strip('"')
            if gene_name == gene:
                if logfc < -0.00005:
                    return 1
                elif logfc > 0.00005:
                    return -1
                else:
                    return 0
    return f"{gene} not found in the TSV file."

def check_gene_pval(gene, tsv_path):
    with open(tsv_path, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        next(reader)  # Skip the header row
        for row in reader:
            pval = float(row[1])
            gene_name = row[6].strip('"')
            if gene_name == gene:
                return pval
    return f"{gene} not found in the TSV file."

def get_logFC_counts(studies, gene):
    up_counts = 0
    down_counts = 0
    nut_counts = 0
    for study in studies:
        count = check_gene_logfc(gene, study.path)
        if count == 1:
            up_counts += 1
        if count == -1:
            down_counts += 1
        if count == 0:
            nut_counts += 1
        
    return (up_counts, down_counts, nut_counts)

def get_correlation_studies(gene1, gene2, factor, studies):
    corelated = []
    for study in studies:
        count = check_gene_logfc(gene2, study.path)
        if count*factor == 1:
            if check_gene_pval(gene1, study.path) < 0.05:
                corelated.append(study)
    return corelated

print_fc_counts = lambda gene, counts: print(f"{gene} up: {counts[0]}, down: {counts[1]}, low: {counts[2]}")
print_gene = lambda gene: print(f"Gene: {gene[0]}, count: {gene[1]}")

def print_correlation_studies(gene1, studies):
    print(f"{gene1} qualifies as indicator (P.val < 0.05, and an expected FC), in: {len(studies)} datasets:")
    for study in studies:
        print(f"Study: {study.name}: {study.country}, {study.tissue}, {study.n} samples.")
        fc = check_gene_logfc("TMED9", study.path)
        pv = check_gene_pval("TMED9", study.path)
        if fc == 1:
            print(f"In this dataset TMED9 exhibits an UP regulation, with a P value of: {pv}")
        elif fc == -1:
            print(f"In this dataset TMED9 exhibits a DOWN regulation, with a P value of: {pv}")
        else:
            print(f"In this dataset TMED9 exhibits a NO change, with a P value of: {pv}")
        print()

def get_report(studies):
    # get tmed Fc info
    tmed_fc_counts = get_logFC_counts(studies, "TMED9")
    print_fc_counts("TMED9", tmed_fc_counts)
    print()
    # get NFKB1 Fc info
    nfkb1_fc_counts = get_logFC_counts(studies, "NFKB1")
    print_fc_counts("NFKB1", nfkb1_fc_counts)
    print()
    # get PDE4B Fc info
    pde4b_fc_counts = get_logFC_counts(studies, "PDE4B")
    print_fc_counts("PDE4B", pde4b_fc_counts)
    print()
    # get SSR Fc info
    ssr_fc_counts = get_logFC_counts(studies, "SRR")
    print_fc_counts("SRR", ssr_fc_counts)
    print()
    # get GAPDH Fc info
    gapdh_fc_counts = get_logFC_counts(studies, "GAPDH")
    print_fc_counts("GAPDH", gapdh_fc_counts)
    print()

    # get correlation info
    # TMED9 and NFKB1
    tmed_nfkb1 = get_correlation_studies("TMED9", "NFKB1", 1, studies)
    print_correlation_studies("NFKB1", tmed_nfkb1)
    print()
    # TMED9 and PDE4B
    tmed_pde4b = get_correlation_studies("TMED9", "PDE4B", 1, studies)
    print_correlation_studies("PDE4B", tmed_pde4b)

def get_freq_genes(path):
    genes = []
    with open(path, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        next(reader)  # Skip the header row
        for row in reader:
            gene_name = row[6].strip('"')
            pval = float(row[1])
            fc = float(row[5])
            if pval < 0.1 and abs(fc) > 0.25:
                if gene_name not in genes:
                    genes.append((gene_name, 0))
                else:
                    for gene in genes:
                        if gene[0] == gene_name:
                            gene[1] += 1
    return genes




print("Reading config file")
studies = read_config_file("Configs_report.txt")
print("Finished reading config file")
get_report(studies)

# genes = get_freq_genes('top1000_combined_table.tsv')
# print("Overall top frequent genes:")
# for gene in genes:
#     print_gene(gene)

# # print_studies(studies)
# print(f'Number of studies: {len(studies)}')
# # Load the R script using the source function
# r = ro.r
# r.source('EDA.R')
# print("performing analysis")
# for study in studies:
#     print(f"Running analysis for {study.name}")
#     print(f"Country: {study.country}")
#     print(f"Tissue: {study.tissue}")
#     print(f"N: {study.n}")
#     name = study.name
#     format = study.format
#     group_string = study.group_string
#     print()
#     # Call a function from the loaded script and pass arguments
#     run_study_analysis = ro.globalenv['run_study_analysis']
#     run_study_analysis(name, group_string, format)
#     print("Finished running analysis for " + study.name)
# print("Finished running all analysis")



