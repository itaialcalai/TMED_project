import numpy as np
from scipy.stats import chi2
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.mixture import GaussianMixture
from scipy.stats import pearsonr
import statsmodels.api as sm
from typing import List, Tuple, Dict

def fisher_method(p_values):
    # Calculate the test statistic Fg for each gene g
    Fg = -2 * np.sum(np.log(p_values))
    # Calculate the combined p-value using the χ2 distribution
    combined_p_value = 1 - chi2.cdf(Fg, df=2 * len(p_values))
    return combined_p_value

class Study:
    def __init__(self, name, path, country, tissue, n):
        self.name = name
        self.path = path
        self.country = country
        self.tissue = tissue
        self.n = int(n)  # Convert n to an integer

    # print method for debugging
    def __repr__(self):
        return f"Study(name={self.name}, path={self.path}, country={self.country}, tissue={self.tissue}, n={self.n})"

def read_tsv(file_path):
    df = pd.read_csv(file_path, sep='\t')
    return df[['normalized_logFC', 'P.Value', 'Gene.symbol']]

def perform_meta_analysis(studies):
    values_dict = {}
    logfc_values_all = []  # collect all logFC values
    for study in studies:
        df = read_tsv(study.path)
        df = df.drop_duplicates(subset=['Gene.symbol'], keep='first')  # drop duplicate genes
        df = df.dropna(subset=['normalized_logFC', 'P.Value'])  # Drop rows with NaN in logFC or P.Value column
        logfc_values = df['normalized_logFC'].values.astype(float)
        p_values = df['P.Value'].values.astype(float)
        gene_symbols = df['Gene.symbol'].values
        logfc_values_all.extend(logfc_values)

        for symbol, logfc_value, p_value in zip(gene_symbols, logfc_values, p_values):
            if symbol not in values_dict:
                values_dict[symbol] = {'logfc': [], 'pval': []}
            values_dict[symbol]['logfc'].append(logfc_value)
            values_dict[symbol]['pval'].append(p_value)  # no longer weight p-value

    # Standardize logFC values
    logfc_values_all = np.array(logfc_values_all).reshape(-1, 1)
    logfc_values_all = (logfc_values_all - np.mean(logfc_values_all)) / np.std(logfc_values_all)

    # Fit GMM
    gmm = GaussianMixture(n_components=2)
    gmm.fit(logfc_values_all)

    combined_scores = []
    for symbol, values in values_dict.items():
        logfc_values = np.array(values['logfc']).reshape(-1, 1)
        logfc_values = (logfc_values - np.mean(logfc_values_all)) / np.std(logfc_values_all)  # standardize
        weights = gmm.predict_proba(logfc_values)[:, 1]  # take the weights for the 2nd Gaussian
        weighted_p_values = np.array(values['pval']) * weights  # weight p-values
        combined_p_value = fisher_method(weighted_p_values)
        avg_abs_logfc = np.abs(np.mean(np.array(values['logfc'])))
        combined_p_value = max(combined_p_value, 1e-300)
        # print(f"{symbol}: {avg_abs_logfc}, {combined_p_value}")
        score = avg_abs_logfc / combined_p_value  # Calculate the composite score
        combined_scores.append(score)

    combined_scores = np.array(combined_scores)

    return values_dict, combined_scores

def get_top_500_genes(combined_scores):
    # Sort the genes based on the combined scores
    sorted_genes_indices = np.argsort(combined_scores)[::-1]  # [::-1] to sort in descending order
    top_500_genes_indices = sorted_genes_indices[:500]
    return top_500_genes_indices

def read_config_file(filename):
    studies = []

    with open(filename, 'r') as file:
        for line in file:
            line = line.strip()
            if line != "---":
                elements = line.split(',')

                if len(elements) >= 4:
                    name, country, tissue, n = elements[:4]
                    path = f"BeforeNormDataSetsOutputs/{name}_output_table_logFCNorm.tsv"
                    study = Study(name, path, country, tissue, n)
                    studies.append(study)
                else:
                    print(f"Ignoring invalid line: {line}")
            else:
                break

    return studies

def compile_data(studies: List[Study]) -> pd.DataFrame:
    all_data = []
    for study in studies:
        data = pd.read_csv(study.path, sep='\t')
        data = data.groupby('Gene.symbol').mean(numeric_only=True).reset_index()  # Take the mean of duplicates
        data['Study'] = study.name
        all_data.append(data)
    return pd.concat(all_data, ignore_index=True)

def correlation_analysis(data: pd.DataFrame, genes: List[str]) -> Dict[str, Tuple[float, float]]:
    results = {}
    for gene in genes:
        tmed9_logfc = data[data['Gene.symbol'] == 'GAPDH']['normalized_logFC'].values
        gene_logfc = data[data['Gene.symbol'] == gene]['normalized_logFC'].values
        # Fill missing value in gene_logfc with mean of existing values
        if len(tmed9_logfc) > len(gene_logfc):
            diff_len = len(tmed9_logfc) - len(gene_logfc)
            gene_logfc = np.append(gene_logfc, [np.mean(gene_logfc)]*diff_len)
        
        elif len(tmed9_logfc) < len(gene_logfc):
            diff_len = len(gene_logfc) - len(tmed9_logfc)
            tmed9_logfc = np.append(tmed9_logfc, [np.mean(tmed9_logfc)]*diff_len)
        
        # print(tmed9_logfc)
        # print("#######################")
        # print(gene_logfc)
        correlation, p_value = pearsonr(tmed9_logfc, gene_logfc)
        results[gene] = (correlation, p_value)
    return results

def regression_analysis(data: pd.DataFrame, genes: List[str]) -> Dict[str, sm.regression.linear_model.RegressionResultsWrapper]:
    results = {}
    for gene in genes:
        tmed9_logfc = data[data['Gene.symbol'] == 'GAPDH']['normalized_logFC'].values
        gene_logfc = data[data['Gene.symbol'] == gene]['normalized_logFC'].values
        # Fill missing value in gene_logfc with mean of existing values
        if len(tmed9_logfc) > len(gene_logfc):
            diff_len = len(tmed9_logfc) - len(gene_logfc)
            gene_logfc = np.append(gene_logfc, [np.mean(gene_logfc)]*diff_len)
        
        elif len(tmed9_logfc) < len(gene_logfc):
            diff_len = len(gene_logfc) - len(tmed9_logfc)
            tmed9_logfc = np.append(tmed9_logfc, [np.mean(tmed9_logfc)]*diff_len)

        X = sm.add_constant(tmed9_logfc)
        model = sm.OLS(gene_logfc, X)
        results[gene] = model.fit()
    return results

def plot_regression_simple(data: pd.DataFrame, genes: List[str], studies: List[Study]):
    for gene in genes:
        # Extract logFC values for TMED9 and the gene of interest
        tmed9_logfc = data[data['Gene.symbol'] == 'TMED9']['normalized_logFC'].values
        gene_logfc = data[data['Gene.symbol'] == gene]['normalized_logFC'].values

        # Create a DataFrame for seaborn
        df = pd.DataFrame({
            'TMED9 Expression': tmed9_logfc,
            'Gene Expression': gene_logfc
        })

        # Create the plot
        sns.lmplot(x='TMED9 Expression', y='Gene Expression', data=df)

        plt.title(f'Regression of {gene} Expression on TMED9 Expression')
        plt.show()


    for gene in genes:
        # Create the regression plot with seaborn
        plot = sns.lmplot(x='TMED9 Expression', y='Gene Expression', data=pd.DataFrame(), line_kws={'color': 'red'}, scatter=False)

        ax = plot.axes[0,0]  # Access the underlying axes object

        # Define the color map for studies - 'tab20b' has darker colors
        cmap = plt.get_cmap('tab20b')
        study_colors = {study.name: cmap(i) for i, study in enumerate(studies)}

        for study in studies:
            study_data = data[data['Study'] == study.name]
            tmed9_logfc = study_data[study_data['Gene.symbol'] == 'TMED9']['normalized_logFC'].values
            gene_logfc = study_data[study_data['Gene.symbol'] == gene]['normalized_logFC'].values

            # Ensure there are non-empty arrays for both tmed9_logfc and gene_logfc
            if len(tmed9_logfc) == 0 or len(gene_logfc) == 0:
                print(f"Skipping {gene} for {study.name} due to insufficient data")
                continue

            # Equalize lengths
            if len(tmed9_logfc) > len(gene_logfc):
                diff_len = len(tmed9_logfc) - len(gene_logfc)
                gene_logfc = np.append(gene_logfc, [np.mean(gene_logfc)]*diff_len)
            elif len(tmed9_logfc) < len(gene_logfc):
                diff_len = len(gene_logfc) - len(tmed9_logfc)
                tmed9_logfc = np.append(tmed9_logfc, [np.mean(tmed9_logfc)]*diff_len)

            ax.scatter(tmed9_logfc, gene_logfc, alpha=0.7, edgecolors='w', label=study.name, color=study_colors[study.name])

        ax.set_title(f'Regression plot of TMED9 vs {gene}')
        ax.set_xlabel('TMED9 LogFC')
        ax.set_ylabel(f'{gene} LogFC')
        ax.legend(title='Study')
        plt.show()

def plot_regression(data: pd.DataFrame, genes: List[str], studies: List[Study]):
    for gene in genes:
        # Extract logFC values for TMED9 and the gene of interest across all studies
        tmed9_logfc = data[data['Gene.symbol'] == 'TMED9']['normalized_logFC'].values
        gene_logfc = data[data['Gene.symbol'] == gene]['normalized_logFC'].values

        # Equalize lengths
        if len(tmed9_logfc) > len(gene_logfc):
            diff_len = len(tmed9_logfc) - len(gene_logfc)
            gene_logfc = np.append(gene_logfc, [np.mean(gene_logfc)]*diff_len)
        elif len(tmed9_logfc) < len(gene_logfc):
            diff_len = len(gene_logfc) - len(tmed9_logfc)
            tmed9_logfc = np.append(tmed9_logfc, [np.mean(tmed9_logfc)]*diff_len)

        # Create a DataFrame for seaborn
        df = pd.DataFrame({
            'TMED9 Expression': tmed9_logfc,
            'Gene Expression': gene_logfc
        })

        # Create the regression plot with seaborn
        plot = sns.lmplot(x='TMED9 Expression', y='Gene Expression', data=df, line_kws={'color': 'red'}, scatter=False)

        ax = plot.axes[0,0]  # Access the underlying axes object

        # Define the color map for studies - 'tab20b' has darker colors
        cmap = plt.get_cmap('tab20b')
        study_colors = {study.name: cmap(i) for i, study in enumerate(studies)}

        for study in studies:
            study_data = data[data['Study'] == study.name]
            tmed9_logfc = study_data[study_data['Gene.symbol'] == 'TMED9']['normalized_logFC'].values
            gene_logfc = study_data[study_data['Gene.symbol'] == gene]['normalized_logFC'].values

            # Ensure there are non-empty arrays for both tmed9_logfc and gene_logfc
            if len(tmed9_logfc) == 0 or len(gene_logfc) == 0:
                print(f"Skipping {gene} for {study.name} due to insufficient data")
                continue

            # Equalize lengths
            if len(tmed9_logfc) > len(gene_logfc):
                diff_len = len(tmed9_logfc) - len(gene_logfc)
                gene_logfc = np.append(gene_logfc, [np.mean(gene_logfc)]*diff_len)
            elif len(tmed9_logfc) < len(gene_logfc):
                diff_len = len(gene_logfc) - len(tmed9_logfc)
                tmed9_logfc = np.append(tmed9_logfc, [np.mean(tmed9_logfc)]*diff_len)

            ax.scatter(tmed9_logfc, gene_logfc, alpha=0.7, edgecolors='w', label=study.name, color=study_colors[study.name])

        ax.set_title(f'Regression plot of TMED9 vs {gene}')
        ax.set_xlabel('TMED9 LogFC')
        ax.set_ylabel(f'{gene} LogFC')
        ax.legend(title='Study')
        plt.show()


    study_dict = {}
    data = compile_data(studies)
    regression_results = regression_analysis(data, genes)
    for study in studies:
        study_data = data[data['Study'] == study.name]

        for gene, result in regression_results.items():
            correlation = np.sqrt(result.rsquared) if 'TMED9' in study_data['Gene.symbol'].values and gene in study_data['Gene.symbol'].values else 0
            p_value = result.pvalues[1] if 'TMED9' in study_data['Gene.symbol'].values and gene in study_data['Gene.symbol'].values else 1
            study_dict[study.name] = (correlation, p_value)

    return study_dict

def plot_forest_plot(results: Dict[str, Tuple[float, float]]):

    fig, ax = plt.subplots()

    # Sort the studies by correlation for better visualization
    sorted_studies = sorted(results.items(), key=lambda x: x[1][0])

    study_names = [study[0] for study in sorted_studies]
    correlations = [result[0] for _, result in sorted_studies]
    p_values = [result[1] for _, result in sorted_studies]

    y_pos = range(len(study_names))

    ax.errorbar(correlations, y_pos, xerr=p_values, fmt='o')
    ax.set_yticks(y_pos)
    ax.set_yticklabels(study_names)
    ax.set_xlabel('Correlation')
    ax.set_title('Forest plot of correlation')

    plt.show()

def compile_data_new(studies: List[Study]) -> pd.DataFrame:
    all_data = []
    for study in studies:
        data = pd.read_csv(study.path, sep='\t')
        data = data.groupby('Gene.symbol').mean(numeric_only=True).reset_index()  # Take the mean of duplicates
        data['Study'] = study.name
        data['Country'] = study.country  # Assuming `study.country` returns the country of the study
        data['Tissue'] = study.tissue  # Assuming `study.tissue` returns the tissue type studied
        all_data.append(data)
    return pd.concat(all_data, ignore_index=True)

def correlation_analysis_new(data: pd.DataFrame, genes: List[str]) -> Dict[str, Dict[str, Tuple[float, float]]]:
    results = {}
    for gene in genes:
        results[gene] = {}
        gene_data = data[data['Gene.symbol'] == gene]
        for study in gene_data['Study'].unique():
            study_data = gene_data[gene_data['Study'] == study]
            country = study_data['Country'].iloc[0]
            tissue = study_data['Tissue'].iloc[0]
            tmed9_logfc = study_data[study_data['Gene.symbol'] == 'TMED9']['normalized_logFC'].values
            gene_logfc = study_data['normalized_logFC'].values
            
            if len(tmed9_logfc) > len(gene_logfc):
                diff_len = len(tmed9_logfc) - len(gene_logfc)
                gene_logfc = np.append(gene_logfc, [np.mean(gene_logfc)]*diff_len)
            
            elif len(tmed9_logfc) < len(gene_logfc):
                diff_len = len(gene_logfc) - len(tmed9_logfc)
                tmed9_logfc = np.append(tmed9_logfc, [np.mean(tmed9_logfc)]*diff_len)
            
            correlation, p_value = pearsonr(tmed9_logfc, gene_logfc)
            results[gene][(study, country, tissue)] = (correlation, p_value)
    return results

def main():
    print("Performing meta-analysis...")
    filename = 'Configs.txt'
    studies = read_config_file(filename)
    genes =  ["NFKB1", "SRR", "PDE4B", "IGF1R"]
    o_genes =  ["TMED9"]
    # data = compile_data_new(studies)
    # results = correlation_analysis_new(data, genes)

    # # Flatten results for visualization
    # flat_results = []
    # for gene, gene_results in results.items():
    #     for (study, country, tissue), (correlation, p_value) in gene_results.items():
    #         flat_results.append({
    #             'Gene': gene,
    #             'Study': study,
    #             'Country': country,
    #             'Tissue': tissue,
    #             'Correlation': correlation,
    #             'P-value': p_value
    #         })
    
    # df = pd.DataFrame(flat_results)
    
    # # Plotting correlation and p-value for each gene, per country and tissue
    # fig, axs = plt.subplots(4, 1, figsize=(10, 20))

    # # Correlation per Gene per Country
    # sns.barplot(data=df, x='Gene', y='Correlation', hue='Country', ax=axs[0])
    # axs[0].set_title('Correlation per Gene per Country')

    # # P-value per Gene per Country
    # sns.barplot(data=df, x='Gene', y='P-value', hue='Country', ax=axs[1])
    # axs[1].set_title('P-value per Gene per Country')

    # # Correlation per Gene per Tissue
    # sns.barplot(data=df, x='Gene', y='Correlation', hue='Tissue', ax=axs[2])
    # axs[2].set_title('Correlation per Gene per Tissue')

    # # P-value per Gene per Tissue
    # sns.barplot(data=df, x='Gene', y='P-value', hue='Tissue', ax=axs[3])
    # axs[3].set_title('P-value per Gene per Tissue')

    # plt.tight_layout()
    # plt.show()
 
    all_data = compile_data(studies)
    correlation_results = correlation_analysis(all_data, o_genes)
    regression_results = regression_analysis(all_data, o_genes)

    # Print the correlation and regression results for this study
    for gene in o_genes:
        print(f"For gene {gene}, correlation: {correlation_results[gene][0]}, p-value: {correlation_results[gene][1]}")
        print(f"For gene {gene}, regression:\n{regression_results[gene].summary()}")

    plot_regression(all_data, o_genes, studies)


    # values_dict, combined_scores = perform_meta_analysis(studies)
    # top_500_genes_indices = get_top_500_genes(combined_scores)
    
    # gene_list = []
    # counter = 1
    # for i, index in enumerate(top_500_genes_indices, 1):
    #     gene_symbol = str(list(values_dict.keys())[index])
    #     if gene_symbol == "nan" or gene_symbol.startswith("LOC") or ("///" in gene_symbol):
    #         continue
    #     print(f"{counter}. {gene_symbol}")
    #     counter += 1
    #     gene_list.append(gene_symbol)

if __name__ == "__main__":
    main()
