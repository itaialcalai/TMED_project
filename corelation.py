import pandas as pd
import numpy as np
import statsmodels.api as sm
from scipy.stats import pearsonr

def load_gene_expression_data(path):
    df = pd.read_csv(path)
    return df

def correlation_test(df, gene1, gene2):
    expression_gene1 = df[gene1]
    expression_gene2 = df[gene2]
    corr, _ = pearsonr(expression_gene1, expression_gene2)
    print(f'Pearsons correlation between {gene1} and {gene2}: %.3f' % corr)

def regression_analysis(df, dependent_var, independent_var):
    X = df[independent_var]
    y = df[dependent_var]

    X = sm.add_constant(X)  # adding a constant

    model = sm.OLS(y, X).fit()
    print_model = model.summary()
    print(print_model)

def main():
    # Load your gene expression data, which should have one row per sample and one column per gene
    # Replace this path with the path to your data
    path_to_gene_expression_data = "gene_expression_data.csv"  
    df = load_gene_expression_data(path_to_gene_expression_data)

    # Run correlation test
    correlation_test(df, "SigGene1", "TMED9")

    # Run regression analysis
    regression_analysis(df, "SigGene1", "TMED9")

if __name__ == "__main__":
    main()
