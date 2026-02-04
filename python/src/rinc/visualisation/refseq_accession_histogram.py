'''
Generate a refseq 
Created on Feb 3, 2026

@author: pleyte
'''

import pandas as pd
import matplotlib.pyplot as plt

def visualisation_accession_historgram(df: pd.DataFrame):
    # 1. Split 'NM_123.1' into 'NM_123' and '1'
    df['base_ac'] = df['accession'].str.split('.').str[0]
    
    # 2. Count how many versions exist for each base accession
    counts = df.groupby('base_ac').size()
    
    # 3. Create the histogram
    plt.figure(figsize=(10, 6))
    plt.hist(counts, bins=range(1, counts.max() + 2), align='left', rwidth=0.8, color='skyblue', edgecolor='black')
    
    plt.title('Distribution of How many RefSeq Versions Each Base Accession Has (hg19)')
    plt.xlabel('Number of Versions Found in UTA')
    plt.ylabel('Count of Base Accessions')
    
    plt.xticks(range(1, counts.max() + 1))
    plt.grid(axis='y', alpha=0.3)
    
    stats_text = f'Mean: {counts.mean():.2f}\nMedian: {int(counts.median())}'
    # Place it in the 'Upper Right' (0.95, 0.95) using 'Axes Coordinates'
    # transform=plt.gca().transAxes ensures (0,0) is bottom-left and (1,1) is top-right
    plt.text(0.95, 0.95, stats_text, transform=plt.gca().transAxes, 
             fontsize=12, verticalalignment='top', horizontalalignment='right',
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))
    
    
    plt.show()    

def visualisation_accession_version_counts(df: pd.DataFrame):
    """
    How many accessions have a .1 version? A .2? etc 
    """
    # 1. Discard any rows where 'accession' contains a slash (eg "NM_001202404.1/111..1620")
    df = df[~df['accession'].str.contains('/', na=False)].copy()

    # 2. Extract the version number (the part after the '.')
    # We convert to int so the X-axis sorts numerically (1, 2, 3...) rather than alphabetically
    df['version_num'] = df['accession'].str.split('.').str[1].astype(int)

    # 3. Count how many times each version number appears
    version_distribution = df['version_num'].value_counts().sort_index()

    # 4. Plotting
    plt.figure(figsize=(10, 6))
    version_distribution.plot(kind='bar', color='salmon', edgecolor='black')

    plt.title('Frequency of RefSeq Version Suffixes in hg19')
    plt.xlabel('Version Suffix (e.g., .1, .2, .3)')
    plt.ylabel('Number of Accessions')  
    plt.yscale('log')

    plt.grid(axis='y', alpha=0.3)
    plt.show()    

def main():
    # Load your data
    df = pd.read_csv('/tmp/monday/hg19_refseq_transcripts/uta_20240523b_nm_refseq_transcripts_splign_hg19.csv') 
    
    visualisation_accession_historgram(df)
    
    visualisation_accession_version_counts(df)

if __name__ == '__main__':
    main()