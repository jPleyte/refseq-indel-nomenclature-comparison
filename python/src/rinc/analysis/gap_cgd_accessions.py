'''
See which CGD transcripts are known to have refseq mismatch. 

1) Compare file1 with file2

2) Because a lot of the transcripts (1918 of them) exported from CGD are accession-only (no version) I also want to do this coparison with the ".x" version trimmed
Created on Jan 23, 2026

This tool is not currently part of the workflow. 
@author: pleyte
'''
import argparse
import pandas as pd

def compare_transcripts(gap_file, cgd_transcript_file, cgd_variants_file, output_xlsx):
    # 1. Load the data
    gap_df = pd.read_csv(gap_file)
    cgd_tx_df = pd.read_csv(cgd_transcript_file)
    variants_df = pd.read_csv(cgd_variants_file)
    
    # Identify accession columns
    gap_col = 'accession' 
    assert 'accession' in gap_df.columns 
    cgd_tx_col = 'accession' 
    assert 'accession' in cgd_tx_df.columns 
    
    # 2. Prepare lookup sets for speed
    # Exact: 'NM_000123.4'
    gap_set_exact = set(gap_df[gap_col].unique())
    # Versionless: 'NM_000123'
    gap_set_versionless = set(gap_df[gap_col].str.split('.').str[0].unique())

    # 3. Build Transcript Comparison (Sheet 1)
    tx_comparison = pd.DataFrame()
    tx_comparison['CGD_Transcript'] = cgd_tx_df[cgd_tx_col]
    
# Create the check flags
    cgd_tx_stripped = tx_comparison['CGD_Transcript'].str.split('.').str[0]
    tx_comparison['Matches_Gap_Base_Accession'] = cgd_tx_stripped.isin(gap_set_versionless)
    tx_comparison['Matches_Gap_Exact_Version'] = tx_comparison['CGD_Transcript'].isin(gap_set_exact)
    
    # 4. Filter for Gapped Transcripts
    # We want any transcript where EITHER match is True to trigger the variant pull
    gapped_tx_list = tx_comparison[
        (tx_comparison['Matches_Gap_Base_Accession']) | 
        (tx_comparison['Matches_Gap_Exact_Version'])
    ]['CGD_Transcript'].unique()

    # 5. Filter the Variants Dataframe (Sheet 2)
    # Using 'cdna_transcript' as the join key from the clinical file
    affected_variants_df = variants_df[variants_df['cdna_transcript'].isin(gapped_tx_list)].copy()

    # 6. Save to Excel
    with pd.ExcelWriter(output_xlsx, engine='openpyxl') as writer:
        # Sheet 1: The mapping logic
        tx_comparison.to_excel(writer, sheet_name='Transcript_Check', index=False)
        
        # Sheet 2: All original rows for variants on those gapped transcripts
        affected_variants_df.to_excel(writer, sheet_name='Affected_Variants', index=False)

    # Summary Stats
    print(f"--- Variant Analysis Complete ---")
    print(f"Total Unique Transcripts checked: {len(tx_comparison)}")
    print(f"Transcripts found with Gaps:      {len(gapped_tx_list)}")
    print(f"Total Variant Rows Filtered:      {len(affected_variants_df)}")
    print(f"Results written to: {output_xlsx}")    

def main():
    
    parser = argparse.ArgumentParser(description='Identify CGD''s reported transcripts that are known to have reference gaps')
    parser.add_argument('--gap_accessions', help="Accessions known to have reference mismatch (csv)", required=True)
    parser.add_argument('--cgd_accessions', help="Accessions associated with reported variants in CGD (csv)", required=True)
    parser.add_argument('--cgd_variant_transcripts', help="Variants and transcript accessions exported from CGD (csv)", required=True)
    parser.add_argument('--out_xlsx', help="Save results to file (xlsx)", required=True)
    
    args = parser.parse_args()
    
    compare_transcripts(args.gap_accessions, 
                        args.cgd_accessions,
                        args.cgd_variant_transcripts, 
                        args.out_xlsx)
    
    
    

    
if __name__ == '__main__':
    main()