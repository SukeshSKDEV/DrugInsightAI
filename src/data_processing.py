import numpy as np
import pandas as pd

def validate_dataset(df):
    """
    Validate the dataset before preprocessing.
    """

    print("Validating dataset...")
    # 1. Dataset is not empty
    if df.empty:
        ValueError("Dataset is empty.")
    
    # 2. Required columns exist
    required_columns=['molecule_chembl_id','canonical_smiles','standard_value','logp','hba','hbd','molecular_weight']
    missing_columns=[column for column in required_columns if column not in df.columns]
    if missing_columns:
        raise ValueError(f"Missing columns: {missing_columns}")

    # 3. Duplicate column names? 
    if df.columns.duplicated().any():
        raise ValueError("Duplicated column names found.")

    # 4. standard_value is numeric
    if not pd.api.types.is_numeric_dtype(df['standard_value']):
        raise ValueError("standard_value must be numeric")       
    print("Dataset validation passed.")
    return df
    
def select_columns(df):
    """
    Select only the columns required for
    machine learning.
    """
    print("Selecting required columns...")    
    columns=['molecule_chembl_id','canonical_smiles','standard_value','logp','hba','hbd','molecular_weight']
    df=df[columns].copy()
    print(f"Selecting Columns : {len(df.columns)}")
    print(f"Selecting Rows : {len(df)}")
    return df

def remove_missing_value(df):
    
    """
    Remove rows containing missing values.
    """
    
    print(f"Rows before removing missing values: {len(df)}")
    row_before=len(df)
    df = df.dropna().copy()
    print(f"Rows after removing missing values: {len(df)}")
    row_after=len(df)
    print(f"Number of missing value : {row_before-row_after}")
    
    return df

def remove_exact_duplicates(df):
    """
    Remove completely identical rows.
    """

    print(f"Rows before removing duplicate values: {len(df)}")
    df = df.drop_duplicates().copy()
    print(f"Rows after removing duplicate values: {len(df)}")
    
    return df

def aggregate_duplicates(df):
    """
    Aggregate bioactivity measurements for
    the same molecule.
    """
    print(f"Rows before aggregating standard_values: {len(df)}")
    df=df.groupby(['molecule_chembl_id','canonical_smiles','logp','hba','hbd','molecular_weight'])['standard_value'].mean().reset_index()
    print(f"Rows after aggregating standard_values: {len(df)}")  
    return df

def convert_to_pic50(df):
    """
    Convert IC50 (nM) to pIC50.
    """
    # Keep only positive IC50 values
    df=df[df['standard_value']>0].copy()
    # Convert nM → M
    ic50_molar=df['standard_value']*1e-9
    # Calculate pIC50
    df["pIC50"]=-np.log10(ic50_molar)
    df.drop(columns=['standard_value'], inplace=True)
    return df

def clean_dataset(df):
    df=validate_dataset(df)

    df=select_columns(df)
    
    df=remove_missing_value(df)

    df=remove_exact_duplicates(df)

    df=aggregate_duplicates(df)

    print("Converting IC50 to pIC50...")
    df=convert_to_pic50(df)
    print("Done.")

    print("Cleaning completed successfully.")    
    return df
    
    
    