import numpy as np
import pandas as pd
import requests
import os

# os.path.exists(cache_file)

Base_url="https://www.ebi.ac.uk/chembl/api/data"

def search_target(target):
    """
        Search a human SINGLE PROTEIN target in ChEMBL.
        Returns the target CHEMBL ID.
    """
    
    url=f"{Base_url}/target/search.json?q={target}"

    response=requests.get(url)

    response.raise_for_status()

    data=response.json()

    for target in data["targets"]:
        if (
            target["organism"]=="Homo sapiens" 
            and target["target_type"]=="SINGLE PROTEIN"
        ) :
            return {
                'target_chembl_id': target['target_chembl_id'],
                'organism': target['organism'],
                'pref_name':target['pref_name']          
            }

    return None


def get_bioactivity(target_chembl_id,standard_type='IC50',limit=1000):

    """
            Retrieve bioactivity records for a given ChEMBL target.
        
            Parameters:
                target_chembl_id (str): Target ChEMBL ID.
                standard_type (str): Bioactivity measurement type (default: IC50).
                limit (int): Number of records per API page.
        
            Returns:
                pandas.DataFrame: Bioactivity records.
    """
    
    url=(
            f"{Base_url}/activity.json"
            f"?target_chembl_id={target_chembl_id}"
            f"&standard_type={standard_type}"
            f"&limit={limit}"
        )
    
    
    
    all_records=[]

    while url:
        
        print(url)
        
        response=requests.get(url)
        response.raise_for_status()
        
        data=response.json()      

        all_records.extend(data['activities'])

        next_page=data['page_meta']['next']

        print("Next page:", next_page)

        if next_page:
            url='https://www.ebi.ac.uk'+next_page
        else:
            url=None

    df=pd.DataFrame(all_records)
    print("Bioactivity shape : ",df.shape)
    print("Bioactivity unique moelcules : ",df["molecule_chembl_id"].nunique())

    return df


def get_molecule_properties(molecule_chembl_id):
    url=f"{Base_url}/molecule/{molecule_chembl_id}.json"

    response=requests.get(url)
    response.raise_for_status()

    data=response.json()
    
    molecule=data
    properties=molecule['molecule_properties']

    return {
        'molecule_chembl_id':molecule['molecule_chembl_id'],
        'logp':properties['alogp'],
        'hba':properties['hba'],
        'hbd':properties['hbd'],
        'molecular_weight':properties['full_mwt']
    }


def get_all_molecule_properties(bioactivity_df, cache_file='../data/cache/descriptor_cache.csv', max_molecules=None):
    
    # Create cache folder if it doesn't exist
    os.makedirs(os.path.dirname(cache_file), exist_ok=True)
    
    # Load cache if it exists
    if os.path.exists(cache_file) and os.path.getsize(cache_file) > 0 :
        print("Loading cache...")
        cache_df=pd.read_csv(cache_file)
    else:
        print("No cache found. Creating new cache.")
        cache_df=pd.DataFrame(columns=['molecule_chembl_id',
                                        'logp',
                                        'hba',
                                        'hbd',
                                        'molecular_weight'])

    # Create a set for fast lookup
    cached_ids=set(cache_df['molecule_chembl_id'])
    print(f"Cached molecules : {len(cached_ids)}")

    unique_molecules=(bioactivity_df['molecule_chembl_id'].dropna().unique())
    

    if max_molecules is not None:
        unique_molecules=unique_molecules[:max_molecules]
    print(f"Unique molecules : {len(unique_molecules)}")

    missing_molecules=[molecule for molecule in unique_molecules if molecule not in cached_ids]
    print(f"Molecules to download : {len(missing_molecules)}")

    new_descriptors=[]

    for i, molecule_id in enumerate(missing_molecules, start=1):
        try:
            compound=get_molecule_properties(molecule_id)
            new_descriptors.append(compound)
        except Exception as e:
            print(f"Failed : {molecule_id}")

    if new_descriptors:
        df=pd.DataFrame(new_descriptors)
        cache_df=pd.concat([cache_df,df], ignore_index=True)
        cache_df.to_csv(cache_file, index=False)

    descriptor_df = cache_df[
        cache_df["molecule_chembl_id"].isin(unique_molecules)
        ].copy()
    
    return descriptor_df
        
        
    
def merge_dataset(bioactivity_df, descriptor_df):
    final_df=pd.merge(bioactivity_df,descriptor_df,on='molecule_chembl_id',how='inner')
    return final_df

def build_dataset(target_name, max_molecules=None):
    
    print("Searching target...")
    target=search_target(target_name)

    if target is None:
        print("Target not found....")
        return None

    print(f"Target Found : {target['pref_name']}")
    print(f"CHEMBL ID    : {target['target_chembl_id']}")

    print("Downloading bioactivity...")
    bioactivity_df=get_bioactivity(target['target_chembl_id'])

    print("Bioactivity shape :", bioactivity_df.shape)
    print("Unique molecules :", bioactivity_df["molecule_chembl_id"].nunique())

    print("\nDownloading molecular properties...")
    descriptor_df=get_all_molecule_properties(bioactivity_df, max_molecules=max_molecules)

    print("Bioactivity Shape:", bioactivity_df.shape)
    print("Descriptor Shape:", descriptor_df.shape)

    print("Merging datasets...")
    final_df=merge_dataset(bioactivity_df,descriptor_df)
    
    print("Done!")
    
    return final_df


        
        

    

    

    
    