import pandas as pd
from chembl_webresource_client.new_client import new_client
import time
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from reportlab.platypus import SimpleDocTemplate, Image, Paragraph, Spacer
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib import colors
from reportlab.lib.units import inch
from reportlab.lib.pagesizes import A4
import os

# Step 1: Connect to ChEMBL API and fetch compounds
def fetch_chembl_compounds(limit=100):
    """
    Fetch compounds from ChEMBL database
    limit: number of compounds to fetch (start small for testing)
    """
    print(f"Fetching {limit} compounds from ChEMBL...")
    
    molecule = new_client.molecule
    compounds = molecule.filter(prodrug=0)[:limit]
    
    # Convert to list of dictionaries
    compounds_list = []
    for comp in compounds:
        compound_data = {
            'chembl_id': comp.get('molecule_chembl_id'),
            'pref_name': comp.get('pref_name'),
            'molecule_type': comp.get('molecule_type'),
            'max_phase': comp.get('max_phase'),
            'molecular_weight': comp.get('molecule_properties', {}).get('mw_freebase'),
            'alogp': comp.get('molecule_properties', {}).get('alogp'),
            'hbd': comp.get('molecule_properties', {}).get('hbd'),
            'hba': comp.get('molecule_properties', {}).get('hba'),
            'ro3_pass': comp.get('molecule_properties', {}).get('ro3_pass'),
            'num_ro5_violations': comp.get('molecule_properties', {}).get('num_ro5_violations'),
            'smiles': comp.get('molecule_structures', {}).get('canonical_smiles'),
            'inchi': comp.get('molecule_structures', {}).get('standard_inchi'),
            'inchi_key': comp.get('molecule_structures', {}).get('standard_inchi_key'),
        }
        compounds_list.append(compound_data)
    
    df = pd.DataFrame(compounds_list)
    print(f"Fetched {len(df)} compounds")
    return df

# Step 2: Store data in Excel
def store_in_excel(df, filename="chembl_compounds.xlsx"):
    """
    Store DataFrame in Excel file
    """
    print(f"Storing {len(df)} compounds in Excel...")
    
    df.to_excel(filename, index=False)
    
    print(f"Data stored successfully in {filename}!")

# Step 3: Printing 2D structures
def save_structures_to_pdf(df, filename="chembl_structures.pdf", max_molecules=100):
    print("Starting PDF generation...")

    doc = SimpleDocTemplate(filename, pagesize=A4)
    elements = []
    styles = getSampleStyleSheet()

    os.makedirs("temp_images", exist_ok=True)

    count = 0

    for idx, row in df.iterrows():

        smiles = row.get("smiles")
        name = row.get("pref_name") or row.get("chembl_id")

        if pd.isna(smiles):
            print("Skipped: Missing SMILES")
            continue

        mol = Chem.MolFromSmiles(smiles)

        if mol is None:
            print(f"Invalid SMILES: {name}")
            continue
        rdDepictor.Compute2DCoords(mol)
         # High resolution canvas
        width, height = 1500, 1500
        drawer = rdMolDraw2D.MolDraw2DCairo(width, height)
        drawer.drawOptions().addAtomIndices = False
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()

        img_path = f"temp_images/{name}.png"

        Draw.MolToFile(mol, img_path, size=(1500, 1500))

        elements.append(Paragraph(f"<b>{name}</b>", styles["Heading3"]))
        elements.append(Image(img_path, width=3.75*inch, height=3.75*inch))
        
        count += 1

        if count >= max_molecules:
            break

    if len(elements) == 0:
        print("No valid molecules found. PDF not created.")
        return

    print(f"Adding {count} molecules to PDF...")
    doc.build(elements)
    print("PDF successfully created!")

# Main execution
def main():
    print("Starting program...")
    
    df_compounds = fetch_chembl_compounds(limit=100)
    
    if df_compounds.empty:
        print("DataFrame is empty!")
        return
    
    print("Fetched:", df_compounds.shape)
    
    df_compounds = df_compounds.dropna(subset=['chembl_id'])
    
    print("After cleaning:", df_compounds.shape)
    
    print("Saving to Excel...")
    df_compounds.to_excel("chembl_compounds.xlsx", index=False)
    save_structures_to_pdf(df_compounds)
    print("Done! Check your project folder.")

main()
