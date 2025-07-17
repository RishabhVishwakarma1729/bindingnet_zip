# app.py

import streamlit as st
import pandas as pd
from io import BytesIO
import os
import zipfile

from rdkit import Chem
from rdkit.Chem import rdmolfiles

# ---- Helper Functions ----

def get_protein_id(protein_filename):
    # Strip extension; you may want to use just first token (before _ etc.) if needed
    return os.path.splitext(protein_filename)[0]

def ligand_to_mol2(ligand_file, ext):
    """Convert ligand to MOL2 format, supports mol2, sdf, pdb using RDKit."""
    # Accept: ext (string with .mol2/.sdf/.pdb)
    # Returns: bytes of .mol2 file, or None if failed
    if ext.lower() == ".mol2":
        return ligand_file.read()
    elif ext.lower() == ".sdf":
        try:
            suppl = Chem.ForwardSDMolSupplier(ligand_file)
            mols = [m for m in suppl if m is not None]
        except Exception:
            return None
    elif ext.lower() == ".pdb":
        try:
            pdb_block = ligand_file.read().decode("utf-8")
            mol = Chem.MolFromPDBBlock(pdb_block, sanitize=True, removeHs=False)
            mols = [mol] if mol is not None else []
        except Exception:
            return None
    else:
        mols = []

    if not mols:
        return None
    out = BytesIO()
    writer = rdmolfiles.Mol2Writer(out)
    for m in mols:
        if m is not None:
            writer.write(m)
    writer.flush()
    return out.getvalue()

def make_zip(filename, content):
    mem_zip = BytesIO()
    with zipfile.ZipFile(mem_zip, "w", zipfile.ZIP_DEFLATED) as zf:
        zf.writestr(filename, content)
    mem_zip.seek(0)
    return mem_zip

# ---- Streamlit UI ----

st.title("Protein-Ligand Zip File Preparer")

st.markdown(
    """
    **Instructions:**  
    1. Upload a **protein file** (.pdb or .mol2)  
    2. Upload a **ligand file** (.mol2, .sdf, or .pdb)  
    3. Download the packaged zip files and the metadata CSV  
    """
)

protein_up = st.file_uploader(
    "Protein file (.pdb, .mol2)", type=["pdb", "mol2"], key="protein"
)
ligand_up = st.file_uploader(
    "Ligand file (.pdb, .mol2, .sdf)", type=["pdb", "mol2", "sdf"], key="ligand"
)

if protein_up and ligand_up:
    protein_filename = protein_up.name
    protein_id = get_protein_id(protein_filename)
    ligand_ext = os.path.splitext(ligand_up.name)[1]

    # Read protein file (as is, you may want to do conversion if necessary)
    protein_content = protein_up.read()

    # Ligand: convert to mol2
    ligand_mol2 = ligand_to_mol2(ligand_up, ligand_ext)
    if ligand_mol2 is None:
        st.error("Failed to convert ligand to MOL2. Please check file format/content.")
        st.stop()

    # "Pocket" content: here just duplicating protein as placeholder
    # Replace with real pocket extraction if you have it
    pocket_content = protein_content

    # Naming convention
    pocket_name = f"{protein_id}_Pocket.mol2"
    ligand_name = f"{protein_id}_Ligand.mol2"
    protein_name = f"{protein_id}.mol2"

    # Create zips
    zip_pocket = make_zip(pocket_name, pocket_content)
    zip_ligand = make_zip(ligand_name, ligand_mol2)
    zip_protein = make_zip(protein_name, protein_content)

    # CSV
    selected_df = pd.DataFrame([{
        "proteinID": protein_id,
        "ligandName": os.path.splitext(ligand_up.name)[0],
        "status": "selected"
    }])
    csv_bytes = selected_df.to_csv(index=False).encode("utf-8")

    st.success("Files created! Please download below:")

    st.download_button(
        label=f"Download {pocket_name}.zip",
        data=zip_pocket,
        file_name=f"{protein_id}_Pocket.zip"
    )
    st.download_button(
        label=f"Download {ligand_name}.zip",
        data=zip_ligand,
        file_name=f"{protein_id}_Ligand.zip"
    )
    st.download_button(
        label=f"Download {protein_name}.zip",
        data=zip_protein,
        file_name=f"{protein_id}.zip"
    )
    st.download_button(
        label="Download selected.csv",
        data=csv_bytes,
        file_name="selected.csv"
    )
else:
    st.info("Please upload both a protein and a ligand file to proceed.")
