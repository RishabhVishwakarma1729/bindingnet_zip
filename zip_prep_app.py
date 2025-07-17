import streamlit as st
import pandas as pd
from io import BytesIO
import os
import zipfile

from rdkit import Chem
from rdkit.Chem import rdmolfiles, AllChem

# ---- Helper Functions ----

def get_protein_id(protein_file):
    """Extracts protein ID from header or filename."""
    try:
        header = protein_file.read(2048).decode('utf-8', errors='ignore').splitlines()
        for line in header:
            if line.startswith("HEADER") or line.startswith("COMPND") or line.startswith("TITLE"):
                tokens = line.strip().split()
                for token in tokens:
                    if len(token) == 4 and token.isalnum():
                        protein_file.seek(0)
                        return token
    except Exception:
        pass
    protein_file.seek(0)
    return os.path.splitext(protein_file.name)[0]

def ligand_to_mol2(ligand_file, ext):
    """Convert ligand to MOL2 format using RDKit with better diagnostics."""
    mols = []

    try:
        if ext.lower() == ".mol2":
            # Just read it directly
            content = ligand_file.read()
            if not content:
                st.error("Uploaded MOL2 file is empty.")
                return None
            return content

        elif ext.lower() == ".sdf":
            suppl = Chem.ForwardSDMolSupplier(ligand_file)
            mols = [m for m in suppl if m is not None]
            if not mols:
                st.error("RDKit could not read any valid molecule from the SDF file.")
                return None

        elif ext.lower() == ".pdb":
            pdb_block = ligand_file.read().decode("utf-8")
            mol = Chem.MolFromPDBBlock(pdb_block, sanitize=True, removeHs=False)
            if mol is None:
                st.error("RDKit could not parse the PDB block.")
                return None
            mol = Chem.AddHs(mol)
            mols = [mol]

        else:
            st.error(f"Unsupported ligand file type: {ext}")
            return None

        if not mols:
            st.error("No valid molecule found in the ligand file.")
            return None

        # Write to MOL2 format
        out = BytesIO()
        writer = rdmolfiles.Mol2Writer(out)
        for m in mols:
            if m is not None:
                writer.write(m)
        writer.flush()
        return out.getvalue()

    except Exception as e:
        st.error(f"Exception during ligand conversion: {e}")
        return None

def make_zip(filename, content):
    mem_zip = BytesIO()
    with zipfile.ZipFile(mem_zip, "w", zipfile.ZIP_DEFLATED) as zf:
        zf.writestr(filename, content)
    mem_zip.seek(0)
    return mem_zip

# ---- Streamlit UI ----

st.title("Protein-Ligand Zip File Preparer")

st.markdown("""
**Instructions:**  
1. Upload a **protein file** (.pdb or .mol2)  
2. Upload a **ligand file** (.mol2, .sdf, or .pdb)  
3. Download the packaged zip files and the metadata CSV  
""")

protein_up = st.file_uploader("Protein file (.pdb, .mol2)", type=["pdb", "mol2"], key="protein")
ligand_up = st.file_uploader("Ligand file (.pdb, .mol2, .sdf)", type=["pdb", "mol2", "sdf"], key="ligand")

if protein_up and ligand_up:
    protein_id = get_protein_id(protein_up)
    ligand_ext = os.path.splitext(ligand_up.name)[1]

    protein_content = protein_up.read()
    ligand_mol2 = ligand_to_mol2(ligand_up, ligand_ext)
    if ligand_mol2 is None:
        st.error("Failed to convert ligand to MOL2. Please check file format/content.")
        st.stop()

    # Extract ligand mol from converted MOL2
    mol2_reader = rdmolfiles.Mol2MolSupplier(BytesIO(ligand_mol2))
    ligand_mol = next((m for m in mol2_reader if m is not None), None)
    if ligand_mol is None:
        st.error("Could not parse converted MOL2 ligand.")
        st.stop()

    pocket_content = extract_pocket(protein_content, ligand_mol)
    if pocket_content is None:
        st.warning("Pocket extraction failed. Using full protein as fallback.")
        pocket_content = protein_content

    # Naming
    pocket_name = f"{protein_id}_Pocket.mol2"
    ligand_name = f"{protein_id}_Ligand.mol2"
    protein_name = f"{protein_id}.mol2"

    # Create ZIPs
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

    st.success("✅ Files created! Please download below:")

    st.download_button(f"📥 Download {pocket_name}.zip", zip_pocket, f"{protein_id}_Pocket.zip")
    st.download_button(f"📥 Download {ligand_name}.zip", zip_ligand, f"{protein_id}_Ligand.zip")
    st.download_button(f"📥 Download {protein_name}.zip", zip_protein, f"{protein_id}.zip")
    st.download_button("📄 Download selected.csv", csv_bytes, "selected.csv")

else:
    st.info("Please upload both a protein and a ligand file to proceed.")
