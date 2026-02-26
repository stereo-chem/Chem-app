import streamlit as st
import requests
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import AllChem, EnumerateStereoisomers
from rdkit.Chem.Draw import rdMolDraw2D
from stmol import showmol
import py3Dmol

st.set_page_config(page_title="Professional Isomer Analyzer", layout="wide")

def get_smiles_smart(name):
    try:
        opsin_url = f"https://opsin.ch.cam.ac.uk/opsin/{name}.json"
        opsin_res = requests.get(opsin_url)
        if opsin_res.status_code == 200: return opsin_res.json()['smiles']
    except: pass
    try:
        pcp_res = pcp.get_compounds(name, 'name')
        if pcp_res: return pcp_res[0].isomeric_smiles
    except: pass
    return None

st.title("🧪 Professional Isomer Analysis")
name_input = st.text_input("Enter Molecule Name:", "1,3-diphenylallene")

if st.button("Generate Analysis"):
    smiles = get_smiles_smart(name_input)
    if smiles:
        base_mol = Chem.MolFromSmiles(smiles)
        opts = EnumerateStereoisomers.StereoEnumerationOptions(tryEmbedding=True, onlyUnassigned=False)
        isomers = list(EnumerateStereoisomers.EnumerateStereoisomers(base_mol, options=opts))
        
        cols = st.columns(len(isomers))
        for i, iso in enumerate(isomers):
            with cols[i]:
                # --- 2D Professional Rendering ---
                d2d = rdMolDraw2D.MolDraw2DCairo(400, 400)
                d_opts = d2d.drawOptions()
                d_opts.bondLineWidth = 3
                d_opts.addAtomIndices = False
                
                Chem.AssignStereochemistry(iso, force=True, cleanIt=True)
                d2d.DrawMolecule(iso)
                d2d.FinishDrawing()
                
                st.image(d2d.GetDrawingText(), caption=f"Isomer {i+1}", use_container_width=True)
                
                # --- 3D Interactive ---
                iso_3d = Chem.AddHs(iso)
                AllChem.EmbedMolecule(iso_3d, AllChem.ETKDG())
                view = py3Dmol.view(width=400, height=300)
                view.addModel(Chem.MolToMolBlock(iso_3d), 'mol')
                view.setStyle({'stick': {'radius': 0.15}, 'sphere': {'scale': 0.25}})
                view.zoomTo()
                showmol(view, height=300, width=400)
    else:
        st.error("Molecule not found.")
