import streamlit as st
import requests
from rdkit import Chem
from rdkit.Chem import AllChem, EnumerateStereoisomers, Draw
from stmol import showmol
import py3Dmol

st.set_page_config(page_title="IUPAC Isomer Analyzer", layout="wide")
st.title("🧪 IUPAC Name to 3D Isomer Analysis")

# دالة لتحويل الاسم إلى SMILES باستخدام OPSIN
def name_to_smiles(name):
    try:
        url = f"https://opsin.ch.cam.ac.uk/opsin/{name}.json"
        response = requests.get(url)
        if response.status_code == 200:
            return response.json()['smiles']
        else:
            return None
    except:
        return None

name_input = st.text_input("Enter IUPAC Name (e.g., 1,3-diphenylpropadiene):", "1,3-diphenylallene")

if st.button("Analyze Molecule"):
    # تحويل الاسم إلى SMILES أولاً
    smiles = name_to_smiles(name_input)
    
    if not smiles:
        st.error("Could not interpret this name. Try a more formal IUPAC name.")
    else:
        base_mol = Chem.MolFromSmiles(smiles)
        
        # توليد الأيزومرات
        opts = EnumerateStereoisomers.StereoEnumerationOptions(tryEmbedding=True, onlyUnassigned=False)
        isomers = list(EnumerateStereoisomers.EnumerateStereoisomers(base_mol, options=opts))
        
        st.success(f"Successfully interpreted '{name_input}' and found {len(isomers)} isomers.")

        # عرض النتائج في شبكة (Grid)
        cols = st.columns(len(isomers))
        
        for i, iso in enumerate(isomers):
            with cols[i]:
                st.subheader(f"Isomer #{i+1}")
                
                # 1. رسم 2D
                # ملاحظة: الألينات كيرالية محورياً، RDKit يظهرها بوضوح في الـ 2D
                img = Draw.MolToImage(iso, size=(400, 400))
                st.image(img, caption=f"2D Structure for Isomer {i+1}")
                
                # 2. عرض 3D
                iso_3d = Chem.AddHs(iso)
                AllChem.EmbedMolecule(iso_3d, AllChem.ETKDG())
                AllChem.MMFFOptimizeMolecule(iso_3d)
                
                view = py3Dmol.view(width=400, height=300)
                view.addModel(Chem.MolToMolBlock(iso_3d), 'mol')
                view.setStyle({'stick': {'radius': 0.2}, 'sphere': {'scale': 0.3}})
                view.zoomTo()
                showmol(view, height=300, width=400)

