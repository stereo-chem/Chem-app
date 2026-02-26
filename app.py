import streamlit as st
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import AllChem, EnumerateStereoisomers, Draw
from stmol import showmol
import py3Dmol
from io import BytesIO

st.set_page_config(page_title="Chemical Isomer Explorer", layout="wide")

# تصميم العنوان بشكل جذاب
st.title("🧪 Comprehensive Isomer Analysis")
st.markdown("---")

name = st.text_input("Enter Molecule Name:", "2-butanol")

if st.button("Generate & Analyze All Isomers"):
    try:
        results = pcp.get_compounds(name, 'name')
        if not results:
            st.error("Molecule not found!")
        else:
            smiles = results[0].isomeric_smiles
            base_mol = Chem.MolFromSmiles(smiles)
            
            # 1. تحليل العلاقات الأيزومرية (Isomeric Relationships)
            st.subheader("1. Isomeric Relationships")
            opts = EnumerateStereoisomers.StereoEnumerationOptions(tryEmbedding=True, onlyUnassigned=False)
            isomers = list(EnumerateStereoisomers.EnumerateStereoisomers(base_mol, options=opts))
            
            st.info(f"💡 Relationships Analysis: {len(isomers)} Stereoisomeric relationship(s) detected.")

            # 2. عرض الـ 2D Structure Grid
            st.subheader("2. 2D Structure Grid")
            cols_2d = st.columns(len(isomers))
            
            for i, iso in enumerate(isomers):
                with cols_2d[i]:
                    # استخراج التسمية (R/S)
                    Chem.AssignStereochemistry(iso, force=True, cleanIt=True)
                    centers = Chem.FindMolChiralCenters(iso)
                    label = centers[0][1] if centers else "N/A"
                    
                    # توليد صورة 2D
                    img = Draw.MolToImage(iso, size=(300, 300))
                    st.image(img, caption=f"Isomer {i+1}: {label}", use_container_width=True)

            st.divider()

            # 3. عرض الـ 3D Visualization
            st.subheader("3. Interactive 3D Visualization")
            cols_3d = st.columns(len(isomers))

            for i, iso in enumerate(isomers):
                with cols_3d[i]:
                    # استخراج التسمية للـ 3D
                    centers = Chem.FindMolChiralCenters(iso)
                    label = centers[0][1] if centers else ""
                    st.markdown(f"**Isomer {i+1}: {label}**")
                    
                    # تحضير الشكل ثلاثي الأبعاد
                    iso_3d = Chem.AddHs(iso)
                    AllChem.EmbedMolecule(iso_3d, AllChem.ETKDG())
                    AllChem.MMFFOptimizeMolecule(iso_3d)
                    
                    view = py3Dmol.view(width=400, height=300)
                    mol_block = Chem.MolToMolBlock(iso_3d)
                    view.addModel(mol_block, 'mol')
                    view.setStyle({'stick': {'radius': 0.2}, 'sphere': {'scale': 0.3}})
                    view.zoomTo()
                    showmol(view, height=300, width=400)

    except Exception as e:
        st.error(f"Error during analysis: {e}")
