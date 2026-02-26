import streamlit as st
import pubchempy as pcp
import requests  # إضافة مكتبة الـ API
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers
from stmol import showmol
import py3Dmol

# 1. إعدادات الصفحة
st.set_page_config(page_title="Advanced Chemical Isomer Analysis", layout="wide")

# --- 🟢 إضافة دالة البحث الذكي (تدمج IUPAC و Common) ---
def get_smiles_smart(name):
    # أولاً: البحث عبر محرك OPSIN (للأسماء العلمية الدقيقة)
    try:
        opsin_url = f"https://opsin.ch.cam.ac.uk/opsin/{name}.json"
        opsin_res = requests.get(opsin_url)
        if opsin_res.status_code == 200:
            return opsin_res.json()['smiles']
    except:
        pass
    
    # ثانياً: البحث عبر PubChem (للأسماء الشائعة)
    try:
        pcp_res = pcp.get_compounds(name, 'name')
        if pcp_res:
            return pcp_res[0].isomeric_smiles
    except:
        pass
    return None
# ----------------------------------------------------

# 2. تصميم الواجهة (CSS لضمان Light Mode)
st.markdown("""
<style>
    .stApp { background-color: white; color: black; }
    .reportview-container { background: white; }
</style>
<h2 style='color: #800000; font-family: serif; border-bottom: 2px solid #dcdde1;'>Chemical Isomer Analysis System 2.0</h2>
<div style="background-color: #f9f9f9; padding: 15px; border: 1px solid #e1e1e1; border-left: 4px solid #800000; margin-bottom: 20px; font-family: sans-serif;">
    <strong style="color: #800000;">Stereoisomerism Reference Guide:</strong><br>
    1. <b style="color: #b22222;">Cis / Trans:</b> Identical groups on same/opposite sides.<br>
    2. <b style="color: #b22222;">E / Z (CIP System):</b> <b>Z (Zusammen)</b> together, <b>E (Entgegen)</b> opposite.<br>
    3. <b style="color: #b22222;">R / S (Optical):</b> Absolute configuration of chiral centers.
</div>
""", unsafe_allow_html=True)

# دالة لعرض المركب 3D
def render_3d(mol, title):
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    mblock = Chem.MolToMolBlock(mol)
    view = py3Dmol.view(width=400, height=300)
    view.addModel(mblock, 'mol')
    view.setStyle({'stick': {}, 'sphere': {'scale': 0.3}})
    view.zoomTo()
    st.write(f"**{title}**")
    showmol(view, height=300, width=400)

# 3. مدخلات المستخدم
compound_name = st.text_input("Enter Structure Name (e.g., Tartaric acid, Threonine, or 2-butene):", "")

if st.button("Analyze & Visualize Isomers"):
    if not compound_name:
        st.warning("Please enter a compound name first.")
    else:
        try:
            # --- 🔵 تبديل سطر البحث القديم بالبحث الذكي الجديد ---
            smiles = get_smiles_smart(compound_name)
            
            if not smiles:
                st.error(f"❌ No compound found for: {compound_name}")
            else:
                mol = Chem.MolFromSmiles(smiles)
                # -----------------------------------------------
                
                mol_no_stereo = Chem.Mol(mol)
                for bond in mol_no_stereo.GetBonds():
                    bond.SetStereo(Chem.BondStereo.STEREONONE)
                for atom in mol_no_stereo.GetAtoms():
                    atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
                
                isomers = list(EnumerateStereoisomers(mol_no_stereo))
                
                # --- القسم 1: تحليل العلاقات ---
                st.subheader("1. Isomeric Relationships")
                if len(isomers) > 1:
                    st.info("💡 Relationships Analysis:")
                    for i in range(len(isomers)):
                        for j in range(i + 1, len(isomers)):
                            st.write(f"• **Isomer {i+1}** & **Isomer {j+1}**: Stereoisomeric relationship detected.")
                else:
                    st.info("The compound is Achiral (No stereoisomers).")

                # --- القسم 2: العرض الثنائي الأبعاد 2D ---
                st.subheader("2. 2D Structure Grid")
                labels = []
                for i, iso in enumerate(isomers):
                    Chem.AssignStereochemistry(iso, force=True, cleanIt=True)
                    stereo_info = []
                    for bond in iso.GetBonds():
                        stereo = bond.GetStereo()
                        if stereo == Chem.BondStereo.STEREOE: stereo_info.append("E")
                        elif stereo == Chem.BondStereo.STEREOZ: stereo_info.append("Z")
                    
                    centers = Chem.FindMolChiralCenters(iso)
                    for c in centers: stereo_info.append(f"{c[1]}")
                    
                    label = f"Isomer {i+1}: " + (", ".join(stereo_info) if stereo_info else "Achiral")
                    labels.append(label)

                img = Draw.MolsToGridImage(isomers, molsPerRow=3, subImgSize=(300, 300), legends=labels)
                st.image(img, use_container_width=True)

                # --- القسم 3: العرض الثلاثي الأبعاد 3D ---
                st.subheader("3. Interactive 3D Visualization")
                cols = st.columns(3)
                for i, iso in enumerate(isomers):
                    with cols[i % 3]:
                        render_3d(iso, labels[i])

        except Exception as e:
            st.error(f"Error: {e}")

st.markdown("---")
st.caption("Advanced Mode: 3D Rendering & Stereochemical Mapping Active.")
