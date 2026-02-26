import streamlit as st
import pubchempy as pcp
import requests
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
from stmol import showmol
import py3Dmol

# 1. إعدادات الصفحة
st.set_page_config(page_title="Advanced Chemical Isomer Analysis", layout="wide")

# دالة البحث الذكي
def get_smiles_smart(name):
    try:
        opsin_url = f"https://opsin.ch.cam.ac.uk/opsin/{name}.json"
        opsin_res = requests.get(opsin_url)
        if opsin_res.status_code == 200:
            return opsin_res.json()['smiles']
    except: pass
    try:
        pcp_res = pcp.get_compounds(name, 'name')
        if pcp_res:
            return pcp_res[0].isomeric_smiles
    except: pass
    return None

# 2. تصميم الواجهة
st.markdown("""
<style>
    .stApp { background-color: white; color: black; }
</style>
<h2 style='color: #800000; font-family: serif; border-bottom: 2px solid #dcdde1;'>Chemical Isomer Analysis System 2.0</h2>
""", unsafe_allow_html=True)

# دالة عرض 3D
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
compound_name = st.text_input("Enter Structure Name:", "1-Cyclohexyl-3-phenylpropadiene")

if st.button("Analyze & Visualize Isomers"):
    if not compound_name:
        st.warning("Please enter a name.")
    else:
        try:
            smiles = get_smiles_smart(compound_name)
            
            if not smiles:
                st.error(f"❌ No compound found.")
            else:
                mol = Chem.MolFromSmiles(smiles)
                
                # --- 🔴 إضافة معالجة الألين (هنا التعديل المطلوب) ---
                allene_pattern = Chem.MolFromSmarts("C=C=C")
                if mol.HasSubstructMatch(allene_pattern):
                    for match in mol.GetSubstructMatches(allene_pattern):
                        # تعيين علامة كيرالية لذرات أطراف الألين لتوليد الأيزومرات
                        mol.GetAtomWithIdx(match[0]).SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CW)
                # --------------------------------------------------

                mol_no_stereo = Chem.Mol(mol)
                # مسح أي معلومات فراغية مسبقة لضمان توليد الكل
                for bond in mol_no_stereo.GetBonds():
                    bond.SetStereo(Chem.BondStereo.STEREONONE)
                for atom in mol_no_stereo.GetAtoms():
                    atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
                
                # تفعيل خيارات التوليد الشامل
                opts = StereoEnumerationOptions(tryEmbedding=True, onlyUnassigned=False)
                isomers = list(EnumerateStereoisomers(mol, options=opts))
                
                # عرض النتائج
                st.subheader(f"Found {len(isomers)} Stereoisomers")
                
                labels = []
                for i, iso in enumerate(isomers):
                    Chem.AssignStereochemistry(iso, force=True, cleanIt=True)
                    labels.append(f"Isomer {i+1}")

                # عرض 2D Grid
                img = Draw.MolsToGridImage(isomers, molsPerRow=3, subImgSize=(300, 300), legends=labels)
                st.image(img, use_container_width=True)

                # عرض 3D Columns
                cols = st.columns(3)
                for i, iso in enumerate(isomers):
                    with cols[i % 3]:
                        render_3d(iso, labels[i])

        except Exception as e:
            st.error(f"Error: {e}")
