import streamlit as st
import requests
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import AllChem, EnumerateStereoisomers, Draw
from stmol import showmol
import py3Dmol

st.set_page_config(page_title="Ultimate Molecule Explorer", layout="wide")
st.title("🧪 Smart Chemical Isomer Analyzer")
st.markdown("Search by **IUPAC**, **Common Name**, or **SMILES**")

def get_smiles_smart(name):
    # الخطوة 1: محاولة استخدام OPSIN للأسماء المنهجية
    try:
        opsin_url = f"https://opsin.ch.cam.ac.uk/opsin/{name}.json"
        opsin_res = requests.get(opsin_url)
        if opsin_res.status_code == 200:
            return opsin_res.json()['smiles']
    except:
        pass

    # الخطوة 2: لو فشل، نبحث في PubChem للأسماء الشائعة
    try:
        pcp_res = pcp.get_compounds(name, 'name')
        if pcp_res:
            return pcp_res[0].isomeric_smiles
    except:
        pass
    
    return None

name_input = st.text_input("Enter Molecule Name (Common or IUPAC):", "1,3-diphenylallene")

if st.button("Generate & Analyze"):
    with st.spinner('Searching for molecule...'):
        smiles = get_smiles_smart(name_input)
    
    if not smiles:
        st.error("Couldn't find this molecule. Please check the spelling or try a more formal name.")
    else:
        # تحويل الـ SMILES لجزيء
        base_mol = Chem.MolFromSmiles(smiles)
        base_mol = Chem.AddHs(base_mol)
        
        # توليد الأيزومرات الشاملة
        opts = EnumerateStereoisomers.StereoEnumerationOptions(tryEmbedding=True, onlyUnassigned=False)
        isomers = list(EnumerateStereoisomers.EnumerateStereoisomers(base_mol, options=opts))
        
        st.success(f"Found structure for '{name_input}' with {len(isomers)} possible isomer(s).")
        
        # عرض العلاقات الأيزومرية في البداية
        st.info(f"💡 Relationship: { 'Stereoisomers detected' if len(isomers) > 1 else 'Single structure found' }")

        # إنشاء الشبكة (Grid) كما طلبتِ في الصورة الأولى
        cols = st.columns(len(isomers))
        
        for i, iso in enumerate(isomers):
            with cols[i]:
                # تحديد الكيمياء الفراغية
                Chem.AssignStereochemistry(iso, force=True, cleanIt=True)
                centers = Chem.FindMolChiralCenters(iso, includeUnassigned=True)
                
                label = ""
                if centers:
                    label = f"({centers[0][1]})" # R or S
                
                st.subheader(f"Isomer {i+1} {label}")
                
                # عرض 2D
                img = Draw.MolToImage(iso, size=(300, 300))
                st.image(img, use_container_width=True)
                
                # عرض 3D
                iso_3d = Chem.AddHs(iso)
                AllChem.EmbedMolecule(iso_3d, AllChem.ETKDG())
                AllChem.MMFFOptimizeMolecule(iso_3d)
                
                view = py3Dmol.view(width=400, height=300)
                view.addModel(Chem.MolToMolBlock(iso_3d), 'mol')
                view.setStyle({'stick': {'radius': 0.15}, 'sphere': {'scale': 0.25}})
                view.zoomTo()
                showmol(view, height=300, width=400)
