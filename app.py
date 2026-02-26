import streamlit as st
import pubchempy as pcp
import requests
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
from rdkit.Chem.Draw import rdMolDraw2D
from stmol import showmol
import py3Dmol

# --- 1. إعدادات الصفحة والستايل ---
st.set_page_config(page_title="Professional Isomer System", layout="wide")

st.markdown("""
<style>
    .stApp { background-color: white; color: black; }
</style>
<div style="background-color: #fdf2f2; padding: 15px; border-radius: 10px; border-left: 5px solid #800000; margin-bottom: 20px;">
    <strong style="color: #800000; font-size: 1.2em;">Stereoisomerism Reference Guide:</strong><br>
    <ul style="list-style-type: none; padding-left: 0; margin-top: 10px; color: black;">
        <li>1. <b>Cis / Trans:</b> Identical groups on same/opposite sides.</li>
        <li>2. <b>E / Z (CIP System):</b> <b>Z (Zusammen)</b> together, <b>E (Entgegen)</b> opposite.</li>
        <li>3. <b>R / S (Optical):</b> Absolute configuration of chiral centers.</li>
        <li>4. <b>Ra / Sa (Axial):</b> Stereochemistry of Allenes (C=C=C).</li>
    </ul>
</div>
""", unsafe_allow_html=True)

st.markdown("<h2 style='color: #800000; font-family: serif; border-bottom: 2px solid #dcdde1;'>Chemical Isomer Analysis System 2.0</h2>", unsafe_allow_html=True)

# --- 2. الدوال المساعدة ---
def get_smiles_smart(name):
    try:
        opsin_url = f"https://opsin.ch.cam.ac.uk/opsin/{name}.json"
        res = requests.get(opsin_url)
        if res.status_code == 200: return res.json()['smiles']
    except: pass
    try:
        pcp_res = pcp.get_compounds(name, 'name')
        if pcp_res: return pcp_res[0].isomeric_smiles
    except: pass
    return None

def render_pro_2d(mol):
    # محرك رسم عالي الجودة لضمان وضوح الـ Wedges
    drawer = rdMolDraw2D.MolDraw2DCairo(450, 450)
    opts = drawer.drawOptions()
    opts.bondLineWidth = 3.5  # روابط سميكة واضحة
    opts.addStereoAnnotation = True
    opts.addAtomIndices = False
    
    mc = Chem.Mol(mol)
    AllChem.Compute2DCoords(mc)
    Chem.WedgeMolBonds(mc, mc.GetConformer()) # تفعيل الروابط المجسمة
    drawer.DrawMolecule(mc)
    drawer.FinishDrawing()
    return drawer.GetDrawingText()

def render_3d(mol, title):
    m3d = Chem.AddHs(mol)
    AllChem.EmbedMolecule(m3d, AllChem.ETKDG())
    mblock = Chem.MolToMolBlock(m3d)
    view = py3Dmol.view(width=400, height=300)
    view.addModel(mblock, 'mol')
    view.setStyle({'stick': {'radius': 0.15}, 'sphere': {'scale': 0.25}})
    view.zoomTo()
    st.write(f"**{title}**")
    showmol(view, height=300, width=400)

# --- 3. المعالجة والعرض ---
compound_name = st.text_input("Enter Structure Name:", "1-Cyclohexyl-3-phenylpropadiene")

if st.button("Analyze & Visualize Isomers"):
    smiles = get_smiles_smart(compound_name)
    if smiles:
        mol = Chem.MolFromSmiles(smiles)
        
        # كود اكتشاف الألين وإجبار البرنامج على إنتاج أيزومرات
        allene_pattern = Chem.MolFromSmarts("C=C=C")
        if mol.HasSubstructMatch(allene_pattern):
            for match in mol.GetSubstructMatches(allene_pattern):
                mol.GetAtomWithIdx(match[0]).SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CW)

        # توليد الأيزومرات الشاملة
        opts = StereoEnumerationOptions(tryEmbedding=True, onlyUnassigned=False)
        isomers = list(EnumerateStereoisomers(mol, options=opts))
        
        st.subheader(f"Found {len(isomers)} Stereoisomer(s)")
        
        # القسم 1: العلاقات الأيزومرية
        st.subheader("1. Isomeric Relationships")
        if len(isomers) > 1:
            st.info("💡 Relationships Analysis: Stereoisomeric relationship detected.")
        else:
            st.info("The compound is achiral or only one isomer was identified.")

        st.write("---")

        # القسم 2 و 3: العرض في أعمدة منظمة
        cols = st.columns(len(isomers))
        for i, iso in enumerate(isomers):
            with cols[i]:
                Chem.AssignStereochemistry(iso, force=True, cleanIt=True)
                centers = Chem.FindMolChiralCenters(iso, includeUnassigned=True)
                label = f"Isomer {i+1}"
                if centers: label += f" ({centers[0][1]})"
                
                st.markdown(f"### {label}")
                st.image(render_pro_2d(iso), use_container_width=True) # عرض 2D المحسن
                render_3d(iso, "Interactive 3D View") # عرض 3D
                
    else:
        st.error("Compound not found.")
