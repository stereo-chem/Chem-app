import streamlit as st
import pubchempy as pcp
import requests
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
from rdkit.Chem.Draw import rdMolDraw2D
from stmol import showmol
import py3Dmol

# --- 1. إعدادات الصفحة والتصميم الأساسي ---
st.set_page_config(page_title="Professional Isomer System", layout="wide")

# دليل المراجع العلمي (Reference Guide) بنفس الستايل المطلوب
st.markdown("""
<div style="background-color: #fdf2f2; padding: 15px; border-radius: 10px; border-left: 5px solid #800000; margin-bottom: 20px;">
    <strong style="color: #800000; font-size: 1.2em;">Stereoisomerism Reference Guide:</strong><br>
    <ul style="list-style-type: none; padding-left: 0; margin-top: 10px; color: black;">
        <li>1. <b>Cis / Trans (Relative):</b> Identical groups on same/opposite sides.</li>
        <li>2. <b>E / Z (Absolute - CIP System):</b> High-priority groups together (Z) or opposite (E).</li>
        <li>3. <b>R / S (Optical):</b> Absolute configuration of chiral centers.</li>
        <li>4. <b>Ra / Sa (Axial):</b> Stereochemistry of Allenes (C=C=C).</li>
    </ul>
</div>
""", unsafe_allow_html=True)

st.markdown("<h2 style='color: #800000; font-family: serif; border-bottom: 2px solid #dcdde1;'>Chemical Isomer Analysis System 2.0</h2>", unsafe_allow_html=True)

# --- 2. الدوال المساعدة ---

def get_smiles_smart(name):
    """البحث المزدوج عن المركب"""
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
    """تحسين جودة رسم الـ 2D لتكون واضحة واحترافية"""
    drawer = rdMolDraw2D.MolDraw2DCairo(400, 400)
    opts = drawer.drawOptions()
    opts.bondLineWidth = 3 # روابط سميكة
    opts.addStereoAnnotation = True # إظهار علامات R/S و E/Z
    opts.addAtomIndices = False
    
    # تحضير الجزيء للرسم
    mc = Chem.Mol(mol)
    Chem.AssignStereochemistry(mc, force=True, cleanIt=True)
    drawer.DrawMolecule(mc)
    drawer.FinishDrawing()
    return drawer.GetDrawingText()

def render_3d(mol):
    """عرض 3D تفاعلي"""
    m3d = Chem.AddHs(mol)
    AllChem.EmbedMolecule(m3d, AllChem.ETKDG())
    AllChem.MMFFOptimizeMolecule(m3d)
    mblock = Chem.MolToMolBlock(m3d)
    view = py3Dmol.view(width=400, height=300)
    view.addModel(mblock, 'mol')
    view.setStyle({'stick': {'radius': 0.15}, 'sphere': {'scale': 0.25}})
    view.zoomTo()
    showmol(view, height=300, width=400)

# --- 3. المنطق البرمجي الرئيسي ---

compound_name = st.text_input("Enter Structure Name:", "1-Cyclohexyl-3-phenylpropadiene")

if st.button("Analyze & Visualize Isomers"):
    if not compound_name:
        st.warning("Please enter a name.")
    else:
        try:
            smiles = get_smiles_smart(compound_name)
            if not smiles:
                st.error("❌ Compound not found.")
            else:
                base_mol = Chem.MolFromSmiles(smiles)
                
                # إجبار البرنامج على تفعيل كيرالية الألين (C=C=C)
                allene_pattern = Chem.MolFromSmarts("C=C=C")
                if base_mol.HasSubstructMatch(allene_pattern):
                    for match in base_mol.GetSubstructMatches(allene_pattern):
                        base_mol.GetAtomWithIdx(match[0]).SetChiralTag(
