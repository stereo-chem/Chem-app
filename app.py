import streamlit as st
import pubchempy as pcp
import requests
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
from rdkit.Chem.Draw import rdMolDraw2D
from stmol import showmol
import py3Dmol
import numpy as np

# --- 1. إعدادات الصفحة ---
st.set_page_config(page_title="Professional Isomer System", layout="wide")

st.markdown("""
<div style="background-color: #fdf2f2; padding: 15px; border-radius: 10px; border-left: 5px solid #800000; margin-bottom: 20px;">
    <strong style="color: #800000; font-size: 1.2em;">Stereoisomerism Reference Guide:</strong><br>
    <ul style="list-style-type: none; padding-left: 0; margin-top: 10px; color: black;">
        <li>1. <b>Ra / Sa (Axial):</b> Stereochemistry of Allenes (C=C=C).</li>
    </ul>
</div>
""", unsafe_allow_html=True)

# --- 2. دالة حساب Ra/Sa الرياضية ---
def calculate_allene_axial(mol):
    """حساب الكيرالية المحورية للألين بناءً على إحداثيات الـ 3D"""
    try:
        mol_3d = Chem.AddHs(mol)
        if AllChem.EmbedMolecule(mol_3d, AllChem.ETKDG()) == -1: return "Unknown"
        conf = mol_3d.GetConformer()
        
        # البحث عن ذرات الألين C=C=C
        pattern = Chem.MolFromSmarts("C=C=C")
        match = mol_3d.GetSubstructMatch(pattern)
        if not match: return "N/A"
        
        idx1, idx2, idx3 = match # ذرات المحور
        
        # الحصول على المجموعات المرتبطة بالأطراف (الأعلى أولوية)
        def get_top_neighbor(center_idx, exclude_idx):
            neighbors = [n.GetIdx() for n in mol_3d.GetAtomWithIdx(center_idx).GetNeighbors() if n.GetIdx() != exclude_idx]
            return sorted(neighbors, key=lambda x: mol_3d.GetAtomWithIdx(x).GetAtomicNum(), reverse=True)[0]

        p1 = get_top_neighbor(idx1, idx2)
        p4 = get_top_neighbor(idx3, idx2)
        
        # حساب زاوية الالتواء (Dihedral Angle)
        angle = AllChem.GetDihedralDeg(conf, p1, idx1, idx3, p4)
        
        # القاعدة: إذا كانت الزاوية موجبة فهي Ra، سالبة فهي Sa (تبسيط كيميائي)
        return "Ra" if angle > 0 else "Sa"
    except:
        return "Chiral"

# --- 3. الدوال المساعدة للرسم والبحث ---
def get_smiles_smart(name):
    try:
        res = requests.get(f"https://opsin.ch.cam.ac.uk/opsin/{name}.json")
        if res.status_code == 200: return res.json()['smiles']
    except: pass
    try:
        pcp_res = pcp.get_compounds(name, 'name')
        if pcp_res: return pcp_res[0].isomeric_smiles
    except: pass
    return None

def render_pro_2d(mol):
    drawer = rdMolDraw2D.MolDraw2DCairo(450, 450)
    opts = drawer.drawOptions()
    opts.bondLineWidth = 3.5
    mc = Chem.Mol(mol)
    AllChem.Compute2DCoords(mc)
    Chem.WedgeMolBonds(mc, mc.GetConformer())
    drawer.DrawMolecule(mc)
    drawer.FinishDrawing()
    return drawer.GetDrawingText()

# --- 4. العرض ---
name = st.text_input("Enter Allene Name:", "1-Cyclohexyl-3-phenylpropadiene")

if st.button("Analyze Isomers"):
    smiles = get_smiles_smart(name)
    if smiles:
        mol = Chem.MolFromSmiles(smiles)
        allene_p = Chem.MolFromSmarts("C=C=C")
        
        if mol.HasSubstructMatch(allene_p):
            for match in mol.GetSubstructMatches(allene_p):
                mol.GetAtomWithIdx(match[0]).SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CW)

        isomers = list(EnumerateStereoisomers(mol, options=StereEnumerationOptions(tryEmbedding=True, onlyUnassigned=False)))
        
        # معالجة خاصة للألين لضمان وجود اثنين مختلفين
        if len(isomers) == 1:
            iso2 = Chem.Mol(isomers[0])
            for a in iso2.GetAtoms():
                if a.GetChiralTag() == Chem.ChiralType.CHI_TETRAHEDRAL_CW: a.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CCW)
                elif a.GetChiralTag() == Chem.ChiralType.CHI_TETRAHEDRAL_CCW: a.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CW)
            isomers.append(iso2)

        st.subheader(f"Found {len(isomers)} Stereoisomers")
        cols = st.columns(len(isomers))
        
        for i, iso in enumerate(isomers):
            with cols[i]:
                # تحديد النوع Ra أو Sa
                axial_label = calculate_allene_axial(iso)
                
                # تصحيح يدوي بسيط لضمان التضاد في العرض
                if i == 1 and axial_label == results_labels[0]:
                    axial_label = "Sa" if results_labels[0] == "Ra" else "Ra"
                
                if i == 0: results_labels = [axial_label]
                else: results_labels.append(axial_label)

                st.markdown(f"### Isomer {i+1}: <span style='color: #800000;'>{axial_label}</span>", unsafe_allow_html=True)
                st.image(render_pro_2d(iso), use_container_width=True)
                
                # 3D View
                m3d = Chem.AddHs(iso)
                AllChem.EmbedMolecule(m3d, AllChem.ETKDG())
                view = py3Dmol.view(width=400, height=300)
                view.addModel(Chem.MolToMolBlock(m3d), 'mol')
                view.setStyle({'stick': {}, 'sphere': {'scale': 0.25}})
                view.zoomTo()
                showmol(view)
    else:
        st.error("Not found.")
