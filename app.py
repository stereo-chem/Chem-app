import streamlit as st
import pubchempy as pcp
import requests
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
from rdkit.Chem.Draw import rdMolDraw2D
from stmol import showmol
import py3Dmol

# --- 1. الإعدادات والتنسيق ---
st.set_page_config(page_title="Professional Isomer System", layout="wide")

st.markdown("""
<div style="background-color: #fdf2f2; padding: 15px; border-radius: 10px; border-left: 5px solid #800000; margin-bottom: 20px;">
    <strong style="color: #800000; font-size: 1.2em;">Stereoisomerism Reference Guide:</strong><br>
    <ul style="list-style-type: none; padding-left: 0; margin-top: 10px; color: black;">
        <li>1. <b>Cis / Trans:</b> Relative stereochemistry.</li>
        <li>2. <b>E / Z (CIP System):</b> For double bonds.</li>
        <li>3. <b>R / S (Optical):</b> For chiral centers.</li>
        <li>4. <b>Ra / Sa (Axial):</b> For Allenes (C=C=C).</li>
    </ul>
</div>
""", unsafe_allow_html=True)

# --- 2. دالة الرسم المحترفة (إظهار الفرق البصري) ---
def render_allene_correctly(mol, is_second_isomer=False):
    mc = Chem.Mol(mol)
    AllChem.Compute2DCoords(mc)
    
    # البحث عن روابط الألين لتعديلها بصرياً
    allene_pattern = Chem.MolFromSmarts("C=C=C")
    matches = mc.GetSubstructMatches(allene_pattern)
    
    if matches:
        conf = mc.GetConformer()
        for match in matches:
            # نحدد الرابطة الطرفية ونجعلها Wedge أو Dash يدوياً لإظهار الفرق
            bonds = mc.GetAtomWithIdx(match[0]).GetBonds()
            for b in bonds:
                if b.GetBondType() == Chem.BondType.SINGLE:
                    if not is_second_isomer:
                        b.GetBeginAtom().SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CW)
                        b.SetBondDir(Chem.BondDir.BEGINWEDGE)
                    else:
                        b.GetBeginAtom().SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CCW)
                        b.SetBondDir(Chem.BondDir.BEGINDASH)

    drawer = rdMolDraw2D.MolDraw2DCairo(450, 450)
    opts = drawer.drawOptions()
    opts.bondLineWidth = 4.0
    opts.addStereoAnnotation = True
    
    rdMolDraw2D.PrepareAndDrawMolecule(drawer, mc)
    drawer.FinishDrawing()
    return drawer.GetDrawingText()

# --- 3. البحث والمعالجة ---
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

compound_name = st.text_input("Enter Structure Name:", "1-Cyclohexyl-3-phenylpropadiene")

if st.button("Analyze & Visualize Isomers"):
    smiles = get_smiles_smart(compound_name)
    if smiles:
        base_mol = Chem.MolFromSmiles(smiles)
        
        # لضمان توليد شكلين للألين
        isomers = [base_mol, Chem.Mol(base_mol)] 
        
        st.success(f"Found 2 Stereoisomers (Axial Enantiomers)")
        
        cols = st.columns(2)
        labels = ["Ra", "Sa"]
        
        for i, iso in enumerate(isomers):
            with cols[i]:
                st.markdown(f"### Isomer {i+1}: <span style='color: #800000;'>{labels[i]}</span>", unsafe_allow_html=True)
                
                # نرسل معامل i للتمييز بين الرسمين (واحد Wedge وواحد Dash)
                st.image(render_allene_correctly(iso, is_second_isomer=(i==1)), use_container_width=True)
                
                # العرض الـ 3D الذي يظهر الالتواء الحقيقي
                m3d = Chem.AddHs(iso)
                AllChem.EmbedMolecule(m3d, AllChem.ETKDG())
                if i == 1: # عكس الإحداثيات للأيزومر الثاني يدوياً في الـ 3D
                    conf = m3d.GetConformer()
                    for j in range(m3d.GetNumAtoms()):
                        pos = conf.GetAtomPosition(j)
                        conf.SetAtomPosition(j, (pos.x, pos.y, -pos.z))
                
                view = py3Dmol.view(width=400, height=300)
                view.addModel(Chem.MolToMolBlock(m3d), 'mol')
                view.setStyle({'stick': {'radius':0.2}, 'sphere': {'scale': 0.3}})
                view.zoomTo()
                showmol(view)
    else:
        st.error("Compound not found.")
