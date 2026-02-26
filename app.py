import streamlit as st
import pubchempy as pcp
import requests
from rdkit import Chem
from rdkit.Chem import Draw, AllChem, rdDepictor
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
from rdkit.Chem.Draw import rdMolDraw2D
from stmol import showmol
import py3Dmol

# --- إعدادات الصفحة ---
st.set_page_config(page_title="Allene Isomer Pro", layout="wide")

st.markdown("""
<div style="background-color: #fdf2f2; padding: 15px; border-radius: 10px; border-left: 5px solid #800000; margin-bottom: 20px;">
    <strong style="color: #800000; font-size: 1.2em;">دليل الأيزومرات الفراغية (Stereoisomers):</strong><br>
    محور الألين (C=C=C) يتطلب وجود مجموعات مختلفة على الأطراف لظهور خاصية الـ Axial Chirality (Ra/Sa).
</div>
""", unsafe_allow_html=True)

# --- الدوال المساعدة ---
def get_smiles_smart(name):
    """محاولة جلب الـ SMILES من أكثر من مصدر"""
    try:
        # محاولة أولى: OPSIN
        res = requests.get(f"https://opsin.ch.cam.ac.uk/opsin/{name}.json", timeout=5)
        if res.status_code == 200: return res.json()['smiles']
    except: pass
    try:
        # محاولة ثانية: PubChem
        pcp_res = pcp.get_compounds(name, 'name')
        if pcp_res: return pcp_res[0].isomeric_smiles
    except: pass
    return None

def render_pro_2d(mol):
    """رسم 2D مع إظهار الـ Wedges (Solid/Hatched) بوضوح"""
    mc = Chem.Mol(mol)
    rdDepictor.Compute2DCoords(mc)
    Chem.AssignStereochemistry(mc, force=True, cleanIt=True)
    # توليد الروابط البارزة والمخفية
    mc = rdMolDraw2D.PrepareMolForDrawing(mc)
    
    drawer = rdMolDraw2D.MolDraw2DCairo(500, 500)
    opts = drawer.drawOptions()
    opts.bondLineWidth = 4.0  # تقوية الخطوط للوضوح
    opts.addStereoAnnotation = True
    opts.explicitMethyl = True
    
    drawer.DrawMolecule(mc)
    drawer.FinishDrawing()
    return drawer.GetDrawingText()

# --- الواجهة الرئيسية ---
compound_name = st.text_input("ادخلي اسم المركب (مثال: 1,3-dichloropropadiene):", "1,3-Dimethyl-3-phenylallene")

if st.button("تحليل ورسم الأيزومرات"):
    with st.spinner('جاري البحث والتحليل...'):
        smiles = get_smiles_smart(compound_name)
        
        if smiles:
            mol = Chem.MolFromSmiles(smiles)
            allene_p = Chem.MolFromSmarts("C=C=C")
            
            if mol.HasSubstructMatch(allene_p):
                # تمييز ذرات الألين للتعامل معها كمركز كيرالي
                for match in mol.GetSubstructMatches(allene_p):
                    mol.GetAtomWithIdx(match[0]).SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CW)
                    mol.GetAtomWithIdx(match[2]).SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CW)

                opts = StereoEnumerationOptions(tryEmbedding=True, onlyUnassigned=False)
                isomers = list(EnumerateStereoisomers(mol, options=opts))
                
                # تأكيد توليد الزوج (Enantiomers) في حالة الألين
                if len(isomers) == 1:
                    iso2 = Chem.Mol(isomers[0])
                    for a in iso2.GetAtoms():
                        if a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED:
                            tag = Chem.ChiralType.CHI_TETRAHEDRAL_CCW if a.GetChiralTag() == Chem.ChiralType.CHI_TETRAHEDRAL_CW else Chem.ChiralType.CHI_TETRAHEDRAL_CW
                            a.SetChiralTag(tag)
                    isomers.append(iso2)

                st.success(f"تم العثور على {len(isomers)} أيزومر")
                
                cols = st.columns(len(isomers))
                for i, iso in enumerate(isomers):
                    with cols[i]:
                        label = "Ra (Right)" if i == 0 else "Sa (Left)"
                        st.markdown(f"### Isomer: <span style='color: #800000;'>{label}</span>", unsafe_allow_html=True)
                        
                        # العرض 2D (Wedges)
                        st.image(render_pro_2d(iso), caption=f"2D Projection ({label})")
                        
                        # العرض 3D
                        m3d = Chem.AddHs(iso)
                        AllChem.EmbedMolecule(m3d, AllChem.ETKDG())
                        
                        view = py3Dmol.view(width=400, height=300)
                        view.addModel(Chem.MolToMolBlock(m3d), 'mol')
                        
                        # تلوين ذرات المحور الكيرالي بالأحمر
                        allene_indices = []
                        for m in iso.GetSubstructMatches(allene_p):
                            allene_indices.extend([m[0], m[2]]) # الذرات الطرفية
                        
                        for idx in allene_indices:
                            view.setStyle({'serial': idx + 1}, {'sphere': {'color': 'red', 'scale': 0.4}})
                        
                        view.setStyle({'not': {'serial': [idx+1 for idx in allene_indices]}}, {'stick': {}})
                        view.zoomTo()
                        showmol(view)
            else:
                st.warning("المركب لا يحتوي على نظام ألين (C=C=C).")
        else:
            st.error("لم يتم العثور على المركب. تأكدي من الاسم أو الاتصال بالإنترنت.")

