import streamlit as st
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
from stmol import showmol
import py3Dmol
import numpy as np

# 1. إعدادات الصفحة والتصميم (CSS) لجعل الواجهة مطابقة للصورة
st.set_page_config(page_title="Chemical Isomer Analysis", layout="wide")

st.markdown("""
<style>
    /* تنسيق الحاوية الخاصة بالصور (الخلفية الرمادية والحدود) */
    .isomer-card {
        background-color: #f8f9fa;
        border: 1px solid #e0e0e0;
        border-radius: 8px;
        padding: 10px;
        text-align: center;
        margin-bottom: 20px;
    }
    .stImage {
        border-radius: 5px;
        background-color: white; /* لجعل الرسم نفسه على أبيض داخل الإطار الرمادي */
    }
</style>
""", unsafe_allow_html=True)

# 2. النوت العلمية
st.markdown("""
<div style="background-color: #fdf2f2; padding: 15px; border-radius: 10px; border-left: 5px solid #800000; margin-bottom: 20px;">
    <strong style="color: #800000; font-size: 1.2em;">Stereoisomerism Reference Guide:</strong>
    <ul style="list-style-type: none; padding-left: 0; margin-top: 10px; color: black;">
        <li>1. <b>Cis / Trans:</b> Identical groups on same/opposite sides.</li>
        <li>2. <b>E / Z:</b> CIP System priorities.</li>
        <li>3. <b>R / S:</b> Optical configuration.</li>
    </ul>
</div>
""", unsafe_allow_html=True)

# دالة الرسم المعدلة لتطابق الستايل المطلوب
def render_smart_2d(mol):
    # تحسين الجزيء للرسم
    m = Chem.Mol(mol)
    Chem.AssignStereochemistry(m, force=True, cleanIt=True)
    AllChem.Compute2DCoords(m)
    
    # إعدادات الرسم لجعل الخطوط واضحة والـ Wedges مطابقة للصورة
    d_opts = Draw.MolDrawOptions()
    d_opts.addStereoAnnotation = True # إظهار R/S تحت الذرة
    d_opts.bondLineWidth = 2.5
    d_opts.minFontSize = 20
    d_opts.annotationFontScale = 0.8
    d_opts.clearBackground = True # خلفية شفافة لدمجها مع CSS
    
    # توليد الصورة
    img = Draw.MolToImage(m, size=(400, 400), options=d_opts)
    return img

# المدخلات
name = st.text_input("Enter Structure Name:", "2-butanol")

if st.button("Analyze & Visualize"):
    try:
        results = pcp.get_compounds(name, 'name')
        if results:
            smiles = results[0].isomeric_smiles
            base_mol = Chem.MolFromSmiles(smiles)
            
            opts = StereoEnumerationOptions(tryEmbedding=True, onlyUnassigned=False)
            isomers = list(EnumerateStereoisomers(base_mol, options=opts))
            
            st.markdown("### 2. 2D Structure Grid")
            cols = st.columns(len(isomers))
            
            for i, iso in enumerate(isomers):
                with cols[i]:
                    # حساب الـ R/S لعرضه في العنوان كما بالصورة
                    Chem.AssignStereochemistry(iso, force=True)
                    centers = Chem.FindMolChiralCenters(iso)
                    label = centers[0][1] if centers else ""
                    
                    # وضع الصورة داخل حاوية (Card) رمادية
                    st.markdown(f'<div class="isomer-card">', unsafe_allow_html=True)
                    st.image(render_smart_2d(iso), use_container_width=True)
                    st.markdown(f"<b>Isomer {i+1}: {label}</b>", unsafe_allow_html=True)
                    st.markdown('</div>', unsafe_allow_html=True)

            st.markdown("---")
            st.markdown("### 3. Interactive 3D Visualization")
            cols_3d = st.columns(len(isomers))
            for i, iso in enumerate(isomers):
                with cols_3d[i]:
                    m3d = Chem.AddHs(iso)
                    AllChem.EmbedMolecule(m3d, AllChem.ETKDG())
                    mblock = Chem.MolToMolBlock(m3d)
                    view = py3Dmol.view(width=300, height=300)
                    view.addModel(mblock, 'mol')
                    view.setStyle({'stick': {'radius': 0.2}, 'sphere': {'scale': 0.3}})
                    view.setBackgroundColor('#f8f9fa') # نفس لون خلفية الكارد
                    view.zoomTo()
                    showmol(view)
        else:
            st.error("Compound not found.")
    except Exception as e:
        st.error(f"Error: {e}")
