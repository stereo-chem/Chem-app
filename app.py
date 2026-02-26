import streamlit as st
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import AllChem, EnumerateStereoisomers
from stmol import showmol
import py3Dmol

st.set_page_config(page_title="Chemical Structure Viewer", layout="wide")
st.title("🧪 مستعرض وتحليل المركبات الكيميائية")

# إدخال الاسم
name = st.text_input("أدخل اسم المركب (مثلاً: 1-Phenyl-2-methylallene):", "1-Phenyl-2-methylallene")

if st.button("عرض وتحليل المركب"):
    try:
        # 1. جلب البيانات من PubChem
        results = pcp.get_compounds(name, 'name')
        
        if not results:
            st.error(f"❌ لم يتم العثور على مركب باسم '{name}' في قاعدة البيانات. تأكد من كتابة الاسم بشكل صحيح (الإنجليزي).")
        else:
            comp = results[0]
            smiles = comp.isomeric_smiles
            st.success(f"✅ تم العثور على المركب: {comp.iupac_name if comp.iupac_name else name}")
            
            # تحويل SMILES إلى جزيء
            base_mol = Chem.MolFromSmiles(smiles)
            base_mol = Chem.AddHs(base_mol)
            
            # 2. توليد الأيزومرات (حتى لو واحد فقط)
            opts = EnumerateStereoisomers.StereoEnumerationOptions(tryEmbedding=True, onlyUnassigned=False)
            isomers = list(EnumerateStereoisomers.EnumerateStereoisomers(base_mol, options=opts))
            
            num_isomers = len(isomers)
            if num_isomers > 1:
                st.info(f"💡 وجدنا {num_isomers} أيزومرات فراغية مختلفة لهذا المركب.")
            else:
                st.warning("⚠️ هذا المركب لا يمتلك أيزومرات فراغية (أو له شكل واحد فقط ثابت).")

            # 3. عرض النتائج
            for i, iso in enumerate(isomers):
                # تحضير الشكل ثلاثي الأبعاد
                iso = Chem.AddHs(iso)
                AllChem.EmbedMolecule(iso, AllChem.ETKDG())
                AllChem.MMFFOptimizeMolecule(iso)
                
                label = f"الأيزومر الوحيد" if num_isomers == 1 else f"أيزومر رقم {i+1}"
                st.subheader(label)
                
                col1, col2 = st.columns([1, 2])
                
                with col1:
                    st.write("**البيانات الكيميائية:**")
                    Chem.AssignStereochemistry(iso, force=True, cleanIt=True)
                    
                    # البحث عن مراكز R/S
                    chiral_centers = Chem.FindMolChiralCenters(iso, includeUnassigned=True)
                    if chiral_centers:
                        for atom_idx, type_label in chiral_centers:
                            st.info(f"الذرة {iso.GetAtomWithIdx(atom_idx).GetSymbol()}{atom_idx} نوعها: **{type_label}**")
                    else:
                        st.write("• لا توجد مراكز كيرالية تقليدية (R/S).")
                    
                    # البحث عن روابط E/Z
                    found_ez = False
                    for bond in iso.GetBonds():
                        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                            stereo = bond.GetStereo()
                            if stereo in [Chem.rdchem.BondStereo.STEREOE, Chem.rdchem.BondStereo.STEREOZ]:
                                found_ez = True
                                type_ez = "E" if stereo == Chem.rdchem.BondStereo.STEREOE else "Z"
                                st.warning(f"• رابطة ثنائية نوع: **{type_ez}**")
                    
                    if not chiral_centers and not found_ez:
                        st.write("• هذا الهيكل متماثل فراغياً (Achiral).")

                with col2:
                    # العرض الثلاثي الأبعاد
                    view = py3Dmol.view(width=600, height=400)
                    mol_block = Chem.MolToMolBlock(iso)
                    view.addModel(mol_block, 'mol')
                    view.setStyle({'stick': {'radius': 0.2}, 'sphere': {'scale': 0.3}})
                    view.setBackgroundColor('#ffffff')
                    view.zoomTo()
                    showmol(view)
                
                st.divider()

    except Exception as e:
        st.error(f"حدث خطأ أثناء المعالجة: {e}")

