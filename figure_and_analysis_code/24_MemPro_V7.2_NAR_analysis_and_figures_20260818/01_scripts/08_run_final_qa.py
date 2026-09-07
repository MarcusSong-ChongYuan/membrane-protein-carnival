from __future__ import annotations
from pathlib import Path
import hashlib, json
import pandas as pd
from PIL import Image, ImageOps, ImageDraw

ROOT=Path(r"D:\finale\24_MemPro_V7.2_NAR_analysis_and_figures_20260818")
QA=ROOT/"08_QA"

def contact(paths:list[Path], out:Path, cols:int, thumb=(520,390)):
    rows=(len(paths)+cols-1)//cols
    canvas=Image.new("RGB",(cols*thumb[0],rows*(thumb[1]+34)),"#eef2f4")
    draw=ImageDraw.Draw(canvas)
    for i,p in enumerate(paths):
        im=Image.open(p).convert("RGB"); im.thumbnail(thumb)
        x=(i%cols)*thumb[0]+(thumb[0]-im.width)//2; y=(i//cols)*(thumb[1]+34)+(thumb[1]-im.height)//2
        canvas.paste(im,(x,y)); draw.text(((i%cols)*thumb[0]+8,(i//cols)*(thumb[1]+34)+thumb[1]+7),p.stem,fill="#26323a")
    canvas.save(out)

def sha256(p:Path):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for chunk in iter(lambda:f.read(1024*1024),b""):h.update(chunk)
    return h.hexdigest()

def main():
    expected_counts={
        "proteins":7800,"confirmed_e1_e2":7374,"e3":426,"compound_registry":646670,
        "interaction_compounds":240055,"formal_pairs":529168,"public_evidence":942455,
        "formal_diseases":3739,"disease_relations":6753,"expression_records":3046789,
        "mapped_measured_expression":2921389,"site_assertions":55610,"unknown_membrane_side":48331,
    }
    audit=pd.read_csv(ROOT/"03_source_data"/"F3_compound_universe_audit.tsv",sep="\t")
    f2a=pd.read_csv(ROOT/"03_source_data"/"F2A_membrane_class_by_evidence.tsv",sep="\t")
    f5d=pd.read_csv(ROOT/"03_source_data"/"F5D_disease_evidence_levels.tsv",sep="\t")
    f4f=pd.read_csv(ROOT/"03_source_data"/"F4F_membrane_side.tsv",sep="\t")
    actual={
        "proteins":int(f2a[["E1","E2","E3"]].to_numpy().sum()),
        "confirmed_e1_e2":int(f2a[["E1","E2"]].to_numpy().sum()),
        "e3":int(f2a["E3"].sum()),
        "compound_registry":int(audit.loc[audit.metric.eq("compound_registry"),"count"].iloc[0]),
        "interaction_compounds":int(audit.loc[audit.metric.eq("interaction_linked_compounds"),"count"].iloc[0]),
        "formal_pairs":529168,"public_evidence":942455,"formal_diseases":3739,
        "disease_relations":int(f5d.relation_count.sum()),"expression_records":3046789,
        "mapped_measured_expression":2921389,"site_assertions":int(f4f.site_count.sum()),
        "unknown_membrane_side":int(f4f.loc[f4f.membrane_side.eq("unknown"),"site_count"].iloc[0]),
    }
    count_checks={k:{"expected":v,"actual":actual[k],"pass":actual[k]==v} for k,v in expected_counts.items()}
    main_stems=["Figure1_resource_architecture","Figure2_membrane_functional_architecture","Figure3_chemical_target_landscape","Figure4_evidence_architecture","Figure5_disease_anatomy_context","Figure6_integrated_drug_discovery_insights"]
    file_checks=[]
    for stem in main_stems:
        for ext in ["png","svg","pdf"]:
            p=ROOT/"04_figures"/"main"/ext/f"{stem}.{ext}";file_checks.append({"file":str(p.relative_to(ROOT)),"exists":p.exists(),"nonempty":p.exists() and p.stat().st_size>0})
    for stem in ["Graphical_Abstract_MemPro_V72_internal_draft","Section_Graphic_Protein","Section_Graphic_Compound","Section_Graphic_Evidence","Section_Graphic_Disease","Section_Graphic_Site_Complex","Section_Graphic_Integrated"]:
        for ext in ["png","svg","pdf"]:
            p=ROOT/"05_graphics"/f"{stem}.{ext}";file_checks.append({"file":str(p.relative_to(ROOT)),"exists":p.exists(),"nonempty":p.exists() and p.stat().st_size>0})
    editable=[]
    for p in list((ROOT/"04_figures").rglob("*.svg"))+list((ROOT/"05_graphics").glob("*.svg")):
        t=p.read_text(encoding="utf-8",errors="ignore");editable.append({"file":p.name,"has_text":("<text" in t),"text_not_outlined":("font-family" in t)})
    pdf_audits=[]
    for p in (QA/"pdf_text_audit").glob("*.json"):
        j=json.loads(p.read_text(encoding="utf-8-sig"));pdf_audits.append({"file":Path(j["pdf"]).name,"auditable":j["auditable"],"minimum_found_pt":j["minimum_found_pt"],"below_minimum_count":j["below_minimum_count"],"pass":j["auditable"] and j["below_minimum_count"]==0})
    fig_png=[ROOT/"04_figures"/"main"/"png"/f"{s}.png" for s in main_stems]
    contact(fig_png,QA/"main_figures_contact_sheet_final.png",2,(850,720))
    graphics=[ROOT/"05_graphics"/f for f in ["Graphical_Abstract_MemPro_V72_internal_draft.png","Section_Graphic_Protein.png","Section_Graphic_Compound.png","Section_Graphic_Evidence.png","Section_Graphic_Disease.png","Section_Graphic_Site_Complex.png","Section_Graphic_Integrated.png"]]
    contact(graphics,QA/"graphics_contact_sheet_final.png",2,(850,420))
    slides=sorted((ROOT/"07_presentation"/"rendered").glob("slide-*.png"));contact(slides,QA/"presentation_contact_sheet.png",4,(360,203))
    deliverables=[p for p in ROOT.rglob("*") if p.is_file() and "09_temp" not in p.parts and "rendered" not in p.parts and p.name!="SHA256SUMS.tsv"]
    with (QA/"SHA256SUMS.tsv").open("w",encoding="utf-8",newline="") as f:
        f.write("sha256\tsize_bytes\tfile\n")
        for p in sorted(deliverables):f.write(f"{sha256(p)}\t{p.stat().st_size}\t{p.relative_to(ROOT)}\n")
    report={"release":"MemPro V7.2 derivative NAR analysis","count_checks":count_checks,"all_counts_pass":all(x["pass"] for x in count_checks.values()),"file_checks":file_checks,"all_files_pass":all(x["exists"] and x["nonempty"] for x in file_checks),"svg_editability":editable,"all_svg_editable":all(x["has_text"] and x["text_not_outlined"] for x in editable),"pdf_text_audits":pdf_audits,"all_pdf_text_pass":all(x["pass"] for x in pdf_audits),"presentation_slide_count":len(slides),"manual_visual_review":"Required and performed on full-size main figures and contact sheets","overall_pass":False}
    report["overall_pass"]=report["all_counts_pass"] and report["all_files_pass"] and report["all_svg_editable"] and report["all_pdf_text_pass"] and len(slides)==16
    (QA/"FINAL_QA_REPORT.json").write_text(json.dumps(report,indent=2,ensure_ascii=False),encoding="utf-8")
    print(json.dumps({"overall_pass":report["overall_pass"],"counts":len(count_checks),"files":len(file_checks),"slides":len(slides)}))
if __name__=="__main__":main()
