from pathlib import Path

p = Path(r"C:\Users\Administrator\MemPro_redraw_v2\build_clean_m1_m8.py")
s = p.read_text(encoding="utf-8")
s = s.replace(
    'anatomy(ax,organ,"detected_proteins","A  Broad tissue RNA coverage spans major organs","YlGnBu","proteins",8)',
    'anatomy(ax,organ,"detected_proteins","YlGnBu","A  Broad tissue RNA coverage spans major organs","proteins",8)',
)
s = s.replace(
    'anatomy(ax,organ,"protein_disease_pairs","A  Disease burden is distributed across organ systems","YlOrRd","pairs",8)',
    'anatomy(ax,organ,"protein_disease_pairs","YlOrRd","A  Disease burden is distributed across organ systems","pairs",8)',
)
s = s.replace('height_ratios=[1.08,1],hspace=.43', 'height_ratios=[1.08,1],hspace=.52')
s = s.replace('ia=ax.inset_axes([0,-.48,1,.40])', 'ia=ax.inset_axes([0,-.40,1,.32])')
old = 'channels={k.replace("_flag","").replace("_"," "):v for k,v in S["disease"].get("channels",{}).items()} if "channels" in S["disease"] else S["disease"]["sources"]'
new = 'chdf=pd.read_csv(Path(r"C:\\Users\\Administrator\\Desktop\\gkg sjf summer intern\\MemPro最终展示\\M1-M8新叙事终版_20260814\\02_plotting_data\\M4_evidence_channel_counts_v711.tsv"),sep="\\t"); channels=dict(zip(chdf["channel"],chdf["pairs"]))'
s = s.replace(old, new)
exec(compile(s, str(p), "exec"), {"__name__": "__main__"})
