#!/usr/bin/env python3
from pathlib import Path
src=Path(r"C:\Users\Administrator\integrate_v701_rescues.py").read_text(encoding='utf-8')
src=src.replace('binding_site_coordinate_rescued_v701.tsv.gz','binding_site_coordinate_rescued_v702.tsv.gz')
src=src.replace('binding_site_coordinate_remaining_v701.tsv.gz','binding_site_coordinate_remaining_v702.tsv.gz')
src=src.replace('binding_site_coordinate_verified_v701.tsv.gz','binding_site_coordinate_verified_v702.tsv.gz')
src=src.replace("'rescue_status_v701'","'rescue_status_v702'")
src=src.replace('V7.0.1 deterministic rescue candidate','V7.0.1 deterministic streaming rescue candidate')
exec(compile(src,'integrate_v702_rescues.generated.py','exec'))
