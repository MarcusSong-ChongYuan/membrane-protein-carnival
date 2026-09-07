# Leiden structural-map trial

This is a separate, all-protein graph-community analysis, not a replacement for the ESM-only HDBSCAN baseline.

The k-nearest-neighbour graph combines equally weighted blocks: ESM sequence features, membrane topology and structural-family features. Primary membrane role, GO molecular function/biological process, specialist classification, ligand evidence and diseases are excluded from graph construction. They are used after community discovery for interpretation only.

Leiden assigns every protein to a graph community. Two granularities are reported: macro modules (target 8–14) and finer communities (target 30–60). The requested scale only chooses among resolution scan candidates; it does not use functional labels.
