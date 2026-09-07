from pipeline_lib import load_config, get_unique_proteins, save_json, write_tsv

def main() -> None:
    cfg = load_config(); out = cfg["output_dir"]
    proteins = get_unique_proteins(cfg)
    write_tsv(proteins, out / "protein_unique_table.tsv")
    save_json({"n_unique_proteins": int(len(proteins)), "deduplication_key": "canonical UniProt accession",
               "columns": list(proteins.columns), "target_class_source": "formal primary membrane role"},
              out / "input_recognition.json")
    print(f"Prepared {len(proteins)} unique canonical proteins")

if __name__ == "__main__": main()
