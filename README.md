# drugresistance
Command 
python run_pipeline.py \
  --r1 /home/muthukumarb/drugresistance/ERR12307390_1.fastq.gz \
  --r2 /home/muthukumarb/drugresistance/ERR12307390_2.fastq.gz \
  --ref /home/muthukumarb/ml/reference/h37rv.fa \
  --sample ERR12307390 \
  --threads 12 \
  --call-variants \
  --normalize \
  --dr-catalog tbdb.compat.json \
  --min-alt-af 0.10 \
  --min-dp 10 \
  --logging DEBUG
                                                                                   
                                                                                        
                                                                                        
######## LineageXpress_new ########
```                                                                                        
python lineageexpress_new_version4.py --sample_list samples.txt --ref_genome /home/muthukumarb/lineageXpress/data/h37rv.fa   --snp_file /home/muthukumarb/lineageXpress/data/lineage_snp_updated_au13.tsv   --bed_file /home/muthukumarb/lineageXpress/data/targeted_modified_regions_au13.bed --threads 12 --mutations_csv mutations.csv --snpeff_cmd "snpEff" --snpeff_db Mycobacterium_tuberculosis_h37rv --output_dir ERR12307390
```

```markdown
![DrugResistance banner](./assets/banner.svg)

# DrugResistance

[![License: MIT](https://img.shields.io/badge/license-MIT-007ec6.svg)](LICENSE)
[![Python](https://img.shields.io/badge/python-3.8%2B-3776AB.svg)](#)
[![Status](https://img.shields.io/badge/status-Research-orange.svg)](#)

A concise toolkit for analyzing drug resistance data and models — now with a cleaner, more readable README and a high-contrast banner for improved visual hierarchy.

## Why this look?
This update focuses on:
- A high-contrast banner to immediately orient visitors.
- Clean layout and readable headings.
- Accessible color palette with good contrast for headings and badges.

Color palette used:
- Primary blue: #0B57A4
- Accent teal: #1FA2A6
- Soft violet: #6A5ACD
- Background gradient (banner): #0B57A4 → #1FA2A6

## Features
- Data parsing and preprocessing utilities
- Resistance metric calculations
- Plotting helpers and visualizations
- Example notebooks demonstrating workflows

## Quickstart

1. Clone the repo:
   ```bash
   git clone https://github.com/muthubio/drugresistance.git
   cd drugresistance
   ```

2. Create a venv and install dependencies:
   ```bash
   python -m venv .venv
   source .venv/bin/activate
   pip install -r requirements.txt
   ```

3. Run an example:
   ```bash
   python examples/run_analysis.py
   ```

## Usage examples
- See the `examples/` folder for demonstration scripts and notebooks.
- `scripts/` contains utility scripts for preprocessing and plotting.

## Contributing
Contributions are welcome. For small changes (typos, docs), open a PR against the default branch or against `improve-readme-color` if you want the visual changes isolated.

When contributing:
- Follow the code style in the repo.
- Add tests for new functionality where applicable.
- Keep PRs focused and document rationale in the description.

## Accessibility notes
- Colors chosen to meet contrast recommendations for headings and badges.
- If you prefer a different palette or lighter/darker banner, I can provide variants.

## License
This project is licensed under the MIT License — see the LICENSE file for details.

## Contact
Author: muthubio — https://github.com/muthubio
```
