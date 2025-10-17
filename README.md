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
                                                                                   
                                                                                        
                                                                                        
                                                                                        ###### LineageXpress_new ####################
                                                                                        
python lineageexpress_new_version4.py --sample_list samples.txt --ref_genome /home/muthukumarb/lineageXpress/data/h37rv.fa   --snp_file /home/muthukumarb/lineageXpress/data/lineage_snp_updated_au13.tsv   --bed_file /home/muthukumarb/lineageXpress/data/targeted_modified_regions_au13.bed --threads 12 --mutations_csv mutations.csv --snpeff_cmd "snpEff" --snpeff_db Mycobacterium_tuberculosis_h37rv --output_dir ERR12307390

