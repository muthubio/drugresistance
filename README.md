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
