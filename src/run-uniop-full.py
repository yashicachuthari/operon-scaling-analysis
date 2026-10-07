import os
import csv
import gzip
import glob
import subprocess
import tempfile
from Bio import SeqIO

PROJECT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GBK_DIR = os.path.join(PROJECT, "results", "gbk-annotation")
OUT_DIR = os.path.join(PROJECT, "results", "uniop-all")
UNIOP = os.path.join(PROJECT, "UniOP", "src", "UniOP")
PRODIGAL_DIR = os.environ.get("PRODIGAL_DIR", os.path.expanduser("~/miniforge3/bin"))
CSV_FILE = os.path.join(OUT_DIR, "operon_scaling_results.csv")

os.makedirs(OUT_DIR, exist_ok=True)

gbff_files = sorted(glob.glob(os.path.join(GBK_DIR, "*.gbff.gz")))

# Find genomes already completed in an earlier run
completed = set()

if os.path.exists(CSV_FILE):
    with open(CSV_FILE, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            completed.add(row["genome_file"])

print(f"Found {len(gbff_files)} genome files")
print(f"Already completed: {len(completed)}")
print(f"Remaining: {len(gbff_files) - len(completed)}")

# Create CSV header if needed
if not os.path.exists(CSV_FILE) or os.path.getsize(CSV_FILE) == 0:
    with open(CSV_FILE, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "genome_file",
            "replicon_id",
            "replicon_type",
            "length_bp",
            "operon_count"
        ])

for genome_num, gbff in enumerate(gbff_files, 1):

    genome_name = os.path.basename(gbff)

    if genome_name in completed:
        continue

    print(f"\nGenome {genome_num}/{len(gbff_files)}: {genome_name}")

    genome_rows = []

    try:
        with gzip.open(gbff, "rt") as handle:
            records = list(SeqIO.parse(handle, "genbank"))

        for record in records:

            description = record.description.lower()

            if "plasmid" in description:
                replicon_type = "plasmid"
            else:
                replicon_type = "chromosome"

            length_bp = len(record.seq)

            with tempfile.TemporaryDirectory() as tmp:

                fna = os.path.join(tmp, f"{record.id}.fna")
                result_dir = os.path.join(tmp, "result")
                os.makedirs(result_dir)

                SeqIO.write(record, fna, "fasta")

                cmd = [
                    "python",
                    UNIOP,
                    "-i",
                    fna,
                    "-t",
                    result_dir,
                    "--bin_dir",
                    PRODIGAL_DIR
                ]

                run = subprocess.run(
                    cmd,
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.DEVNULL
                )

                operon_file = os.path.join(
                    result_dir,
                    "uniop.operon"
                )

                if run.returncode == 0 and os.path.exists(operon_file):
                    with open(operon_file) as f:
                        operon_count = max(
                            sum(1 for _ in f) - 1,
                            0
                        )
                else:
                    operon_count = ""

                genome_rows.append([
                    genome_name,
                    record.id,
                    replicon_type,
                    length_bp,
                    operon_count
                ])

                print(
                    f"  {record.id}: {replicon_type}, "
                    f"{length_bp:,} bp, "
                    f"{operon_count} operons"
                )

        # Save this genome immediately
        with open(CSV_FILE, "a", newline="") as f:
            writer = csv.writer(f)
            writer.writerows(genome_rows)

        print("  SAVED")

    except Exception as e:
        print(f"  ERROR: {e}")

print("\nALL AVAILABLE GENOMES FINISHED")
print(f"Results: {CSV_FILE}")
