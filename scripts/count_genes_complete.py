import os
import re
from collections import Counter

folder = "./"

species_list = [
    "Apompejana","Pfijiensis","Anobothrus","Punidentata","Acarldarei",
    "Mpalmata","Skaia","Pectinaria","Phessleri","Pmira","Ppandorae",
    "Hinvalida","Ppalmiformis","Psulfincola","Agunneri","Pgrasslei",
    "Nedwardsi","Acaudata","Pnov"
]

species_set = set(species_list)
total_species = len(species_set)

results = Counter()
files_processed = 0

for filename in os.listdir(folder):
    if not filename.endswith(".fa"):
        continue
    if filename.endswith("_aligned.fa") or filename.endswith("_structure2.fa"):
        continue
    if not filename.startswith("APscaff"):
    	continue

    filepath = os.path.join(folder, filename)
    found_species = set()

    with open(filepath, "r") as f:
        for line in f:
            if line.startswith(">"):
                name = line[1:].strip().split()[0]

                # remove trailing digits (Phessleri3 -> Phessleri)
                base_name = re.sub(r"\d+$", "", name)

                if base_name in species_set:
                    found_species.add(base_name)

    count = len(found_species)
    results[count] += 1
    files_processed += 1

print(f"Total files processed: {files_processed}\n")

for i in range(total_species, -1, -1):
    if results[i] > 0:
        missing = total_species - i
        if missing == 0:
            label = "all species"
        elif missing == 1:
            label = "all but 1 species"
        else:
            label = f"missing {missing} species"
        print(f"{results[i]} files have {label} ({i}/{total_species})")