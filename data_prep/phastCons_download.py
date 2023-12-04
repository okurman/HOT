
import os
import sys
from pathlib import Path
import urllib.request

DATA_PATH = Path(os.environ["HOT_DATA"])
PHASTCONS_DIR = DATA_PATH / "phastCons"
PHASTCONS_DIR.mkdir(exist_ok=True)


def check_files_exist(species="vertebrate"):

    chroms = ["chr%d" % _ for _ in range(1, 23)] + ["chrX", "chrY"]

    if species == "vertebrate":
        files = [PHASTCONS_DIR / species /f"{chrom}.phastCons46way.wigFix.gz" for chrom in chroms]
        url_prefix = "https://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/"
    else:
        files = [PHASTCONS_DIR / species / f"{chrom}.phastCons46way.{species}.wigFix.gz" for chrom in chroms]
        url_prefix = f"https://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/{species}/"
    
    if not all([f.exists() for f in files]):
        print(f" {species} phastCons files (or some of them) are not present. Proceeding to download from:")
        print(url_prefix)
        print("\n")

        species_dir = PHASTCONS_DIR/species
        species_dir.mkdir(exist_ok=True, parents=True)

        for file in files:
            if not file.exists():
                file_url = url_prefix + file.name
                print(f"Downloading: {file.name}. Saving to: {file}")
                urllib.request.urlretrieve(file_url, file)


if __name__ == "__main__":

    species = sys.argv[1]
    check_files_exist(species)

