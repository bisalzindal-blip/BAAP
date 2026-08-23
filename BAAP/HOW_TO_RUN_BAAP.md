# NOTE: Copy and run the entire cell in Google Colab to run BAAP.
# BAAP — CLEAN INSTALL, CLONE AND RUN
# This cell installs the required Python dependencies, downloads the latest
# BAAP source code from GitHub, configures the BAAP source path, and verifies
# that Biopython/SeqIO is available.
# Run this cell once at the beginning of a new Google Colab session.
# After successful setup, upload your genome FASTA and run BAAP.

# 1. Install Python dependencies
!pip install -q biopython pandas numpy requests tqdm openpyxl

# 2. Remove old BAAP source
!rm -rf /content/BAAP

# 3. Remove old BAAP runtime data
!rm -rf /content/protease_pipeline

# 4. Clone the latest BAAP repository
!git clone -q https://github.com/bisalzindal-blip/BAAP.git /content/BAAP

# 5. Add BAAP source directory
import sys

BAAP_SRC = "/content/BAAP/BAAP/src"

if BAAP_SRC not in sys.path:
    sys.path.insert(0, BAAP_SRC)

# 6. Create BAAP runtime directories
from pathlib import Path

BASE_DIR = Path("/content/protease_pipeline")

for directory in [
    BASE_DIR / "input",
    BASE_DIR / "prokka",
    BASE_DIR / "results",
    BASE_DIR / "logs",
    BASE_DIR / "merops",
    BASE_DIR / "interproscan",
]:
    directory.mkdir(
        parents=True,
        exist_ok=True,
    )

print("BAAP directories initialized successfully.")

# 7. IMPORTANT:
# Remove any previously imported BAAP modules from this Python session.
for module_name in list(sys.modules):
    if module_name == "baap" or module_name.startswith("baap."):
        del sys.modules[module_name]

# 8. Import BAAP AFTER cloning and path setup
from baap import run_pipeline
from baap.config import LOG_FILE, DIRS

# 9. Ensure configured directories exist
LOG_FILE.parent.mkdir(
    parents=True,
    exist_ok=True,
)

for path in DIRS.values():
    Path(path).mkdir(
        parents=True,
        exist_ok=True,
    )

print()
print("=" * 80)
print("BAAP SOURCE")
print("=" * 80)
print(f"BAAP source : {BAAP_SRC}")
print(f"BAAP config : {DIRS}")
print(f"Log file    : {LOG_FILE}")

# 10. Verify the actual merops.py loaded by Python
import baap.merops as merops

print()
print("=" * 80)
print("VERIFYING MEROPS MODULE")
print("=" * 80)
print("Loaded merops.py:")
print(merops.__file__)

print("SeqIO available:", hasattr(merops, "SeqIO"))
print("SeqIO object:", merops.SeqIO)

# 11. Run BAAP
run_pipeline()
