# Q-Core (Quantum Coherence in Redox Enzymes) version alpha 0.001


An interactive Python tool to explore macromolecular structures from the PDB: search and filter atoms, compute bounding-box dimensions, detect short-range contacts, subdivide the structure into 3D grids, visualize in PyMOL, and export per-residue coordinates to XYZ/EXTXYZ (including an optional “stack” layout). The interface is bilingual (Português/English) and runs entirely from the terminal.

Project: Advisor — Filipe Dalmatti (Lattes: http://lattes.cnpq.br/9691181918031689
) · Student — Daniel Castilho (Lattes: http://lattes.cnpq.br/9890731963550827
)

# Highlights

The program downloads a PDB or mmCIF file from RCSB, parses it with Bio.PDB, converts atoms to an internal matrix, and offers an interactive menu. You can search by criteria or coordinate ranges, compute dimensions, find atom–atom pairs within a distance using a KD-tree, isolate a chain, subdivide the box into nx × ny × nz grids and list residues per grid, color by B-factor in PyMOL, highlight secondary structures or ligands, and export results to CSV/Excel. The residue-level exporter writes XYZ, annotated XYZ, or EXTXYZ for each residue and can also build a single STACK file where residues are centered and laid out along +X with configurable spacing.

# Installation

Use Python 3.9+ and install the scientific stack. PyMOL is optional; if not present, visualization features are skipped gracefully.

# create & activate a virtual environment (recommended)
python -m venv .venv
source .venv/bin/activate  # on Windows: .venv\Scripts\activate

# core dependencies
pip install biopython numpy pandas scipy requests

# optional: PyMOL (use your OS package or open-source build if available)
# e.g., on some systems: pip install pymol-open-source

# Running

Simply execute the script and follow the prompts. The first prompt asks for the interface language; the second asks for a PDB code.

python pdb_structure_analyzer.py

When prompted for a PDB code, enter something like 1CRN or 7XYZ. The tool attempts download in PDB then CIF format; if successful, it parses and shows header metadata. The menu then offers numbered actions such as searches, grid subdivision, visualization, and exports. Press Enter when asked to continue; type S/N (pt) or Y/N (en) for yes/no questions.

# Menu Overview

The main operations are presented as numbered options:

Search atoms by criterion — choose a field (id, name/nome, residue/residuo, chain/cadeia, sequence/sequencia, x, y, z, atom_type/tipo_atomo) and a value, then print matching rows.

Search by coordinate range — filter by X, Y or Z between start and end values.

Calculate structure dimensions — prints min/max per axis.

Search atoms by distance — returns pairs of atoms within a cutoff, using scipy.spatial.KDTree.

Simultaneous XYZ box search — filter by rectangular box; optionally include whole amino acids for any atom that falls inside. Immediately after, you can export found residues to XYZ / annotated XYZ / EXTXYZ, choose an output folder, and optionally build a STACK placing each residue centroid along +X with a configurable spacing (default 5.0 Å). A master CSV catalog is saved.

Subdivide into grids — choose nx, ny, nz. The tool shows CGO wireframe boxes in PyMOL (if available), creates PyMOL selections naming residues in each grid, and prints their residue lists in the terminal. You can then export residues from a specific grid to XYZ/EXTXYZ (with the same options as above).

Download file — saves the raw .pdb or .cif locally.

Action history — lists your actions in this session.

Export atom table — writes all current atoms (or pairwise distances if that was the latest result set) to CSV or Excel.

About — credits and project info.

Secondary structures — colors helices/sheets/turns in PyMOL.

Ligands and ions — selects non-water HET groups and shows them as sticks in PyMOL.

Color by B-factor — blue-white-red spectrum in PyMOL.

Restart — return to the PDB prompt.

Exit — quits the program.

If PyMOL is not installed or import fails, visualization options print a message and continue without crashing.

Residue Export Details

Three formats are available:

XYZ — standard two-line header followed by element x y z.

Annotated XYZ — adds per-line annotations: element x y z atom_name resname chain resid bfactor atom_id.

EXTXYZ — writes a Properties= header compatible with ASE/OVITO and includes the same annotations as structured fields.

Each residue is written to its own file inside the chosen output folder (default xyz_exports). A master CSV named <PDB>_residue_xyz_index.csv catalogs residue name, chain, sequence number, atom count, the per-residue file path, and, if requested, the path to the generated STACK file. The STACK is a simple XYZ where every residue is translated so its centroid sits at (offset_x, 0, 0) with offset_x increasing by the chosen spacing.

Element inference follows PDB atom naming rules as robustly as possible; if the element is missing, the first letters of the atom name are used to guess it.

Notes on Chains, Grids, and PyMOL

It is possible to filter a specific chain before grid subdivision. When grids are created, the tool builds CGO wireframes and per-grid PyMOL selections that you can show as sticks and zoom into. Colors are assigned deterministically based on (i + j + k). The terminal prints human-readable residue names in Portuguese or English according to the language setting.

Data Sources & Fair Use

Structure files are fetched from RCSB PDB (https://files.rcsb.org/download/<code>.(pdb|cif)). Please make sure your use complies with RCSB PDB’s terms and data usage policies. This tool is intended for educational and research purposes.

Troubleshooting

If the download fails, check your network and the PDB code. If PyMOL features are unavailable, ensure PyMOL can be imported in Python; otherwise continue using the non-visual features. If OpenGL/GUI issues occur on remote servers, try running PyMOL locally and using only the export and analysis options on the server.

# MIT License

Copyright (c) 2025 Daniel Castilho de Oliveira-Neto and supervisor Filipe Camargo Dalmatti Alves Lima

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

# Português — Visão Geral Rápida

O programa é interativo em terminal e permite baixar e analisar estruturas do PDB, buscar por critérios e intervalos, calcular dimensões, detectar pares de átomos por distância, dividir o espaço em grids 3D, visualizar no PyMOL (opcional) e exportar resíduos para XYZ/EXTXYZ, inclusive em um arquivo STACK com resíduos alinhados ao longo de +X e espaçamento configurável. A exportação gera um CSV-mestre com índice dos arquivos. Se o PyMOL não estiver instalado, as funções de visualização são desativadas sem interromper o programa. Para executar, instale Python 3.9+, biopython, numpy, pandas, scipy, requests e, opcionalmente, PyMOL; então rode python pdb_structure_analyzer.py e siga o menu.
