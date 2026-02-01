````md
# Q-Core (Quantum Coherence in Redox Enzymes) — v0.1.0

An interactive Python terminal tool to explore macromolecular structures from the PDB/mmCIF:
search and filter atoms, compute bounding-box dimensions, detect short-range contacts, subdivide the structure into 3D grids, optionally visualize in PyMOL, and export per-residue coordinates to XYZ / annotated XYZ / EXTXYZ (including an optional **STACK** layout). The interface is bilingual (Português/English).

Project: Advisor — Filipe Dalmatti (Lattes: http://lattes.cnpq.br/9691181918031689) · Student — Daniel Castilho (Lattes: http://lattes.cnpq.br/9890731963550827)

## Highlights

- Download PDB or mmCIF files from RCSB and parse with Bio.PDB.
- Atom searches by criteria and coordinate ranges.
- Structure dimensions (min/max per axis).
- Short-range contact search using `scipy.spatial.KDTree`.
- 3D grid subdivision (**nx × ny × nz**) and residue listing per grid.
- Optional PyMOL visualization (selections, grids, secondary structures, ligands/ions, B-factor coloring).
- Exports to CSV/TSV/Excel/JSON and per-residue XYZ / annotated XYZ / EXTXYZ.
- Optional **STACK** file (residues centered and laid along +X with configurable spacing).

If PyMOL is not installed or import fails, visualization options are skipped gracefully.

## Repository layout

- `qcore.py` — main CLI (terminal) entrypoint.
- `qcore_studio.py` — alternate CLI entrypoint (“Studio”) with optional PyMOL integration and optional Tkinter file dialog.
- `qcore_alpha.py` — experimental/alpha version (may change quickly).
- `Main/` — internal project folder.
- `LICENSE` — MIT License.

## Installation

Requires **Python 3.9+**.

### Create & activate a virtual environment

```bash
python -m venv .venv
source .venv/bin/activate  # Windows: .venv\Scripts\activate
```
````

### Core dependencies

```bash
pip install -U pip
pip install biopython numpy pandas scipy requests
```

### Optional dependencies

- PyMOL (visualization): install according to your environment; features are enabled only if `from pymol import cmd` works.
- Tkinter (file dialog): used only to select local PDB/CIF files; if unavailable, the program asks for the path in the terminal.

## Running

### Q-Core (CLI)

```bash
python qcore.py
```

### Q-Core Studio (CLI)

```bash
python qcore_studio.py
```

Typical flow:

1. Choose interface language (Português/English).
2. Choose input mode: PDB code (download) or local file (`.pdb` / `.cif`).
3. Use the numbered menu options for searches, grids, visualization and exports.

Yes/No prompts accept: **S/N** (pt) or **Y/N** (en).

## Menu overview

Main operations include:

- Search atoms by criterion (id, name/nome, residue/residuo, chain/cadeia, sequence/sequencia, x, y, z, atom_type/tipo_atomo).
- Search by coordinate range (X/Y/Z).
- Calculate structure dimensions.
- Search atom–atom pairs within a distance cutoff (KD-tree).
- Simultaneous XYZ box search (optionally include whole residues; per-residue export).
- Subdivide into grids (nx, ny, nz), list residues per grid, and export by grid.
- Download/save raw PDB/CIF.
- Action history.
- Export results (CSV/TSV/Excel/JSON and summary modes).
- Optional PyMOL tools: secondary structures, ligands/ions, color by B-factor.
- Optional tunneling path (A→B) within a cutoff neighborhood.

## Residue export details

Supported formats:

- **XYZ** — standard header + `element x y z`.
- **Annotated XYZ** — adds per-line metadata (atom name, residue, chain, B-factor, atom id, etc.).
- **EXTXYZ** — includes a `Properties=` header compatible with ASE/OVITO and the same metadata fields.

Outputs:

- One file per residue in the chosen output folder (default `xyz_exports`).
- A master CSV `<PDB>_residue_xyz_index.csv` cataloging residue name, chain, sequence number, atom count, and file paths.
- Optional **STACK** file (`STACK_<PDB>.xyz`) where residues are translated so their centroids lie on +X with configurable spacing (default 5.0 Å).

## Data sources & fair use

Structure files are fetched from RCSB PDB downloads (e.g., `https://files.rcsb.org/download/<code>.(pdb|cif)`). Please ensure your use complies with RCSB PDB’s terms and data usage policies. This tool is intended for educational and research purposes.

## Troubleshooting

- Download failures: check network connectivity and the PDB code.
- PyMOL features unavailable: ensure PyMOL is importable in your Python environment; otherwise continue using non-visual features.
- Remote/OpenGL issues: run PyMOL locally and use export/analysis features on the server.

## License

MIT — see `LICENSE`.

---

# Q-Core (Coerência Quântica em Enzimas Redox) — v0.1.0

Uma ferramenta interativa em Python para **terminal** (CLI) destinada à exploração de estruturas macromoleculares do PDB/mmCIF: busca e filtragem de átomos, cálculo das dimensões da caixa delimitadora, detecção de contatos de curto alcance, subdivisão em grids 3D, visualização opcional no PyMOL e exportação por resíduo para XYZ / XYZ anotado / EXTXYZ (incluindo opção de **STACK**). A interface é bilíngue (Português/English).

Projeto: Orientador — Filipe Dalmatti (Lattes: [http://lattes.cnpq.br/9691181918031689](http://lattes.cnpq.br/9691181918031689)) · Discente — Daniel Castilho (Lattes: [http://lattes.cnpq.br/9890731963550827](http://lattes.cnpq.br/9890731963550827))

## Destaques

- Download de arquivos PDB/mmCIF do RCSB e parse com Bio.PDB.
- Buscas por critérios e por intervalos de coordenadas.
- Dimensões da estrutura (min/max por eixo).
- Busca de contatos por cutoff usando `scipy.spatial.KDTree`.
- Subdivisão em grids 3D (**nx × ny × nz**) e listagem de resíduos por grid.
- Visualização opcional no PyMOL (seleções, grids, estruturas secundárias, ligantes/íons, coloração por B-factor).
- Exportação para CSV/TSV/Excel/JSON e exportação por resíduo para XYZ / XYZ anotado / EXTXYZ.
- Arquivo **STACK** opcional (resíduos centralizados e alinhados em +X com espaçamento configurável).

Se o PyMOL não estiver instalado ou o import falhar, as opções de visualização são ignoradas sem travar.

## Estrutura do repositório

- `qcore.py` — CLI principal (terminal).
- `qcore_studio.py` — entrada alternativa (“Studio”) com integração opcional ao PyMOL e seletor de arquivos opcional via Tkinter.
- `qcore_alpha.py` — versão experimental/alpha (pode mudar rapidamente).
- `Main/` — pasta interna do projeto.
- `LICENSE` — Licença MIT.

## Instalação

Requer **Python 3.9+**.

### Criar e ativar um ambiente virtual

```bash
python -m venv .venv
source .venv/bin/activate  # Windows: .venv\Scripts\activate
```

### Dependências principais

```bash
pip install -U pip
pip install biopython numpy pandas scipy requests
```

### Dependências opcionais

- PyMOL (visualização): instale conforme o seu ambiente; os recursos só são ativados se `from pymol import cmd` funcionar.
- Tkinter (seletor de arquivos): usado apenas para escolher arquivo local; se indisponível, o programa pede o caminho no terminal.

## Execução

### Q-Core (CLI)

```bash
python qcore.py
```

### Q-Core Studio (CLI)

```bash
python qcore_studio.py
```

Fluxo típico:

1. Escolha o idioma (Português/English).
2. Escolha o modo: código PDB (download) ou arquivo local (`.pdb` / `.cif`).
3. Use o menu numerado para buscas, grids, visualização e exportações.

Perguntas Sim/Não aceitam: **S/N** (pt) ou **Y/N** (en).

## Visão geral do menu

Operações principais incluem:

- Busca por critério (id, nome/name, resíduo/residue, cadeia/chain, sequência/sequence, x, y, z, tipo_atomo/atom_type).
- Busca por intervalo de coordenadas (X/Y/Z).
- Cálculo das dimensões da estrutura.
- Busca de pares átomo–átomo por cutoff (KD-tree).
- Busca por caixa XYZ simultânea (opção de incluir resíduos completos; exportação por resíduo).
- Subdivisão em grids (nx, ny, nz), listagem de resíduos por grid e exportação por grid.
- Download/salvar PDB/CIF.
- Histórico de ações.
- Exportação de resultados (CSV/TSV/Excel/JSON e modos de resumo).
- PyMOL opcional: estruturas secundárias, ligantes/íons, coloração por B-factor.
- Caminho de tunelamento (A→B) opcional em vizinhança por cutoff.

## Exportação por resíduo

Formatos:

- **XYZ** — cabeçalho padrão + `elemento x y z`.
- **XYZ anotado** — adiciona metadados por linha (nome do átomo, resíduo, cadeia, B-factor, id do átomo, etc.).
- **EXTXYZ** — inclui cabeçalho `Properties=` compatível com ASE/OVITO e os mesmos metadados.

Saídas:

- Um arquivo por resíduo na pasta escolhida (padrão `xyz_exports`).
- Um CSV-mestre `<PDB>_residue_xyz_index.csv` com resíduos e caminhos.
- **STACK** opcional (`STACK_<PDB>.xyz`) com resíduos alinhados em +X por centróide (espaçamento padrão 5,0 Å).

## Fontes de dados & uso responsável

Arquivos são obtidos via download do RCSB PDB (ex.: `https://files.rcsb.org/download/<code>.(pdb|cif)`). Garanta que seu uso esteja em conformidade com as políticas do RCSB. Ferramenta destinada a fins educacionais e de pesquisa.

## Solução de problemas

- Falha no download: verifique rede e o código PDB.
- PyMOL indisponível: confirme que o PyMOL é importável no ambiente Python; caso contrário use recursos não-visuais.
- Problemas de OpenGL/GUI em servidor: rode o PyMOL localmente e use exportação/análise no servidor.

## Licença

MIT — veja `LICENSE`.
