# ================== PDB Structure Analyzer — Extended (Residue XYZ Export) ==================
import os
import typing as t
import math
import logging
import numpy as np
import pandas as pd
from datetime import datetime
import requests
from Bio.PDB import PDBParser, MMCIFParser
from Bio.PDB.Structure import Structure

try:
    from pymol import cmd, finish_launching
    from pymol.cgo import BEGIN, END, COLOR, VERTEX, LINEWIDTH, LINES
    pymol_available = True
except ImportError:
    pymol_available = False
    print("PyMOL não está instalado ou não foi encontrado. Algumas funcionalidades não estarão disponíveis.")

from scipy.spatial import KDTree

# ---------------- Logging ----------------
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# ---------------- Amino acid names ----------------
AMINOACID_CODES_PT = {
    'ALA': 'Alanina', 'ARG': 'Arginina', 'ASN': 'Asparagina',
    'ASP': 'Ácido Aspártico', 'CYS': 'Cisteína', 'GLU': 'Ácido Glutâmico',
    'GLN': 'Glutamina', 'GLY': 'Glicina', 'HIS': 'Histidina',
    'ILE': 'Isoleucina', 'LEU': 'Leucina', 'LYS': 'Lisina',
    'MET': 'Metionina', 'PHE': 'Fenilalanina', 'PRO': 'Prolina',
    'SER': 'Serina', 'THR': 'Treonina', 'TRP': 'Triptofano',
    'TYR': 'Tirosina', 'VAL': 'Valina',
}
AMINOACID_CODES_EN = {
    'ALA': 'Alanine', 'ARG': 'Arginine', 'ASN': 'Asparagine',
    'ASP': 'Aspartic Acid', 'CYS': 'Cysteine', 'GLU': 'Glutamic Acid',
    'GLN': 'Glutamine', 'GLY': 'Glycine', 'HIS': 'Histidine',
    'ILE': 'Isoleucine', 'LEU': 'Leucine', 'LYS': 'Lysine',
    'MET': 'Methionine', 'PHE': 'Phenylalanine', 'PRO': 'Proline',
    'SER': 'Serine', 'THR': 'Threonine', 'TRP': 'Tryptophan',
    'TYR': 'Tyrosine', 'VAL': 'Valine',
}

# ---------------- Language/UI ----------------
language = 'pt'  # default
AMINOACID_CODES = AMINOACID_CODES_PT  # default

TEXTS = {
    'welcome': {'pt': "Bem-vindo ao Analisador de Estruturas PDB!", 'en': "Welcome to the PDB Structure Analyzer!"},
    'enter_pdb_code': {'pt': "Digite o código PDB ou 'sair' para encerrar: ", 'en': "Enter the PDB code or 'exit' to quit: "},
    'exit_message': {'pt': "Encerrando o programa.", 'en': "Exiting the program."},
    'invalid_pdb': {'pt': "Por favor, tente novamente com um código PDB válido.", 'en': "Please try again with a valid PDB code."},
    'fail_read_structure': {'pt': "Falha ao ler a estrutura do arquivo.", 'en': "Failed to read the structure from the file."},
    'fail_extract_atoms': {'pt': "Falha ao extrair átomos do arquivo.", 'en': "Failed to extract atoms from the file."},
    'menu_options': {
        'pt': "\nMenu Principal:\n"
              "1  - Buscar átomos por critério\n"
              "2  - Buscar átomos por intervalo de coordenadas\n"
              "3  - Calcular as dimensões da estrutura\n"
              "4  - Buscar átomos por distância\n"
              "5  - Buscar átomos por intervalo simultâneo (com exportação por resíduo)\n"
              "6  - Subdividir a estrutura em grids (com exportação por grid)\n"
              "7  - Baixar arquivo\n"
              "8  - Ver histórico de ações\n"
              "9  - Exportar resultados (tabela de átomos)\n"
              "10 - Sobre o Projeto\n"
              "11 - Visualizar estruturas secundárias\n"
              "12 - Visualizar ligantes e íons\n"
              "13 - Colorir estrutura por B-factor\n"
              "14 - Voltar ao início e escolher outro PDB\n"
              "15 - Sair\n"
              "Escolha uma opção e pressione Enter: ",
        'en': "\nMain Menu:\n"
              "1  - Search atoms by criterion\n"
              "2  - Search atoms by coordinate range\n"
              "3  - Calculate structure dimensions\n"
              "4  - Search atoms by distance\n"
              "5  - Simultaneous range search (with per-residue export)\n"
              "6  - Subdivide into grids (with per-grid export)\n"
              "7  - Download file\n"
              "8  - View action history\n"
              "9  - Export results (atom table)\n"
              "10 - About the Project\n"
              "11 - Visualize secondary structures\n"
              "12 - Visualize ligands and ions\n"
              "13 - Color structure by B-factor\n"
              "14 - Return to start and choose another PDB\n"
              "15 - Exit\n"
              "Choose an option and press Enter: "
    },
    'language_choice': {'pt': "Selecione o idioma:", 'en': "Select the language:"},
    'language_options': {'pt': "1 - Português\n2 - Inglês", 'en': "1 - Portuguese\n2 - English"},
    'choice': {'pt': "Escolha: ", 'en': "Choice: "},
    'returning_to_start': {'pt': "Retornando ao início para selecionar outro PDB.", 'en': "Returning to the start to select another PDB."},
    'invalid_option': {'pt': "Opção inválida. Por favor, escolha uma opção válida.", 'en': "Invalid option. Please choose a valid option."},
    'no_header': {'pt': "Cabeçalho não disponível.", 'en': "Header not available."},
    'pdb_loaded': {'pt': "Arquivo PDB carregado com sucesso:", 'en': "PDB file loaded successfully:"},
    'press_enter': {'pt': "Pressione Enter para continuar...", 'en': "Press Enter to continue..."},
    'connection_error': {'pt': "Erro na conexão", 'en': "Connection error"},
    'cannot_download_pdb': {'pt': "Não foi possível baixar o arquivo para o código PDB", 'en': "Unable to download file for PDB code"},
    'unknown_format': {'pt': "Formato desconhecido", 'en': "Unknown format"},
    'error_parsing_file': {'pt': "Erro ao parsear o arquivo", 'en': "Error parsing file"},
    'error_converting_atom': {'pt': "Erro ao converter átomo", 'en': "Error converting atom"},
    'invalid_attribute': {'pt': "Atributo inválido", 'en': "Invalid attribute"},
    'structure_dimensions': {'pt': "Dimensões da Estrutura:", 'en': "Structure Dimensions:"},
    'dimension': {'pt': "Dimensão", 'en': "Dimension"},
    'min': {'pt': "Min", 'en': "Min"},
    'max': {'pt': "Max", 'en': "Max"},
    'no_results_found': {'pt': "Nenhum resultado encontrado.", 'en': "No results found."},
    'name': {'pt': "Nome", 'en': "Name"},
    'residue': {'pt': "Resíduo", 'en': "Residue"},
    'chain': {'pt': "Cadeia", 'en': "Chain"},
    'sequence': {'pt': "Sequência", 'en': "Sequence"},
    'atom_type': {'pt': "Tipo de Átomo", 'en': "Atom Type"},
    'atom': {'pt': "Átomo", 'en': "Atom"},
    'distance': {'pt': "Distância", 'en': "Distance"},
    'error_saving_file': {'pt': "Erro ao salvar o arquivo", 'en': "Error saving file"},
    'file_saved': {'pt': "Arquivo salvo como", 'en': "File saved as"},
    'error_exporting_results': {'pt': "Erro ao exportar resultados", 'en': "Error exporting results"},
    'results_exported': {'pt': "Resultados exportados para", 'en': "Results exported to"},
    'unsupported_format': {'pt': "Formato não suportado.", 'en': "Unsupported format."},
    'analyze_specific_chain': {'pt': "Deseja analisar uma cadeia específica? (S/N): ", 'en': "Do you want to analyze a specific chain? (Y/N): "},
    'enter_chain_identifier': {'pt': "Digite o identificador da cadeia (por exemplo, 'A'): ", 'en': "Enter the chain identifier (e.g., 'A'): "},
    'chain_not_found': {'pt': "A cadeia '{cadeia}' não foi encontrada na estrutura.", 'en': "Chain '{cadeia}' was not found in the structure."},
    'chain_selected': {'pt': "Cadeia '{cadeia}' selecionada para análise.", 'en': "Chain '{cadeia}' selected for analysis."},
    'analyzing_all_chains': {'pt': "Analisando todas as cadeias.", 'en': "Analyzing all chains."},
    'invalid_grid_values': {'pt': "Valores de número de grids inválidos. Por favor, insira números inteiros.", 'en': "Invalid grid number values. Please enter integer numbers."},
    'grid_numbers_must_be_positive': {'pt': "O número de grids em cada eixo deve ser maior que zero.", 'en': "The number of grids on each axis must be greater than zero."},
    'total_subboxes_generated': {'pt': "Total de subcaixas geradas", 'en': "Total subboxes generated"},
    'visualizing_chain_in_pymol': {'pt': "Visualizando apenas a cadeia '{cadeia}' no PyMOL.", 'en': "Visualizing only chain '{cadeia}' in PyMOL."},
    'aminoacids_in_grid': {'pt': "Aminoácidos no grid", 'en': "Amino acids in grid"},
    'instructions_pymol_grid': {'pt': "Para visualizar os aminoácidos em um grid específico no PyMOL, utilize os seguintes comandos:", 'en': "To visualize amino acids in a specific grid in PyMOL, use the following commands:"},
    'show_sticks_example': {'pt': "1. Mostrar os aminoácidos do grid desejado:\n   show sticks, grid_i_j_k_sel\n   Exemplo: show sticks, grid_0_0_0_sel", 'en': "1. Show amino acids of the desired grid:\n   show sticks, grid_i_j_k_sel\n   Example: show sticks, grid_0_0_0_sel"},
    'zoom_example': {'pt': "2. Centralizar a visualização nos aminoácidos selecionados:\n   zoom grid_i_j_k_sel", 'en': "2. Center the view on the selected amino acids:\n   zoom grid_i_j_k_sel"},
    'secondary_structures_highlighted': {'pt': "Estruturas secundárias destacadas no PyMOL.", 'en': "Secondary structures highlighted in PyMOL."},
    'no_ligands_found': {'pt': "Nenhum ligante ou íon encontrado na estrutura.", 'en': "No ligands or ions found in the structure."},
    'ligands_highlighted': {'pt': "Ligantes e íons destacados no PyMOL.", 'en': "Ligands and ions highlighted in PyMOL."},
    'structure_colored_by_bfactor': {'pt': "Estrutura colorida por B-factor no PyMOL.", 'en': "Structure colored by B-factor in PyMOL."},
    'file_saved_locally': {'pt': "Arquivo {filename} salvo localmente.", 'en': "File {filename} saved locally."},
    'action_history': {'pt': "Histórico de Ações", 'en': "Action History"},
    'enter_filename_for_export': {'pt': "Digite o nome do arquivo para exportar os resultados (ex: resultados.csv): ", 'en': "Enter the filename to export the results (e.g., results.csv): "},
    'enter_export_format': {'pt': "Digite o formato de exportação ('csv' ou 'excel'): ", 'en': "Enter the export format ('csv' or 'excel'): "},
    'results_exported_to_file': {'pt': "Resultados exportados para {filename}.", 'en': "Results exported to {filename}."},
    'about_project': {'pt': "Sobre o Projeto", 'en': "About the Project"},
    'project_info': {'pt': "Projeto elaborado pelo professor orientador: Filipe Dalmatti (Lattes: http://lattes.cnpq.br/9691181918031689)\nDiscente: Daniel Castilho (Lattes: http://lattes.cnpq.br/9890731963550827)", 'en': "Project developed by the advisor: Filipe Dalmatti (Lattes: http://lattes.cnpq.br/9691181918031689)\nStudent: Daniel Castilho (Lattes: http://lattes.cnpq.br/9890731963550827)"},
    'pymol_opened': {'pt': "PyMOL aberto com o código PDB", 'en': "PyMOL opened with PDB code"},
    'no_atoms_in_range': {'pt': "Nenhum átomo encontrado no intervalo especificado.", 'en': "No atoms found in the specified range."},
    'enter_search_field': {'pt': "Digite o campo para buscar (id, nome, residuo, cadeia, sequencia, x, y, z, tipo_atomo) ou 'voltar' para o menu principal: ", 'en': "Enter the field to search (id, name, residue, chain, sequence, x, y, z, atom_type) or 'back' to the main menu: "},
    'enter_value_for': {'pt': "Digite o valor para {criterio}: ", 'en': "Enter the value for {criterio}: "},
    'invalid_value_for_criterion': {'pt': "Valor inválido para o tipo de critério escolhido. Tente novamente.", 'en': "Invalid value for the chosen criterion. Please try again."},
    'search_by_criterion_done': {'pt': "Busca por critério '{criterio}' com valor '{valor}' realizada.", 'en': "Search by criterion '{criterio}' with value '{valor}' completed."},
    'enter_axis': {'pt': "Digite o eixo para buscar (x, y, z) ou 'voltar' para o menu principal: ", 'en': "Enter the axis to search (x, y, z) or 'back' to the main menu: "},
    'enter_start_value': {'pt': "Digite o valor inicial do intervalo: ", 'en': "Enter the start value of the range: "},
    'enter_end_value': {'pt': "Digite o valor final do intervalo: ", 'en': "Enter the end value of the range: "},
    'enter_numeric_values': {'pt': "Por favor, insira valores numéricos para os intervalos.", 'en': "Please enter numeric values for the ranges."},
    'search_by_range_done': {'pt': "Busca por intervalo no eixo '{eixo}' de {inicio} a {fim} realizada.", 'en': "Search by range on axis '{eixo}' from {inicio} to {fim} completed."},
    'structure_dimensions_calculated': {'pt': "Cálculo das dimensões da estrutura realizado.", 'en': "Structure dimensions calculated."},
    'enter_max_distance': {'pt': "Digite a distância máxima entre átomos (em Ångström) ou 'voltar' para o menu principal: ", 'en': "Enter the maximum distance between atoms (in Ångström) or 'back' to the main menu: "},
    'enter_numeric_value_for_distance': {'pt': "Por favor, insira um valor numérico para a distância.", 'en': "Please enter a numeric value for the distance."},
    'search_by_max_distance_done': {'pt': "Busca por distância máxima de {distancia} Å realizada.", 'en': "Search by maximum distance of {distancia} Å completed."},
    'enter_start_x': {'pt': "Digite o valor inicial de x: ", 'en': "Enter the start value of x: "},
    'enter_end_x': {'pt': "Digite o valor final de x: ", 'en': "Enter the end value of x: "},
    'enter_start_y': {'pt': "Digite o valor inicial de y: ", 'en': "Enter the start value of y: "},
    'enter_end_y': {'pt': "Digite o valor final de y: ", 'en': "Enter the end value of y: "},
    'enter_start_z': {'pt': "Digite o valor inicial de z: ", 'en': "Enter the start value of z: "},
    'enter_end_z': {'pt': "Digite o valor final de z: ", 'en': "Enter the end value of z: "},
    'include_complete_aminoacids': {'pt': "Deseja incluir aminoácidos completos nos resultados? (S/N): ", 'en': "Do you want to include complete amino acids in the results? (Y/N): "},
    'open_pymol_now': {'pt': "Deseja abrir o PyMOL agora com este resultado? (S/N): ", 'en': "Do you want to open PyMOL now with this result? (Y/N): "},
    'simultaneous_range_search_done': {'pt': "Busca por intervalo simultâneo realizada com inclusão de aminoácidos completos: {incluir}", 'en': "Simultaneous range search completed with inclusion of complete amino acids: {incluir}"},
    'structure_subdivided_into_grids': {'pt': "Subdivisão da estrutura em grids realizada.", 'en': "Structure subdivision into grids completed."},
    'secondary_structures_visualized': {'pt': "Visualização das estruturas secundárias realizada.", 'en': "Secondary structures visualization completed."},
    'ligands_visualized': {'pt': "Visualização de ligantes e íons realizada.", 'en': "Visualization of ligands and ions completed."},
    'structure_colored_by_bfactor_done': {'pt': "Estrutura colorida por B-factor no PyMOL.", 'en': "Structure colored by B-factor in PyMOL."},
    'enter_num_grids_x': {'pt': "Digite o número de grids no eixo X: ", 'en': "Enter the number of grids on the X axis: "},
    'enter_num_grids_y': {'pt': "Digite o número de grids no eixo Y: ", 'en': "Enter the number of grids on the Y axis: "},
    'enter_num_grids_z': {'pt': "Digite o número de grids no eixo Z: ", 'en': "Enter the number of grids on the Z axis: "},
    'pymol_not_installed': {'pt': "PyMOL não está instalado ou não foi encontrado. Algumas funcionalidades não estarão disponíveis.", 'en': "PyMOL is not installed or not found. Some functionalities will not be available."},
    'pymol_visualization_unavailable': {'pt': "Visualização no PyMOL não disponível.", 'en': "Visualization in PyMOL not available."},
    'id_code': {'pt': "Código ID", 'en': "ID Code"},
    'title': {'pt': "Título", 'en': "Title"},
    'classification': {'pt': "Classificação", 'en': "Classification"},
    'deposition_date': {'pt': "Data de Depósito", 'en': "Deposition Date"},
    'resolution': {'pt': "Resolução", 'en': "Resolution"},
    'experimental_method': {'pt': "Método Experimental", 'en': "Experimental Method"},
    'descriptor': {'pt': "Descritor", 'en': "Descriptor"},
    'perform_another_search': {'pt': "Deseja realizar outra busca/utilizar a função novamente? (S/N): ", 'en': "Do you want to perform another search/use the function again? (Y/N): "},

    # Novos textos para exportação XYZ
    'export_xyz_prompt': {
        'pt': "Exportar resíduos encontrados? (S/N): ",
        'en': "Export found residues? (Y/N): "
    },
    'export_xyz_format': {
        'pt': "Formato: 1) XYZ  2) XYZ anotado  3) EXTXYZ  [Enter=3]: ",
        'en': "Format: 1) XYZ  2) Annotated XYZ  3) EXTXYZ  [Enter=3]: "
    },
    'export_xyz_folder': {
        'pt': "Pasta de saída (Enter para 'xyz_exports'): ",
        'en': "Output folder (Enter for 'xyz_exports'): "
    },
    'export_stack_prompt': {
        'pt': "Gerar também um STACK (resíduos enfileirados em +X)? (S/N): ",
        'en': "Also generate a STACK (residues along +X)? (Y/N): "
    },
    'export_stack_spacing': {
        'pt': "Espaçamento entre resíduos no STACK (Å, padrão 5.0): ",
        'en': "Residue spacing for the STACK (Å, default 5.0): "
    },
    'export_xyz_done': {
        'pt': "Exportação concluída. CSV-mestre salvo em:",
        'en': "Export finished. Master CSV saved at:"
    },
    'grid_export_which': {
        'pt': "Exportar de um grid específico? Informe i,j,k (ou Enter para pular): ",
        'en': "Export from a specific grid? Enter i,j,k (or press Enter to skip): "
    },
    'grid_not_found': {
        'pt': "Grid não encontrado ou vazio.",
        'en': "Grid not found or empty."
    },
}

def translate(key): return TEXTS.get(key, {}).get(language, '')

# ---------------- Data structures ----------------
class Atomo:
    def __init__(self, id: int, nome: str, residuo: str, cadeia: str,
                 sequencia: int, x: float, y: float, z: float,
                 b_factor: float, tipo_atomo: str):
        self.id = id
        self.nome = nome
        self.residuo = residuo
        self.cadeia = cadeia
        self.sequencia = sequencia
        self.x = x
        self.y = y
        self.z = z
        self.b_factor = b_factor
        self.tipo_atomo = tipo_atomo

    def distancia_ate(self, outro: 'Atomo') -> float:
        dx, dy, dz = self.x - outro.x, self.y - outro.y, self.z - outro.z
        return math.sqrt(dx * dx + dy * dy + dz * dz)

    def __repr__(self):
        return (f"Atomo(id={self.id}, nome={self.nome}, residuo={self.residuo}, "
                f"cadeia={self.cadeia}, sequencia={self.sequencia}, x={self.x}, "
                f"y={self.y}, z={self.z}, b_factor={self.b_factor}, tipo_atomo={self.tipo_atomo})")

MatrizPDB = t.List[Atomo]

# ---------------- Download / parse ----------------
def download_pdb(pdb_code: str) -> t.Optional[t.Tuple[str, str]]:
    base_url = "https://files.rcsb.org/download/"
    formatos = ['pdb', 'cif']
    for formato in formatos:
        url = f"{base_url}{pdb_code}.{formato}"
        try:
            response = requests.get(url, timeout=10)
            if response.status_code == 200:
                return response.text, formato
        except requests.RequestException as e:
            logger.error(f"{translate('connection_error')}: {e}")
            continue
    logger.error(f"{translate('cannot_download_pdb')} {pdb_code}.")
    return None

def salvar_arquivo(pdb_code: str, content: str, formato: str) -> None:
    filename = f"{pdb_code}.{formato}"
    try:
        with open(filename, 'w') as f: f.write(content)
        print(f"{translate('file_saved')} {filename}")
    except IOError as e:
        logger.error(f"{translate('error_saving_file')}: {e}")

def ler_estrutura_pdb(pdb_content: str, formato: str, pdb_code: str) -> t.Optional[Structure]:
    from io import StringIO
    try:
        if formato == 'pdb':
            parser = PDBParser(QUIET=True)
            estrutura = parser.get_structure(pdb_code, StringIO(pdb_content))
        elif formato == 'cif':
            parser = MMCIFParser(QUIET=True)
            estrutura = parser.get_structure(pdb_code, StringIO(pdb_content))
        else:
            logger.error(f"{translate('unknown_format')}: {formato}")
            return None
        return estrutura
    except Exception as e:
        logger.error(f"{translate('error_parsing_file')} {formato.upper()}: {e}")
        return None

def converter_estrutura_para_atomos(estrutura: Structure) -> MatrizPDB:
    matriz_pdb = []
    for model in estrutura:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    try:
                        # Alguns campos podem estar ausentes dependendo do parser
                        serial = getattr(atom, "serial_number", None)
                        if serial is None:
                            # fallback: gerar um id incremental por ordem de leitura
                            serial = len(matriz_pdb) + 1
                        element = atom.element.strip() if getattr(atom, "element", None) else ''
                        res_id = residue.get_id()[1] if residue.get_id() else None
                        x, y, z = atom.get_coord()
                        b = atom.get_bfactor() if hasattr(atom, "get_bfactor") else 0.0

                        matriz_pdb.append(Atomo(
                            id=serial,
                            nome=atom.get_name(),
                            residuo=residue.get_resname(),
                            cadeia=chain.get_id(),
                            sequencia=int(res_id) if res_id is not None else -1,
                            x=float(x), y=float(y), z=float(z),
                            b_factor=float(b),
                            tipo_atomo=element
                        ))
                    except Exception as e:
                        logger.error(f"{translate('error_converting_atom')}: {e}")
    return matriz_pdb

def ler_cabecalho_de_estrutura(estrutura: Structure) -> t.List[str]:
    header_info = []
    header = estrutura.header
    if 'idcode' in header:
        header_info.append(f"{translate('id_code')}: {header.get('idcode', '')}")
    if 'name' in header:
        header_info.append(f"{translate('title')}: {header.get('name', '')}")
    elif 'title' in header:
        header_info.append(f"{translate('title')}: {header.get('title', '')}")
    if 'classification' in header:
        header_info.append(f"{translate('classification')}: {header.get('classification', '')}")
    if 'deposition_date' in header:
        header_info.append(f"{translate('deposition_date')}: {header.get('deposition_date', '')}")
    if 'resolution' in header:
        resolution = header.get('resolution', '')
        if resolution:
            header_info.append(f"{translate('resolution')}: {resolution} Å")
    return header_info

def buscar_ligantes(estrutura: Structure) -> t.List[str]:
    ligantes = set()
    for model in estrutura:
        for chain in model:
            for residue in chain:
                hetfield = residue.id[0]
                if hetfield.strip() != '':
                    lig = residue.get_resname()
                    if lig != 'HOH':
                        ligantes.add(lig)
    return list(ligantes)

# ---------------- Queries and utilities ----------------
def _normalize_choice(s: str) -> str:
    return (s or '').strip().lower()

def _yes(s: str) -> bool:
    return _normalize_choice(s) in ['s', 'y', 'yes', 'sim']

def _back(s: str) -> bool:
    return _normalize_choice(s) in ['voltar', 'back']

FIELD_MAP_PT = {
    'id': 'id', 'nome': 'nome', 'residuo': 'residuo', 'cadeia': 'cadeia',
    'sequencia': 'sequencia', 'x': 'x', 'y': 'y', 'z': 'z',
    'tipo_atomo': 'tipo_atomo'
}
FIELD_MAP_EN = {
    'id': 'id', 'name': 'nome', 'residue': 'residuo', 'chain': 'cadeia',
    'sequence': 'sequencia', 'x': 'x', 'y': 'y', 'z': 'z',
    'atom_type': 'tipo_atomo'
}

def buscar_por_criterio(matriz_pdb: MatrizPDB, criterio: str, valor: t.Union[int, float, str]) -> MatrizPDB:
    try:
        fieldmap = FIELD_MAP_EN if language == 'en' else FIELD_MAP_PT
        attr = fieldmap.get(criterio.lower(), criterio)
        return [a for a in matriz_pdb if str(getattr(a, attr)).lower() == str(valor).lower()]
    except AttributeError as e:
        logger.error(f"{translate('invalid_attribute')}: {e}")
        return []

def calcular_dimensoes(matriz_pdb: MatrizPDB) -> t.Dict[str, t.Tuple[float, float]]:
    if not matriz_pdb: return {"x": (0.0, 0.0), "y": (0.0, 0.0), "z": (0.0, 0.0)}
    x_vals = [a.x for a in matriz_pdb]; y_vals = [a.y for a in matriz_pdb]; z_vals = [a.z for a in matriz_pdb]
    return {"x": (min(x_vals), max(x_vals)), "y": (min(y_vals), max(y_vals)), "z": (min(z_vals), max(z_vals))}

def imprimir_dimensoes(dimensoes: t.Dict[str, t.Tuple[float, float]]) -> None:
    print(translate('structure_dimensions'))
    for eixo, (mn, mx) in dimensoes.items():
        print(f"{translate('dimension')} {eixo.upper()}: {translate('min')} = {mn}, {translate('max')} = {mx}")

def buscar_por_distancia(matriz_pdb: MatrizPDB, distancia_max: float) -> t.List[t.Tuple[Atomo, Atomo, float]]:
    pontos = np.array([(a.x, a.y, a.z) for a in matriz_pdb])
    kdtree = KDTree(pontos)
    resultados = []
    for i, j in kdtree.query_pairs(distancia_max):
        a1, a2 = matriz_pdb[i], matriz_pdb[j]
        resultados.append((a1, a2, a1.distancia_ate(a2)))
    return resultados

def identificar_interacoes(matriz_pdb: MatrizPDB, distancia_max: float) -> t.List[t.Tuple[Atomo, Atomo, float]]:
    return buscar_por_distancia(matriz_pdb, distancia_max)

def imprimir_resultados(resultados: t.List[t.Union[Atomo, t.Tuple[Atomo, Atomo, float]]]) -> None:
    if not resultados:
        print(translate('no_results_found')); return
    for r in resultados:
        if isinstance(r, Atomo):
            nome_res = AMINOACID_CODES.get(r.residuo, r.residuo)
            print(f"ID: {r.id} | {translate('name')}: {r.nome} | {translate('residue')}: {nome_res} ({r.residuo}) | "
                  f"{translate('chain')}: {r.cadeia} | {translate('sequence')}: {r.sequencia} | "
                  f"X={r.x}, Y={r.y}, Z={r.z} | B-factor: {r.b_factor} | {translate('atom_type')}: {r.tipo_atomo}")
        else:
            a1, a2, d = r
            print(f"{translate('atom')} 1: ID {a1.id} ({a1.nome}) | {translate('atom')} 2: ID {a2.id} ({a2.nome}) | "
                  f"{translate('distance')}: {d:.2f} Å")
        print("")

def buscar_por_intervalo(matriz_pdb: MatrizPDB, eixo: str, inicio: float, fim: float) -> MatrizPDB:
    inicio, fim = sorted([inicio, fim])
    return [a for a in matriz_pdb if inicio <= getattr(a, eixo) <= fim]

def buscar_por_intervalo_simultaneo(matriz_pdb: MatrizPDB,
                                    inicio: t.Tuple[float, float, float],
                                    fim: t.Tuple[float, float, float],
                                    incluir_aminoacidos_completos: bool = False
                                    ) -> t.Tuple[MatrizPDB, t.List[int]]:
    ix, fx = sorted([inicio[0], fim[0]]); iy, fy = sorted([inicio[1], fim[1]]); iz, fz = sorted([inicio[2], fim[2]])
    resultados = [a for a in matriz_pdb if ix <= a.x <= fx and iy <= a.y <= fy and iz <= a.z <= fz]
    if not resultados: print(translate('no_results_found')); return [], []
    if incluir_aminoacidos_completos:
        keys = {(a.residuo, a.cadeia, a.sequencia) for a in resultados}
        resultados = [a for a in matriz_pdb if (a.residuo, a.cadeia, a.sequencia) in keys]
    sequencias_unicas = list({a.sequencia for a in resultados})
    return resultados, sequencias_unicas

def abrir_pymol_comando(atoms: t.List[Atomo], pdb_code: str):
    if not pymol_available:
        print(translate('pymol_visualization_unavailable')); return
    try:
        finish_launching()
        cmd.fetch(pdb_code)
        atom_ids = [a.id for a in atoms]
        selection_str = "id " + "+".join(map(str, atom_ids))
        cmd.select("selecionados", selection_str)
        cmd.show("sticks", "selecionados")
        cmd.zoom("selecionados")
        cmd.deselect()
        print(f"{translate('pymol_opened')} {pdb_code}.")
    except Exception as e:
        print(translate('pymol_visualization_unavailable')); logger.error(f"Error in PyMOL visualization: {e}")

def filtrar_por_cadeia(matriz_pdb: MatrizPDB, cadeia: str) -> MatrizPDB:
    return [a for a in matriz_pdb if a.cadeia == cadeia]

# ---------------- Grid model ----------------
class Grid:
    def __init__(self, i: int, j: int, k: int, nxmin: float, nxmax: float, nymin: float, nymax: float, nzmin: float, nzmax: float):
        self.i, self.j, self.k = i, j, k
        self.nxmin, self.nxmax = nxmin, nxmax
        self.nymin, self.nymax = nymin, nymax
        self.nzmin, self.nzmax = nzmin, nzmax
    def __repr__(self):
        return (f"Grid({self.i}, {self.j}, {self.k}): x({self.nxmin}, {self.nxmax}), "
                f"y({self.nymin}, {self.nymax}), z({self.nzmin}, {self.nzmax})")

def calcular_grids(dimensoes: t.Dict[str, t.Tuple[float, float]], nx: int, ny: int, nz: int) -> t.List[Grid]:
    minx, maxx = dimensoes['x']; miny, maxy = dimensoes['y']; minz, maxz = dimensoes['z']
    dx = (maxx - minx) / nx; dy = (maxy - miny) / ny; dz = (maxz - minz) / nz
    grids = []
    for i in range(nx):
        gx0, gx1 = minx + i*dx, minx + (i+1)*dx
        for j in range(ny):
            gy0, gy1 = miny + j*dy, miny + (j+1)*dy
            for k in range(nz):
                gz0, gz1 = minz + k*dz, minz + (k+1)*dz
                grids.append(Grid(i,j,k,gx0,gx1,gy0,gy1,gz0,gz1))
    return grids

def buscar_aminoacidos_no_grid(matriz_pdb: MatrizPDB, grid: Grid, incluir_aminoacidos_completos: bool=False) -> t.Tuple[MatrizPDB, t.List[int]]:
    resultados = [a for a in matriz_pdb
                  if grid.nxmin <= a.x <= grid.nxmax and grid.nymin <= a.y <= grid.nymax and grid.nzmin <= a.z <= grid.nzmax]
    if not resultados: return [], []
    if incluir_aminoacidos_completos:
        keys = {(a.residuo, a.cadeia, a.sequencia) for a in resultados}
        resultados = [a for a in matriz_pdb if (a.residuo, a.cadeia, a.sequencia) in keys]
    seqs = list({a.sequencia for a in resultados})
    return resultados, seqs

def create_cgo_box(grid: Grid) -> t.List:
    return [
        BEGIN, LINES, COLOR, 1.0, 0.0, 0.0, LINEWIDTH, 2.0,
        VERTEX, grid.nxmin, grid.nymin, grid.nzmin, VERTEX, grid.nxmax, grid.nymin, grid.nzmin,
        VERTEX, grid.nxmax, grid.nymin, grid.nzmin, VERTEX, grid.nxmax, grid.nymax, grid.nzmin,
        VERTEX, grid.nxmax, grid.nymax, grid.nzmin, VERTEX, grid.nxmin, grid.nymax, grid.nzmin,
        VERTEX, grid.nxmin, grid.nymax, grid.nzmin, VERTEX, grid.nxmin, grid.nymin, grid.nzmin,
        VERTEX, grid.nxmin, grid.nymin, grid.nzmax, VERTEX, grid.nxmax, grid.nymin, grid.nzmax,
        VERTEX, grid.nxmax, grid.nymin, grid.nzmax, VERTEX, grid.nxmax, grid.nymax, grid.nzmax,
        VERTEX, grid.nxmax, grid.nymax, grid.nzmax, VERTEX, grid.nxmin, grid.nymax, grid.nzmax,
        VERTEX, grid.nxmin, grid.nymax, grid.nzmax, VERTEX, grid.nxmin, grid.nymin, grid.nzmax,
        VERTEX, grid.nxmin, grid.nymin, grid.nzmin, VERTEX, grid.nxmin, grid.nymin, grid.nzmax,
        VERTEX, grid.nxmax, grid.nymin, grid.nzmin, VERTEX, grid.nxmax, grid.nymin, grid.nzmax,
        VERTEX, grid.nxmax, grid.nymax, grid.nzmin, VERTEX, grid.nxmax, grid.nymax, grid.nzmax,
        VERTEX, grid.nxmin, grid.nymax, grid.nzmin, VERTEX, grid.nxmin, grid.nymax, grid.nzmax, END
    ]

# ---------------- Helpers for element/name inference ----------------
def _infer_element(atom_name: str, tipo_atomo: str) -> str:
    """Inferência robusta do elemento químico."""
    if tipo_atomo and tipo_atomo.strip():
        return tipo_atomo.strip().title()
    s = atom_name.strip()
    if not s:
        return 'X'
    if s[0].isdigit() and len(s) > 1:
        c = s[1]
        return (c + (s[2] if len(s) > 2 and s[2].islower() else '')).title()
    if s[0].isalpha():
        if len(s) > 1 and s[1].islower():
            return (s[0] + s[1]).title()
        return s[0].upper()
    for ch in s:
        if ch.isalpha(): return ch.upper()
    return 'X'

def _group_atoms_by_residue(atoms: t.List[Atomo]) -> t.Dict[t.Tuple[str,str,int], t.List[Atomo]]:
    groups: t.Dict[t.Tuple[str,str,int], t.List[Atomo]] = {}
    for a in atoms:
        key = (a.residuo, a.cadeia, a.sequencia)
        groups.setdefault(key, []).append(a)
    return groups

def _centroid(res_atoms: t.List[Atomo]) -> t.Tuple[float,float,float]:
    xs = [a.x for a in res_atoms]; ys = [a.y for a in res_atoms]; zs = [a.z for a in res_atoms]
    n = max(len(res_atoms),1)
    return (sum(xs)/n, sum(ys)/n, sum(zs)/n)

# ---------------- NEW: XYZ writers (XYZ / Annotated / EXTXYZ) ----------------
def _write_xyz_simple(res_atoms: t.List[Atomo], filepath: str) -> None:
    with open(filepath, 'w') as f:
        f.write(f"{len(res_atoms)}\n")
        f.write(f"Generated by PDB Structure Analyzer — residue export\n")
        for a in res_atoms:
            elem = _infer_element(a.nome, a.tipo_atomo)
            f.write(f"{elem:2s} {a.x: .6f} {a.y: .6f} {a.z: .6f}\n")

def _write_xyz_annotated(res_atoms: t.List[Atomo], filepath: str) -> None:
    """XYZ com colunas extras por linha: element x y z atom_name resname chain resid bfactor atom_id"""
    with open(filepath, 'w') as f:
        f.write(f"{len(res_atoms)}\n")
        f.write(f"grid/residue export — element x y z atom_name resname chain resid bfactor atom_id\n")
        for a in res_atoms:
            elem = _infer_element(a.nome, a.tipo_atomo)
            f.write(f"{elem:2s} {a.x: .6f} {a.y: .6f} {a.z: .6f} "
                    f"{a.nome} {a.residuo} {a.cadeia} {a.sequencia:d} {a.b_factor:.2f} {a.id:d}\n")

def _write_extxyz(res_atoms: t.List[Atomo], filepath: str) -> None:
    """
    EXTXYZ com linha de propriedades — compatível com ASE/OVITO.
    Campos: species, pos(3), atom_name, resname, chain, resid, bfactor, atom_id
    """
    with open(filepath, 'w') as f:
        f.write(f"{len(res_atoms)}\n")
        f.write("Properties=species:S:1:pos:R:3:atom_name:S:1:resname:S:1:chain:S:1:resid:I:1:bfactor:R:1:atom_id:I:1\n")
        for a in res_atoms:
            elem = _infer_element(a.nome, a.tipo_atomo)
            f.write(f"{elem:2s} {a.x: .6f} {a.y: .6f} {a.z: .6f} "
                    f"{a.nome} {a.residuo} {a.cadeia} {a.sequencia:d} {a.b_factor:.2f} {a.id:d}\n")

# ---------------- NEW: Per-residue export orchestrator ----------------
def export_residues_to_xyz(result_atoms: t.List[Atomo], pdb_code: str,
                           outdir: str = "xyz_exports",
                           make_stack: bool = False, spacing: float = 5.0,
                           fmt: str = "extxyz"  # 'xyz' | 'annot' | 'extxyz'
                           ) -> str:
    """
    Gera: (1) um arquivo por resíduo (coords originais) no formato escolhido
          (2) opcional: STACK_<pdb>.xyz (resíduos enfileirados ao longo de +X; formato XYZ simples)
    Retorna o caminho do CSV-mestre com o catálogo da exportação.
    """
    fmt = fmt.lower().strip() or "extxyz"
    os.makedirs(outdir, exist_ok=True)
    groups = _group_atoms_by_residue(result_atoms)

    def _write_residue_file(res_atoms: t.List[Atomo], path: str):
        if fmt == "xyz":
            _write_xyz_simple(res_atoms, path)
        elif fmt in ("annot", "annotated", "xyz_annot"):
            _write_xyz_annotated(res_atoms, path)
        else:
            _write_extxyz(res_atoms, path)

    rows = []
    # 1) arquivos por resíduo
    for (resn, chain, seq), atoms in groups.items():
        atoms_sorted = sorted(atoms, key=lambda a: a.id)
        suffix = {"xyz":"xyz", "annot":"annot.xyz", "annotated":"annot.xyz", "xyz_annot":"annot.xyz"}.get(fmt, "extxyz")
        filename = f"{pdb_code}_{chain}_{resn}{seq}.{suffix}"
        path = os.path.join(outdir, filename)
        _write_residue_file(atoms_sorted, path)
        rows.append({
            'PDB': pdb_code, 'Residue': resn, 'Chain': chain, 'Seq': seq,
            'Atoms': len(atoms_sorted), 'XYZ_File': path, 'Stack_File': ''
        })

    # 2) arquivo de stack (XYZ simples), centralizando cada resíduo no (offset_x, 0, 0)
    stack_path = ""
    if make_stack and groups:
        stack_filename = f"STACK_{pdb_code}.xyz"
        stack_path = os.path.join(outdir, stack_filename)
        ordered = sorted(groups.items(), key=lambda kv: (kv[0][1], kv[0][2], kv[0][0]))
        total_atoms = sum(len(v) for _, v in ordered)
        lines = []
        offset_x = 0.0
        for (resn, chain, seq), atoms in ordered:
            atoms_sorted = sorted(atoms, key=lambda a: a.id)
            cx, cy, cz = _centroid(atoms_sorted)
            for a in atoms_sorted:
                elem = _infer_element(a.nome, a.tipo_atomo)
                x = a.x - cx + offset_x
                y = a.y - cy
                z = a.z - cz
                lines.append(f"{elem:2s} {x: .6f} {y: .6f} {z: .6f}\n")
            offset_x += spacing
        with open(stack_path, 'w') as f:
            f.write(f"{total_atoms}\n")
            f.write(f"Stacked residues along +X; spacing = {spacing:.3f} Å\n")
            for L in lines:
                f.write(L)
        for r in rows:
            r['Stack_File'] = stack_path

    # CSV mestre
    csv_path = os.path.join(outdir, f"{pdb_code}_residue_xyz_index.csv")
    pd.DataFrame(rows).sort_values(by=['Chain','Seq','Residue']).to_csv(csv_path, index=False)
    return csv_path

# ---------------- UI: grids + PyMOL ----------------
def subdividir_estrutura_em_grids(matriz_pdb: MatrizPDB, pdb_code: str) -> None:
    # Seleção de cadeia (opcional)
    while True:
        esc = input(translate('analyze_specific_chain')).strip().upper()
        if esc in ['S','Y','YES']:
            cadeia_sel = input(translate('enter_chain_identifier')).strip().upper()
            matriz_pdb_filtrada = filtrar_por_cadeia(matriz_pdb, cadeia_sel)
            if not matriz_pdb_filtrada:
                print(translate('chain_not_found').format(cadeia=cadeia_sel))
                input(translate('press_enter')); continue
            matriz_pdb_para_grids = matriz_pdb_filtrada
            print(translate('chain_selected').format(cadeia=cadeia_sel)); break
        elif esc in ['N','NO','']:
            print(translate('analyzing_all_chains'))
            cadeia_sel = None
            matriz_pdb_para_grids = matriz_pdb
            break
        else:
            print(translate('invalid_option'))

    dimensoes = calcular_dimensoes(matriz_pdb_para_grids)
    try:
        nx = int(input(translate('enter_num_grids_x')))
        ny = int(input(translate('enter_num_grids_y')))
        nz = int(input(translate('enter_num_grids_z')))
    except ValueError:
        print(translate('invalid_grid_values')); return
    if nx<=0 or ny<=0 or nz<=0:
        print(translate('grid_numbers_must_be_positive')); return

    grids = calcular_grids(dimensoes, nx, ny, nz)
    print(f"\n{translate('total_subboxes_generated')}: {len(grids)}\n")

    if not pymol_available:
        print(translate('pymol_visualization_unavailable'))
    else:
        try:
            finish_launching(); cmd.fetch(pdb_code)
            if cadeia_sel:
                cmd.remove(f"not chain {cadeia_sel}")
                print(translate('visualizing_chain_in_pymol').format(cadeia=cadeia_sel))
            colors = ["red","green","yellow","orange","purple","cyan","magenta","blue"]
            for g in grids:
                grid_name = f"grid_{g.i}_{g.j}_{g.k}"
                box = create_cgo_box(g); cmd.load_cgo(box, grid_name)
                resultados, _ = buscar_aminoacidos_no_grid(matriz_pdb_para_grids, g, incluir_aminoacidos_completos=True)
                if resultados:
                    atom_ids = [a.id for a in resultados]
                    selection_str = "id " + "+".join(map(str, atom_ids))
                    sel_name = f"{grid_name}_sel"; cmd.select(sel_name, selection_str)
                    color = colors[(g.i+g.j+g.k)%len(colors)]; cmd.color(color, sel_name)
                    residuos_unicos = {(a.residuo, a.sequencia, a.cadeia) for a in resultados}
                    print(f"{translate('aminoacids_in_grid')} {grid_name}:")
                    for res, seq, ch in sorted(residuos_unicos, key=lambda x: x[1]):
                        nome_res = AMINOACID_CODES.get(res, res)
                        print(f"- {translate('residue')}: {nome_res} ({res}), {translate('sequence')}: {seq}, {translate('chain')}: {ch}")
                    print("")
            cmd.show("cgo", "grid_*"); cmd.show("cartoon", "all"); cmd.zoom("all")
            print(f"\n{translate('instructions_pymol_grid')}")
            print(translate('show_sticks_example')); print(translate('zoom_example'))
        except Exception as e:
            print(translate('pymol_visualization_unavailable')); logger.error(f"Error in PyMOL visualization: {e}")

    # -------- Export XYZ from a specific grid --------
    grid_str = input(translate('grid_export_which'))
    if grid_str.strip():
        try:
            gi, gj, gk = [int(x.strip()) for x in grid_str.split(',')]
            g_candidates = [g for g in grids if g.i==gi and g.j==gj and g.k==gk]
            if not g_candidates:
                print(translate('grid_not_found')); return
            gsel = g_candidates[0]
            res_atoms, _ = buscar_aminoacidos_no_grid(matriz_pdb_para_grids, gsel, incluir_aminoacidos_completos=True)
            if not res_atoms:
                print(translate('grid_not_found')); return
            if _yes(input(translate('export_xyz_prompt'))):
                fmt_in = input(translate('export_xyz_format')).strip()
                if fmt_in == '1': fmt = 'xyz'
                elif fmt_in == '2': fmt = 'annot'
                else: fmt = 'extxyz'
                outdir = input(translate('export_xyz_folder')).strip() or "xyz_exports"
                do_stack = _yes(input(translate('export_stack_prompt')))
                spacing_in = input(translate('export_stack_spacing')).strip()
                spacing = float(spacing_in) if spacing_in else 5.0
                csv_path = export_residues_to_xyz(res_atoms, pdb_code, outdir, do_stack, spacing, fmt)
                print(f"{translate('export_xyz_done')} {csv_path}")
        except Exception as e:
            print(translate('grid_not_found')); logger.error(e)

# ---------------- Export list of atoms (table) ----------------
def exportar_resultados(resultados: t.List[t.Union[Atomo, t.Tuple[Atomo, Atomo, float]]],
                        filename: str, formato: str = 'csv') -> None:
    try:
        if formato == 'csv':
            if not filename.lower().endswith('.csv'): filename += '.csv'
            df = pd.DataFrame([{
                'ID': a.id, 'Nome': a.nome,
                'Resíduo': AMINOACID_CODES.get(a.residuo, a.residuo),
                'Código Resíduo': a.residuo, 'Cadeia': a.cadeia,
                'Sequência': a.sequencia, 'X': a.x, 'Y': a.y, 'Z': a.z,
                'B-factor': a.b_factor, 'Tipo Atomo': a.tipo_atomo
            } if isinstance(a, Atomo) else {
                'Atomo1_ID': a[0].id, 'Atomo1_Nome': a[0].nome,
                'Atomo2_ID': a[1].id, 'Atomo2_Nome': a[1].nome,
                'Distância': a[2]
            } for a in resultados])
            df.to_csv(filename, index=False)
        elif formato == 'excel':
            if not filename.lower().endswith(('.xlsx', '.xls')): filename += '.xlsx'
            df = pd.DataFrame([{
                'ID': a.id, 'Nome': a.nome,
                'Resíduo': AMINOACID_CODES.get(a.residuo, a.residuo),
                'Código Resíduo': a.residuo, 'Cadeia': a.cadeia,
                'Sequência': a.sequencia, 'X': a.x, 'Y': a.y, 'Z': a.z,
                'B-factor': a.b_factor, 'Tipo Atomo': a.tipo_atomo
            } if isinstance(a, Atomo) else {
                'Atomo1_ID': a[0].id, 'Atomo1_Nome': a[0].nome,
                'Atomo2_ID': a[1].id, 'Atomo2_Nome': a[1].nome,
                'Distância': a[2]
            } for a in resultados])
            df.to_excel(filename, index=False)
        else:
            print(translate('unsupported_format')); return
        print(f"{translate('results_exported')} {filename}")
    except (IOError, ValueError) as e:
        logger.error(f"{translate('error_exporting_results')}: {e}")
        print(f"{translate('error_exporting_results')}: {e}")

# ---------------- Visual helpers ----------------
def visualizar_estruturas_secundarias(pdb_code: str) -> None:
    if not pymol_available:
        print(translate('pymol_visualization_unavailable')); return
    try:
        finish_launching(); cmd.fetch(pdb_code); cmd.show("cartoon", "all"); cmd.color("gray", "all")
        cmd.select("helices", "ss h"); cmd.color("red", "helices")
        cmd.select("sheets", "ss s"); cmd.color("yellow", "sheets")
        cmd.select("turns", "ss l+"); cmd.color("green", "turns")
        cmd.zoom("all"); print(translate('secondary_structures_highlighted'))
    except Exception as e:
        print(translate('pymol_visualization_unavailable')); logger.error(f"Error in PyMOL visualization: {e}")

def visualizar_ligantes(pdb_code: str, ligantes: t.List[str]) -> None:
    if not ligantes:
        print(translate('no_ligands_found')); return
    if not pymol_available:
        print(translate('pymol_visualization_unavailable')); return
    try:
        finish_launching(); cmd.fetch(pdb_code); cmd.show("cartoon", "polymer"); cmd.color("gray", "polymer")
        lig_str = "+".join(ligantes); cmd.select("ligantes", f"resn {lig_str}"); cmd.show("sticks","ligantes"); cmd.color("cyan","ligantes")
        cmd.zoom("ligantes"); print(translate('ligands_highlighted'))
    except Exception as e:
        print(translate('pymol_visualization_unavailable')); logger.error(f"Error in PyMOL visualization: {e}")

def colorir_por_bfactor(pdb_code: str) -> None:
    if not pymol_available:
        print(translate('pymol_visualization_unavailable')); return
    try:
        finish_launching(); cmd.fetch(pdb_code); cmd.show("cartoon","all"); cmd.spectrum("b","blue_white_red","all"); cmd.zoom("all")
        print(translate('structure_colored_by_bfactor'))
    except Exception as e:
        print(translate('pymol_visualization_unavailable')); logger.error(f"Error in PyMOL visualization: {e}")

def limpar_tela() -> None:
    os.system('cls' if os.name == 'nt' else 'clear')

def exibir_header(header: t.List[str]) -> None:
    limpar_tela()
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    if header: print("\n".join(header))
    else: print(translate('no_header'))
    print(f"Timestamp: {timestamp}")

# ---------------- Main menu ----------------
def menu_principal() -> None:
    global language, AMINOACID_CODES
    print(translate('language_choice')); print(translate('language_options'))
    lang_choice = input(translate('choice')).strip()
    if lang_choice == '2':
        language = 'en'; AMINOACID_CODES = AMINOACID_CODES_EN
    else:
        language = 'pt'; AMINOACID_CODES = AMINOACID_CODES_PT

    print(translate('welcome'))
    while True:
        header = []
        pdb_code = input(translate('enter_pdb_code')).strip().upper()
        if (pdb_code == 'SAIR' and language == 'pt') or (pdb_code == 'EXIT' and language == 'en'):
            print(translate('exit_message')); return

        resultado_download = download_pdb(pdb_code)
        if resultado_download is None:
            print(translate('invalid_pdb')); continue

        pdb_content, formato = resultado_download
        estrutura = ler_estrutura_pdb(pdb_content, formato, pdb_code)
        if estrutura is None:
            print(translate('fail_read_structure')); continue

        header = ler_cabecalho_de_estrutura(estrutura)
        matriz_pdb = converter_estrutura_para_atomos(estrutura)
        if not matriz_pdb:
            print(translate('fail_extract_atoms')); continue

        ligantes = buscar_ligantes(estrutura)

        historico_acoes = [f"{translate('pdb_loaded')} {pdb_code}."]
        while True:
            exibir_header(header)
            escolha = input(translate('menu_options')).strip()

            if escolha == "1":
                while True:
                    exibir_header(header)
                    criterio = input(translate('enter_search_field')).strip()
                    if _back(criterio): break
                    valor = input(translate('enter_value_for').format(criterio=criterio))
                    try:
                        fieldmap = FIELD_MAP_EN if language=='en' else FIELD_MAP_PT
                        attr = fieldmap.get(criterio.lower(), criterio)
                        if attr in ['x','y','z']: valor = float(valor)
                        elif attr in ['id','sequencia','sequence']: valor = int(valor)
                    except ValueError:
                        print(translate('invalid_value_for_criterion')); continue
                    resultados = buscar_por_criterio(matriz_pdb, criterio, valor)
                    limpar_tela(); exibir_header(header); imprimir_resultados(resultados)
                    historico_acoes.append(translate('search_by_criterion_done').format(criterio=criterio, valor=valor))
                    input(translate('press_enter'))
                    if not _yes(input(translate('perform_another_search'))): break

            elif escolha == "2":
                while True:
                    exibir_header(header)
                    eixo = input(translate('enter_axis')).strip()
                    if _back(eixo): break
                    try:
                        inicio = float(input(translate('enter_start_value')))
                        fim = float(input(translate('enter_end_value')))
                    except ValueError:
                        print(translate('enter_numeric_values')); continue
                    resultados = buscar_por_intervalo(matriz_pdb, eixo, inicio, fim)
                    limpar_tela(); exibir_header(header); imprimir_resultados(resultados)
                    historico_acoes.append(translate('search_by_range_done').format(eixo=eixo, inicio=inicio, fim=fim))
                    input(translate('press_enter'))
                    if not _yes(input(translate('perform_another_search'))): break

            elif escolha == "3":
                while True:
                    exibir_header(header)
                    dimensoes = calcular_dimensoes(matriz_pdb)
                    imprimir_dimensoes(dimensoes)
                    historico_acoes.append(translate('structure_dimensions_calculated'))
                    input(translate('press_enter'))
                    if not _yes(input(translate('perform_another_search'))): break

            elif escolha == "4":
                while True:
                    exibir_header(header)
                    dist_in = input(translate('enter_max_distance')).strip()
                    if _back(dist_in): break
                    try:
                        dmax = float(dist_in)
                    except ValueError:
                        print(translate('enter_numeric_value_for_distance')); continue
                    resultados = identificar_interacoes(matriz_pdb, dmax)
                    limpar_tela(); exibir_header(header); imprimir_resultados(resultados)
                    historico_acoes.append(translate('search_by_max_distance_done').format(distancia=dmax))
                    input(translate('press_enter'))
                    if not _yes(input(translate('perform_another_search'))): break

            elif escolha == "5":
                while True:
                    exibir_header(header)
                    try:
                        ix = float(input(translate('enter_start_x'))); fx = float(input(translate('enter_end_x')))
                        iy = float(input(translate('enter_start_y'))); fy = float(input(translate('enter_end_y')))
                        iz = float(input(translate('enter_start_z'))); fz = float(input(translate('enter_end_z')))
                    except ValueError:
                        print(translate('enter_numeric_values')); continue
                    incluir = _yes(input(translate('include_complete_aminoacids')))
                    resultados, _ = buscar_por_intervalo_simultaneo(matriz_pdb, (ix,iy,iz), (fx,fy,fz), incluir)
                    limpar_tela(); exibir_header(header); imprimir_resultados(resultados)
                    if resultados and _yes(input(translate('export_xyz_prompt'))):
                        fmt_in = input(translate('export_xyz_format')).strip()
                        if fmt_in == '1': fmt = 'xyz'
                        elif fmt_in == '2': fmt = 'annot'
                        else: fmt = 'extxyz'
                        outdir = input(translate('export_xyz_folder')).strip() or "xyz_exports"
                        do_stack = _yes(input(translate('export_stack_prompt')))
                        spacing_in = input(translate('export_stack_spacing')).strip()
                        spacing = float(spacing_in) if spacing_in else 5.0
                        csv_path = export_residues_to_xyz(resultados, pdb_code, outdir, do_stack, spacing, fmt)
                        print(f"{translate('export_xyz_done')} {csv_path}")
                    historico_acoes.append(translate('simultaneous_range_search_done').format(incluir=incluir))
                    input(translate('press_enter'))
                    if not _yes(input(translate('perform_another_search'))): break

            elif escolha == "6":
                # MODIFICAÇÃO AQUI: Adicionado loop para repetir a função de grid
                while True:
                    exibir_header(header)
                    subdividir_estrutura_em_grids(matriz_pdb, pdb_code)
                    historico_acoes.append(translate('structure_subdivided_into_grids'))
                    input(translate('press_enter'))
                    if not _yes(input(translate('perform_another_search'))):
                        break

            elif escolha == "7":
                exibir_header(header); salvar_arquivo(pdb_code, pdb_content, formato)
                historico_acoes.append(translate('file_saved_locally').format(filename=f"{pdb_code}.{formato}"))
                input(translate('press_enter'))

            elif escolha == "8":
                exibir_header(header); print(f"\n{translate('action_history')}:")
                for acao in historico_acoes: print(f"- {acao}")
                input(translate('press_enter'))

            elif escolha == "9":
                exibir_header(header)
                filename = input(translate('enter_filename_for_export'))
                formato_export = input(translate('enter_export_format')).strip().lower()
                exportar_resultados(matriz_pdb, filename, formato_export)
                historico_acoes.append(translate('results_exported_to_file').format(filename=filename))
                input(translate('press_enter'))

            elif escolha == "10":
                exibir_header(header); print(f"\n{translate('about_project')}:"); print(translate('project_info'))
                input(translate('press_enter'))

            elif escolha == "11":
                visualizar_estruturas_secundarias(pdb_code)
                historico_acoes.append(translate('secondary_structures_visualized'))
                input(translate('press_enter'))

            elif escolha == "12":
                visualizar_ligantes(pdb_code, ligantes)
                historico_acoes.append(translate('ligands_visualized'))
                input(translate('press_enter'))

            elif escolha == "13":
                colorir_por_bfactor(pdb_code)
                historico_acoes.append(translate('structure_colored_by_bfactor_done'))
                input(translate('press_enter'))

            elif escolha == "14":
                print(translate('returning_to_start')); break

            elif escolha == "15":
                print(translate('exit_message')); return

            else:
                print(translate('invalid_option'))

if __name__ == "__main__":
    menu_principal()
# ================== end ==================