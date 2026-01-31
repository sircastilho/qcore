"""Q-Core v0.1.0"""
from __future__ import annotations

import os
import math
import logging
import typing as t
import threading
import time
from dataclasses import dataclass
from datetime import datetime

APP_NAME = "Q-Core"
APP_VERSION = "0.1.0"
APP_FULL_NAME = f"{APP_NAME} v{APP_VERSION}"
APP_GUI_NAME = "Q-Core Studio"
__version__ = APP_VERSION
__title__ = APP_FULL_NAME

try:
    import numpy as np
except Exception as e:
    raise ModuleNotFoundError("NumPy não encontrado. Instale com: python -m pip install numpy") from e

try:
    import pandas as pd
except Exception as e:
    raise ModuleNotFoundError("Pandas não encontrado. Instale com: python -m pip install pandas") from e

try:
    import requests
except Exception as e:
    raise ModuleNotFoundError("Requests não encontrado. Instale com: python -m pip install requests") from e

try:
    from scipy.spatial import KDTree
except Exception as e:
    raise ModuleNotFoundError("SciPy não encontrado. Instale com: python -m pip install scipy") from e

try:
    from Bio.PDB.PDBParser import PDBParser
    from Bio.PDB.MMCIFParser import MMCIFParser
    from Bio.PDB.Structure import Structure
except Exception as e:
    raise ModuleNotFoundError("Biopython não encontrado. Instale com: python -m pip install biopython") from e

try:
    from pymol import cmd, finish_launching  # type: ignore[import-not-found]
    from pymol.cgo import BEGIN, END, COLOR, VERTEX, LINEWIDTH, LINES  # type: ignore[import-not-found]
    pymol_available = True
except Exception:
    pymol_available = False
    print("PyMOL não está instalado ou não foi encontrado. Algumas funcionalidades não estarão disponíveis.")

_pymol_click_thread: t.Optional[threading.Thread] = None
_pymol_grid_thread: t.Optional[threading.Thread] = None
_pymol_click_stop: t.Optional[threading.Event] = None
_pymol_last_pick: t.Optional[str] = None
_pymol_last_grid_sel: t.Optional[str] = None
_grid_atom_index: dict[int, set[str]] = {}
_grid_bounds: dict[str, "Grid"] = {}
_grid_phys_stats: dict[str, dict[str, t.Union[float, str]]] = {}
_last_grid_dump: t.Optional[dict[str, t.Any]] = None
_last_pk1: t.Optional[tuple[str, int, str]] = None  # chain, resi, atom name

GRID_DUMP_MAX_INLINE = 120
GRID_DUMP_AUTO_SAVE = False
CLICK_VERBOSE = True

try:
    import tkinter as tk
    from tkinter import filedialog
    tk_available = True
except Exception:
    tk_available = False

logging.basicConfig(
    level=logging.INFO,
    format="%(levelname)s | %(name)s | %(message)s"
)
logger = logging.getLogger("qcore")

AMINOACID_CODES_PT: dict[str, str] = {
    'ALA': 'Alanina', 'ARG': 'Arginina', 'ASN': 'Asparagina',
    'ASP': 'Ácido Aspártico', 'CYS': 'Cisteína', 'GLU': 'Ácido Glutâmico',
    'GLN': 'Glutamina', 'GLY': 'Glicina', 'HIS': 'Histidina',
    'ILE': 'Isoleucina', 'LEU': 'Leucina', 'LYS': 'Lisina',
    'MET': 'Metionina', 'PHE': 'Fenilalanina', 'PRO': 'Prolina',
    'SER': 'Serina', 'THR': 'Treonina', 'TRP': 'Triptofano',
    'TYR': 'Tirosina', 'VAL': 'Valina',
}
AMINOACID_CODES_EN: dict[str, str] = {
    'ALA': 'Alanine', 'ARG': 'Arginine', 'ASN': 'Asparagine',
    'ASP': 'Aspartic Acid', 'CYS': 'Cysteine', 'GLU': 'Glutamic Acid',
    'GLN': 'Glutamine', 'GLY': 'Glycine', 'HIS': 'Histidine',
    'ILE': 'Isoleucine', 'LEU': 'Leucine', 'LYS': 'Lysine',
    'MET': 'Methionine', 'PHE': 'Phenylalanine', 'PRO': 'Proline',
    'SER': 'Serine', 'THR': 'Threonine', 'TRP': 'Tryptophan',
    'TYR': 'Tyrosine', 'VAL': 'Valine',
}

language: str = "pt"
AMINOACID_CODES: dict[str, str] = AMINOACID_CODES_PT

TEXTS: dict[str, dict[str, str]] = {
    'welcome': {'pt': f"Bem-vindo ao {APP_FULL_NAME}!", 'en': f"Welcome to {APP_FULL_NAME}!"},
    'enter_pdb_code': {'pt': "Digite o código PDB ou 'sair' para encerrar: ", 'en': "Enter the PDB code or 'exit' to quit: "},
    'start_mode': {'pt': "Escolha o modo:", 'en': "Choose the mode:"},
    'start_mode_options': {'pt': "1 - Buscar por código PDB\n2 - Carregar arquivo local (.pdb/.cif)",
                           'en': "1 - Search by PDB code\n2 - Load local file (.pdb/.cif)"},
    'enter_local_path': {'pt': "Digite o caminho do arquivo .pdb ou .cif (ou 'voltar'): ",
                         'en': "Enter the path to a .pdb or .cif file (or 'back'): "},
    'file_not_found': {'pt': "Arquivo não encontrado.", 'en': "File not found."},
    'unsupported_file_ext': {'pt': "Extensão não suportada. Use .pdb ou .cif.",
                             'en': "Unsupported extension. Use .pdb or .cif."},
    'file_dialog_unavailable': {'pt': "GUI de seleção de arquivos não disponível. Digite o caminho manualmente.",
                                'en': "File selection GUI not available. Please type the path manually."},
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
              "14 - Trocar arquivo/código (voltar ao início)\n"
              "15 - Fechar PyMOL\n"
              "16 - Exportar último grid clicado\n"
              "17 - Alternar modo de clique (resumo/detalhado)\n"
              "18 - Caminho de tunelamento (A->B)\n"
              "19 - Sair\n"
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
              "14 - Change file/code (return to start)\n"
              "15 - Close PyMOL\n"
              "16 - Export last clicked grid\n"
              "17 - Toggle click mode (summary/detailed)\n"
              "18 - Tunneling path (A->B)\n"
              "19 - Exit\n"
              "Choose an option and press Enter: "
    },
    'language_choice': {'pt': "Selecione o idioma:", 'en': "Select the language:"},
    'language_options': {'pt': "1 - Português\n2 - Inglês", 'en': "1 - Portuguese\n2 - English"},
    'choice': {'pt': "Escolha: ", 'en': "Choice: "},
    'returning_to_start': {'pt': "Retornando ao início para selecionar outro arquivo/código.",
                           'en': "Returning to the start to select another file/code."},
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
    'enter_export_format': {'pt': "Digite o formato de exportação (csv, tsv, excel, json, pdb, cif): ", 'en': "Enter export format (csv, tsv, excel, json, pdb, cif): "},
    'export_mode_prompt': {
        'pt': "Tipo de exportação:\n"
              "1 - Completa (átomos + metadados) — máxima fidelidade\n"
              "2 - Compacta (campos essenciais) — arquivos menores\n"
              "3 - Estatística por resíduo — centroides, bbox, B-factor\n"
              "4 - Estatística por cadeia — visão global + bbox\n"
              "5 - Composição elementar — contagem de elementos\n"
              "6 - Coordenadas XYZ (elemento,x,y,z) — integração simples\n"
              "Escolha: ",
        'en': "Export type:\n"
              "1 - Full (atoms + metadata) — maximum fidelity\n"
              "2 - Compact (essential fields) — smaller files\n"
              "3 - Per-residue stats — centroids, bbox, B-factor\n"
              "4 - Per-chain stats — global view + bbox\n"
              "5 - Elemental composition — element counts\n"
              "6 - XYZ coordinates (element,x,y,z) — simple integration\n"
              "Choice: "
    },
    'export_filter_prompt': {'pt': "Aplicar filtros antes de exportar? (S/N): ", 'en': "Apply filters before export? (Y/N): "},
    'export_filter_chain': {'pt': "Cadeias (ex: A,B) [Enter=ignorar]: ", 'en': "Chains (e.g., A,B) [Enter=skip]: "},
    'export_filter_residue': {'pt': "Resíduos (ex: HIS,GLU) [Enter=ignorar]: ", 'en': "Residues (e.g., HIS,GLU) [Enter=skip]: "},
    'export_filter_seq': {'pt': "Sequência (ex: 10-50) [Enter=ignorar]: ", 'en': "Sequence (e.g., 10-50) [Enter=skip]: "},
    'export_filter_atom': {'pt': "Átomos (ex: CA,CB,OG) [Enter=ignorar]: ", 'en': "Atom names (e.g., CA,CB,OG) [Enter=skip]: "},
    'export_filter_element': {'pt': "Elementos (ex: C,N,O,S) [Enter=ignorar]: ", 'en': "Elements (e.g., C,N,O,S) [Enter=skip]: "},
    'export_filter_bfactor': {'pt': "B-factor (ex: 10-40) [Enter=ignorar]: ", 'en': "B-factor (e.g., 10-40) [Enter=skip]: "},
    'export_filter_coords': {'pt': "Caixa XYZ (ex: x_min,x_max,y_min,y_max,z_min,z_max) [Enter=ignorar]: ", 'en': "XYZ box (e.g., x_min,x_max,y_min,y_max,z_min,z_max) [Enter=skip]: "},
    'export_filter_applied': {'pt': "Filtros aplicados. Átomos restantes: ", 'en': "Filters applied. Remaining atoms: "},
    'results_exported_to_file': {'pt': "Resultados exportados para {filename}.", 'en': "Results exported to {filename}."},
    'about_project': {'pt': "Sobre o Projeto", 'en': "About the Project"},
    'project_info': {
        'pt': f"{APP_FULL_NAME} — Projeto elaborado pelo professor orientador: Filipe Dalmatti (Lattes: http://lattes.cnpq.br/9691181918031689)\nDiscente: Daniel Castilho (Lattes: http://lattes.cnpq.br/9890731963550827)",
        'en': f"{APP_FULL_NAME} — Project developed by advisor Filipe Dalmatti (Lattes: http://lattes.cnpq.br/9691181918031689)\nStudent: Daniel Castilho (Lattes: http://lattes.cnpq.br/9890731963550827)",
    },
    'pymol_opened': {'pt': "PyMOL aberto com o código PDB", 'en': "PyMOL opened with PDB code"},
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
    'distance_ignore_h': {'pt': "Ignorar hidrogênios? (S/N): ", 'en': "Ignore hydrogens? (Y/N): "},
    'distance_exclude_same_residue': {'pt': "Excluir pares do mesmo resíduo? (S/N): ", 'en': "Exclude pairs from same residue? (Y/N): "},
    'enter_start_x': {'pt': "Digite o valor inicial de x: ", 'en': "Enter the start value of x: "},
    'enter_end_x': {'pt': "Digite o valor final de x: ", 'en': "Enter the end value of x: "},
    'enter_start_y': {'pt': "Digite o valor inicial de y: ", 'en': "Enter the start value of y: "},
    'enter_end_y': {'pt': "Digite o valor final de y: ", 'en': "Enter the end value of y: "},
    'enter_start_z': {'pt': "Digite o valor inicial de z: ", 'en': "Enter the start value of z: "},
    'enter_end_z': {'pt': "Digite o valor final de z: ", 'en': "Enter the end value of z: "},
    'include_complete_aminoacids': {'pt': "Deseja incluir aminoácidos completos nos resultados? (S/N): ", 'en': "Do you want to include complete amino acids in the results? (Y/N): "},
    'perform_another_search': {'pt': "Deseja realizar outra busca/utilizar a função novamente? (S/N): ", 'en': "Do you want to perform another search/use the function again? (Y/N): "},
    'simultaneous_range_search_done': {'pt': "Busca por intervalo simultâneo realizada com inclusão de aminoácidos completos: {incluir}", 'en': "Simultaneous range search completed with inclusion of complete amino acids: {incluir}"},
    'structure_subdivided_into_grids': {'pt': "Subdivisão da estrutura em grids realizada.", 'en': "Structure subdivision into grids completed."},
    'secondary_structures_visualized': {'pt': "Visualização das estruturas secundárias realizada.", 'en': "Secondary structures visualization completed."},
    'ligands_visualized': {'pt': "Visualização de ligantes e íons realizada.", 'en': "Visualization of ligands and ions completed."},
    'structure_colored_by_bfactor_done': {'pt': "Estrutura colorida por B-factor no PyMOL.", 'en': "Structure colored by B-factor in PyMOL."},
    'enter_num_grids_x': {'pt': "Digite o número de grids no eixo X: ", 'en': "Enter the number of grids on the X axis: "},
    'enter_num_grids_y': {'pt': "Digite o número de grids no eixo Y: ", 'en': "Enter the number of grids on the Y axis: "},
    'enter_num_grids_z': {'pt': "Digite o número de grids no eixo Z: ", 'en': "Enter the number of grids on the Z axis: "},
    'pymol_visualization_unavailable': {'pt': "Visualização no PyMOL não disponível.", 'en': "Visualization in PyMOL not available."},
    'pymol_closed': {'pt': "PyMOL ocultado/limpo. Você pode continuar usando o programa.", 'en': "PyMOL hidden/cleared. You can keep using the program."},
    'pymol_close_failed': {'pt': "Não foi possível fechar o PyMOL.", 'en': "Unable to close PyMOL."},
    'export_last_grid': {'pt': "Exportar último grid clicado", 'en': "Export last clicked grid"},
    'export_last_grid_none': {'pt': "Nenhum grid foi clicado ainda.", 'en': "No grid has been clicked yet."},
    'export_last_grid_format': {'pt': "Formato (csv/txt) [Enter=csv]: ", 'en': "Format (csv/txt) [Enter=csv]: "},
    'export_last_grid_filename': {'pt': "Nome do arquivo (Enter=padrão): ", 'en': "Filename (Enter=default): "},
    'toggle_click_mode': {'pt': "Alternar modo de clique (resumo/detalhado)", 'en': "Toggle click mode (summary/detailed)"},
    'click_mode_now': {'pt': "Modo de clique atual: {modo}", 'en': "Current click mode: {mode}"},
    'mode_summary': {'pt': "resumo", 'en': "summary"},
    'mode_detailed': {'pt': "detalhado", 'en': "detailed"},
    'tunnel_prompt_start': {'pt': "Ponto A (cadeia,residuo,atomo ex: A,50,CA): ", 'en': "Point A (chain,residue,atom e.g., A,50,CA): "},
    'tunnel_prompt_end': {'pt': "Ponto B (cadeia,residuo,atomo ex: A,120,CA): ", 'en': "Point B (chain,residue,atom e.g., A,120,CA): "},
    'tunnel_not_found': {'pt': "Não foi possível encontrar ambos os pontos.", 'en': "Could not find both points."},
    'tunnel_path_found': {'pt': "Caminho (passos={steps}, peso={weight:.2f}):", 'en': "Path (steps={steps}, weight={weight:.2f}):"},
    'tunnel_no_path': {'pt': "Nenhum caminho encontrado dentro do cutoff.", 'en': "No path found within cutoff."},

    # XYZ export
    'export_xyz_prompt': {'pt': "Exportar resíduos encontrados? (S/N): ", 'en': "Export found residues? (Y/N): "},
    'export_xyz_format': {'pt': "Formato: 1) XYZ  2) XYZ anotado  3) EXTXYZ  [Enter=3]: ",
                          'en': "Format: 1) XYZ  2) Annotated XYZ  3) EXTXYZ  [Enter=3]: "},
    'export_xyz_folder': {'pt': "Pasta de saída (Enter para 'xyz_exports'): ",
                          'en': "Output folder (Enter for 'xyz_exports'): "},
    'export_stack_prompt': {'pt': "Gerar também um STACK (resíduos enfileirados em +X)? (S/N): ",
                            'en': "Also generate a STACK (residues along +X)? (Y/N): "},
    'export_stack_spacing': {'pt': "Espaçamento entre resíduos no STACK (Å, padrão 5.0): ",
                             'en': "Residue spacing for the STACK (Å, default 5.0): "},
    'export_xyz_done': {'pt': "Exportação concluída. CSV-mestre salvo em:",
                        'en': "Export finished. Master CSV saved at:"},

    'grid_export_which': {'pt': "Exportar de um grid específico? Informe i,j,k (ou Enter para pular): ",
                          'en': "Export from a specific grid? Enter i,j,k (or press Enter to skip): "},
    'grid_not_found': {'pt': "Grid não encontrado ou vazio.", 'en': "Grid not found or empty."},
    'grid_inspect_prompt': {'pt': "Inspecionar um grid específico agora? Informe i,j,k (ou Enter para pular): ",
                            'en': "Inspect a specific grid now? Enter i,j,k (or press Enter to skip): "},
    'grid_inspect_depth': {'pt': "Nível de detalhe (1=rápido, 2=detalhado) [Enter=1]: ",
                           'en': "Detail level (1=fast, 2=detailed) [Enter=1]: "},
    'grid_export_stats': {'pt': "Exportar estatísticas de todos os grids? (S/N): ",
                          'en': "Export stats for all grids? (Y/N): "},
    'grid_stats_filename': {'pt': "Nome do arquivo CSV de estatísticas (ex: grid_stats.csv): ",
                            'en': "CSV filename for grid stats (e.g., grid_stats.csv): "},
    'grid_physical_prompt': {'pt': "Aplicar grade fisica (vdW) para maior fidelidade? (S/N): ",
                             'en': "Apply physical (vdW) grid for higher fidelity? (Y/N): "},
    'grid_phys_mode_prompt': {'pt': "Modo fisico (vdW): 1-uniao (sem sobreposicao) 2-soma (com sobreposicao) [Enter=1]: ",
                              'en': "Physical mode (vdW): 1-union (no overlap) 2-sum (with overlap) [Enter=1]: "},
    'grid_expand_bounds_prompt': {'pt': "Expandir limites do grid pelo maior raio vdW? (S/N) [padrao S]: ",
                                  'en': "Expand grid bounds by max vdW radius? (Y/N) [default Y]: "},
    'grid_expand_amount_prompt': {'pt': "Padding adicional em Angstrom (ex: 0.0) [Enter=0.0]: ",
                                  'en': "Extra padding in Angstrom (e.g., 0.0) [Enter=0.0]: "},
    'grid_residue_mode_prompt': {'pt': "Incluir residuos completos no mapeamento do grid? (S/N) [padrao N]: ",
                                 'en': "Include full residues when mapping grids? (Y/N) [default N]: "},
    'grid_voxel_size_prompt': {'pt': "Tamanho do voxel em Angstrom (ex: 0.5) [Enter=0.5]: ",
                               'en': "Voxel size in Angstrom (e.g., 0.5) [Enter=0.5]: "},
    'grid_physical_note': {'pt': "Modo fisico pode ser lento para estruturas grandes. Aguarde...",
                           'en': "Physical mode can be slow for large structures. Please wait..."},
}

def translate(key: str) -> str:
    d = TEXTS.get(key, {})
    return d.get(language, d.get("pt", ""))

class Atomo:
    __slots__ = (
        "id", "nome", "residuo", "cadeia", "sequencia", "icode", "altloc",
        "x", "y", "z", "b_factor", "ocupacao", "tipo_atomo", "hetflag", "modelo"
    )

    def __init__(
        self,
        id: int,
        nome: str,
        residuo: str,
        cadeia: str,
        sequencia: int,
        x: float,
        y: float,
        z: float,
        b_factor: float,
        tipo_atomo: str,
        icode: str = "",
        altloc: str = "",
        ocupacao: float = 1.0,
        hetflag: str = " ",
        modelo: int = 0,
    ):
        self.id = int(id)
        self.nome = str(nome)
        self.residuo = str(residuo)
        self.cadeia = str(cadeia)
        self.sequencia = int(sequencia)
        self.icode = str(icode)
        self.altloc = str(altloc)
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)
        self.b_factor = float(b_factor)
        self.ocupacao = float(ocupacao)
        self.tipo_atomo = str(tipo_atomo)
        self.hetflag = str(hetflag)
        self.modelo = int(modelo)

    def distancia_ate(self, outro: "Atomo") -> float:
        dx, dy, dz = self.x - outro.x, self.y - outro.y, self.z - outro.z
        return math.sqrt(dx * dx + dy * dy + dz * dz)

    def __repr__(self) -> str:
        return (
            f"Atomo(id={self.id}, nome={self.nome}, residuo={self.residuo}, "
            f"cadeia={self.cadeia}, sequencia={self.sequencia}, x={self.x}, "
            f"y={self.y}, z={self.z}, b_factor={self.b_factor}, tipo_atomo={self.tipo_atomo})"
        )

MatrizPDB = t.List[Atomo]

def _infer_element(atom_name: str, tipo_atomo: str) -> str:
    """Inferência robusta do elemento químico."""
    if tipo_atomo and tipo_atomo.strip():
        return tipo_atomo.strip().title()

    s = (atom_name or "").strip()
    if not s:
        return "X"

    # Ex.: "Na", "Cl"
    if len(s) >= 2 and s[0].isalpha() and s[1].islower():
        return (s[0] + s[1]).title()

    # Ex.: "1HG1" -> "H"
    if s[0].isdigit() and len(s) > 1:
        c1 = s[1]
        c2 = s[2] if len(s) > 2 and s[2].islower() else ""
        return (c1 + c2).title()

    if s[0].isalpha():
        return s[0].upper()

    for ch in s:
        if ch.isalpha():
            return ch.upper()

    return "X"

VDW_RADII: dict[str, float] = {
    "H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80, "P": 1.80,
    "F": 1.47, "CL": 1.75, "BR": 1.85, "I": 1.98, "SE": 1.90,
    "NA": 2.27, "K": 2.75, "MG": 1.73, "CA": 2.31,
    "MN": 2.00, "FE": 2.00, "CO": 2.00, "NI": 2.00, "CU": 2.00, "ZN": 2.10,
}

METAL_ELEMENTS: set[str] = {
    "LI", "NA", "K", "RB", "CS", "MG", "CA", "SR", "BA",
    "MN", "FE", "CO", "NI", "CU", "ZN", "CD", "HG",
}

WATER_RESN: set[str] = {"HOH", "WAT", "H2O"}

@dataclass
class TunnelHop:
    i: int
    a_id: int
    b_id: int
    a_spec: str
    b_spec: str
    dist_A: float
    factor: float
    cost: float
    cumulative: float


@dataclass
class TunnelResult:
    ok: bool
    reason: str
    start_spec: str
    end_spec: str
    cutoff_A: float
    direct_dist_A: float
    total_cost: float
    n_nodes: int
    n_hops: int
    nodes: list["Atomo"]
    hops: list[TunnelHop]

def _vdw_radius(element: str) -> float:
    el = (element or "").strip().upper()
    return VDW_RADII.get(el, 1.70 if el else 1.70)

def _max_vdw_radius(atoms: MatrizPDB) -> float:
    if not atoms:
        return 0.0
    return max(_vdw_radius(_infer_element(a.nome, a.tipo_atomo)) for a in atoms)

def _expand_dims(dims: dict[str, tuple[float, float]], pad: float) -> dict[str, tuple[float, float]]:
    if pad <= 0:
        return dict(dims)
    return {
        "x": (dims["x"][0] - pad, dims["x"][1] + pad),
        "y": (dims["y"][0] - pad, dims["y"][1] + pad),
        "z": (dims["z"][0] - pad, dims["z"][1] + pad),
    }

def _is_protein_residue(resn: str) -> bool:
    r = (resn or "").strip().upper()
    return (r in AMINOACID_CODES_EN) or (r in AMINOACID_CODES_PT)

def _group_atoms_by_residue(atoms: t.List[Atomo]) -> dict[tuple[str, str, int, str], t.List[Atomo]]:
    groups: dict[tuple[str, str, int, str], t.List[Atomo]] = {}
    for a in atoms:
        key = (a.residuo, a.cadeia, a.sequencia, a.icode)
        groups.setdefault(key, []).append(a)
    return groups

def _centroid(res_atoms: t.List[Atomo]) -> tuple[float, float, float]:
    n = max(len(res_atoms), 1)
    xs = sum(a.x for a in res_atoms)
    ys = sum(a.y for a in res_atoms)
    zs = sum(a.z for a in res_atoms)
    return (xs / n, ys / n, zs / n)

def _grid_classification_counts(atoms: MatrizPDB) -> dict[str, int]:
    counts = {"protein": 0, "ligand": 0, "water": 0, "ion": 0}
    for a in atoms:
        resn = (a.residuo or "").strip().upper()
        elem = _infer_element(a.nome, a.tipo_atomo).upper()
        if resn in WATER_RESN:
            counts["water"] += 1
        elif elem in METAL_ELEMENTS:
            counts["ion"] += 1
        elif _is_protein_residue(resn):
            counts["protein"] += 1
        else:
            counts["ligand"] += 1
    return counts

def download_pdb(pdb_code: str) -> t.Optional[tuple[str, str]]:
    base_url = "https://files.rcsb.org/download/"
    formatos = ["pdb", "cif"]
    pdb_code = (pdb_code or "").strip().upper()
    if not pdb_code:
        return None

    for formato in formatos:
        url = f"{base_url}{pdb_code}.{formato}"
        try:
            r = requests.get(url, timeout=15)
            if r.status_code == 200 and r.text:
                return r.text, formato
        except requests.RequestException as e:
            logger.error(f"{translate('connection_error')}: {e}")
            continue

    logger.error(f"{translate('cannot_download_pdb')} {pdb_code}.")
    return None

def salvar_arquivo(pdb_code: str, content: str, formato: str) -> None:
    filename = f"{pdb_code}.{formato}"
    try:
        with open(filename, "w", encoding="utf-8") as f:
            f.write(content)
        print(f"{translate('file_saved')} {filename}")
    except IOError as e:
        logger.error(f"{translate('error_saving_file')}: {e}")

def ler_estrutura_pdb(pdb_content: str, formato: str, pdb_code: str) -> t.Optional[Structure]:
    from io import StringIO
    try:
        if formato == "pdb":
            parser = PDBParser(QUIET=True)
            estrutura = parser.get_structure(pdb_code, StringIO(pdb_content))
        elif formato == "cif":
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
    matriz_pdb: MatrizPDB = []
    serial_fallback = 1

    for model in estrutura:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    try:
                        serial = getattr(atom, "serial_number", None)
                        if serial is None:
                            serial = serial_fallback
                            serial_fallback += 1

                        element = ""
                        if getattr(atom, "element", None):
                            element = str(atom.element).strip()

                        if not element:
                            element = _infer_element(atom.get_name(), "")

                        rid = residue.get_id() if residue.get_id() else None
                        res_id = rid[1] if rid else None
                        hetflag = rid[0] if rid else " "
                        icode = rid[2] if rid else ""

                        altloc = ""
                        try:
                            altloc = str(atom.get_altloc() or "").strip()
                        except Exception:
                            altloc = ""

                        x, y, z = atom.get_coord()

                        b = 0.0
                        try:
                            b_raw = atom.get_bfactor()
                            b = float(b_raw) if b_raw is not None else 0.0
                        except Exception:
                            b = 0.0

                        occ = 1.0
                        try:
                            occ_raw = atom.get_occupancy()
                            occ = float(occ_raw) if occ_raw is not None else 1.0
                        except Exception:
                            occ = 1.0

                        try:
                            seq = int(res_id) if res_id is not None else -1
                        except Exception:
                            seq = -1

                        model_id = 0
                        try:
                            model_id = int(getattr(model, "id", 0))
                        except Exception:
                            model_id = 0

                        matriz_pdb.append(
                            Atomo(
                                id=int(serial),
                                nome=atom.get_name(),
                                residuo=residue.get_resname(),
                                cadeia=chain.get_id(),
                                sequencia=seq,
                                x=float(x), y=float(y), z=float(z),
                                b_factor=float(b),
                                tipo_atomo=str(element),
                                icode=str(icode).strip(),
                                altloc=str(altloc).strip(),
                                ocupacao=float(occ),
                                hetflag=str(hetflag),
                                modelo=model_id,
                            )
                        )
                    except Exception as e:
                        logger.error(f"{translate('error_converting_atom')}: {e}")

    return matriz_pdb

def carregar_arquivo_local(path: str) -> t.Optional[tuple[str, str, str]]:
    if not os.path.isfile(path):
        print(translate("file_not_found"))
        return None

    _, ext = os.path.splitext(path)
    ext = ext.lower().lstrip(".")
    if ext not in ("pdb", "cif"):
        print(translate("unsupported_file_ext"))
        return None

    try:
        with open(path, "r", encoding="utf-8", errors="replace") as f:
            content = f.read()
    except IOError as e:
        logger.error(f"{translate('error_parsing_file')} {ext.upper()}: {e}")
        return None

    pdb_code = os.path.splitext(os.path.basename(path))[0].upper()
    return content, ext, pdb_code

def selecionar_arquivo_local_gui() -> t.Optional[str]:
    if not tk_available:
        print(translate("file_dialog_unavailable"))
        return None
    try:
        root = tk.Tk()
        root.withdraw()
        root.attributes("-topmost", True)
        path = filedialog.askopenfilename(
            title="Select PDB/CIF file",
            filetypes=[
                ("PDB/CIF", "*.pdb *.cif"),
                ("PDB", "*.pdb"),
                ("CIF", "*.cif"),
                ("All files", "*.*"),
            ],
        )
        root.destroy()
        return path or None
    except Exception as e:
        logger.error(f"{translate('error_parsing_file')}: {e}")
        return None

def ler_cabecalho_de_estrutura(estrutura: Structure) -> t.List[str]:
    header_info: t.List[str] = []
    header = getattr(estrutura, "header", None) or {}

    if "idcode" in header:
        header_info.append(f"ID Code: {header.get('idcode', '')}")

    if "name" in header:
        header_info.append(f"{translate('title')}: {header.get('name', '')}")
    elif "title" in header:
        header_info.append(f"{translate('title')}: {header.get('title', '')}")

    if "classification" in header:
        header_info.append(f"{translate('classification')}: {header.get('classification', '')}")

    if "deposition_date" in header:
        header_info.append(f"{translate('deposition_date')}: {header.get('deposition_date', '')}")

    if "resolution" in header:
        resolution = header.get("resolution", "")
        if resolution:
            header_info.append(f"{translate('resolution')}: {resolution} Å")

    return header_info

def buscar_ligantes(estrutura: Structure) -> t.List[str]:
    ligantes: set[str] = set()
    for model in estrutura:
        for chain in model:
            for residue in chain:
                hetfield = residue.id[0]
                if str(hetfield).strip() != "":
                    lig = residue.get_resname()
                    if lig != "HOH":
                        ligantes.add(lig)
    return sorted(ligantes)

def _normalize_choice(s: str) -> str:
    return (s or "").strip().lower()

def _yes(s: str) -> bool:
    return _normalize_choice(s) in ["s", "y", "yes", "sim"]

def _yes_default(s: str, default: bool = False) -> bool:
    if s is None or str(s).strip() == "":
        return default
    return _yes(s)

def _back(s: str) -> bool:
    return _normalize_choice(s) in ["voltar", "back"]

def _validate_axis(eixo: str) -> t.Optional[str]:
    eixo = (eixo or "").strip().lower()
    if eixo in ("x", "y", "z"):
        return eixo
    return None

def _safe_input(prompt: str) -> str:
    try:
        return input(prompt)
    except EOFError:
        return ""

def _repeat_prompt() -> bool:
    return _yes(_safe_input(translate("perform_another_search")))

FIELD_MAP_PT: dict[str, str] = {
    "id": "id", "nome": "nome", "residuo": "residuo", "cadeia": "cadeia",
    "sequencia": "sequencia", "x": "x", "y": "y", "z": "z",
    "tipo_atomo": "tipo_atomo"
}
FIELD_MAP_EN: dict[str, str] = {
    "id": "id", "name": "nome", "residue": "residuo", "chain": "cadeia",
    "sequence": "sequencia", "x": "x", "y": "y", "z": "z",
    "atom_type": "tipo_atomo"
}

def buscar_por_criterio(matriz_pdb: MatrizPDB, criterio: str, valor: t.Union[int, float, str]) -> MatrizPDB:
    try:
        fieldmap = FIELD_MAP_EN if language == "en" else FIELD_MAP_PT
        attr = fieldmap.get((criterio or "").lower(), criterio)
        v = str(valor).lower()
        return [a for a in matriz_pdb if str(getattr(a, attr)).lower() == v]
    except AttributeError as e:
        logger.error(f"{translate('invalid_attribute')}: {e}")
        return []

def calcular_dimensoes(matriz_pdb: MatrizPDB) -> dict[str, tuple[float, float]]:
    if not matriz_pdb:
        return {"x": (0.0, 0.0), "y": (0.0, 0.0), "z": (0.0, 0.0)}
    x_vals = [a.x for a in matriz_pdb]
    y_vals = [a.y for a in matriz_pdb]
    z_vals = [a.z for a in matriz_pdb]
    return {"x": (min(x_vals), max(x_vals)), "y": (min(y_vals), max(y_vals)), "z": (min(z_vals), max(z_vals))}

def imprimir_dimensoes(dimensoes: dict[str, tuple[float, float]]) -> None:
    print(translate("structure_dimensions"))
    for eixo, (mn, mx) in dimensoes.items():
        print(f"{translate('dimension')} {eixo.upper()}: {translate('min')} = {mn}, {translate('max')} = {mx}")

def buscar_por_distancia(
    matriz_pdb: MatrizPDB,
    distancia_max: float,
    ignore_h: bool = False,
    exclude_same_residue: bool = False,
) -> t.List[tuple[Atomo, Atomo, float]]:
    if not matriz_pdb or distancia_max <= 0:
        return []

    pontos = np.array([(a.x, a.y, a.z) for a in matriz_pdb], dtype=float)
    kdtree = KDTree(pontos)

    resultados: t.List[tuple[Atomo, Atomo, float]] = []
    for i, j in kdtree.query_pairs(distancia_max):
        a1, a2 = matriz_pdb[i], matriz_pdb[j]

        if ignore_h:
            if _infer_element(a1.nome, a1.tipo_atomo).upper() == "H":
                continue
            if _infer_element(a2.nome, a2.tipo_atomo).upper() == "H":
                continue

        if exclude_same_residue:
            if (a1.cadeia, a1.residuo, a1.sequencia, a1.icode) == (a2.cadeia, a2.residuo, a2.sequencia, a2.icode):
                continue

        resultados.append((a1, a2, a1.distancia_ate(a2)))

    return resultados

def identificar_interacoes(
    matriz_pdb: MatrizPDB,
    distancia_max: float,
    ignore_h: bool = False,
    exclude_same_residue: bool = False,
) -> t.List[tuple[Atomo, Atomo, float]]:
    return buscar_por_distancia(matriz_pdb, distancia_max, ignore_h, exclude_same_residue)

def imprimir_resultados(resultados: t.Sequence[t.Union[Atomo, tuple[Atomo, Atomo, float]]]) -> None:
    if not resultados:
        print(translate("no_results_found"))
        return

    for r in resultados:
        if isinstance(r, Atomo):
            nome_res = AMINOACID_CODES.get(r.residuo, r.residuo)
            print(
                f"ID: {r.id} | {translate('name')}: {r.nome} | {translate('residue')}: {nome_res} ({r.residuo}) | "
                f"{translate('chain')}: {r.cadeia} | {translate('sequence')}: {r.sequencia} | "
                f"X={r.x}, Y={r.y}, Z={r.z} | B-factor: {r.b_factor} | {translate('atom_type')}: {r.tipo_atomo}"
            )
        else:
            a1, a2, d = r
            print(
                f"{translate('atom')} 1: ID {a1.id} ({a1.nome}) | {translate('atom')} 2: ID {a2.id} ({a2.nome}) | "
                f"{translate('distance')}: {d:.2f} Å"
            )
        print("")

def buscar_por_intervalo(matriz_pdb: MatrizPDB, eixo: str, inicio: float, fim: float) -> MatrizPDB:
    eixo_val = _validate_axis(eixo)
    if eixo_val is None:
        logger.error(translate("invalid_attribute"))
        return []
    inicio, fim = sorted([inicio, fim])
    return [a for a in matriz_pdb if inicio <= getattr(a, eixo_val) <= fim]

def buscar_por_intervalo_simultaneo(
    matriz_pdb: MatrizPDB,
    inicio: tuple[float, float, float],
    fim: tuple[float, float, float],
    incluir_aminoacidos_completos: bool = False,
) -> tuple[MatrizPDB, t.List[int]]:
    ix, fx = sorted([inicio[0], fim[0]])
    iy, fy = sorted([inicio[1], fim[1]])
    iz, fz = sorted([inicio[2], fim[2]])

    resultados = [a for a in matriz_pdb if ix <= a.x <= fx and iy <= a.y <= fy and iz <= a.z <= fz]
    if not resultados:
        print(translate("no_results_found"))
        return [], []

    if incluir_aminoacidos_completos:
        keys = {(a.residuo, a.cadeia, a.sequencia, a.icode) for a in resultados}
        resultados = [a for a in matriz_pdb if (a.residuo, a.cadeia, a.sequencia, a.icode) in keys]

    sequencias_unicas = list({a.sequencia for a in resultados})
    return resultados, sequencias_unicas

def abrir_pymol_comando(atoms: t.List[Atomo], pdb_code: str) -> None:
    if not pymol_available:
        print(translate("pymol_visualization_unavailable"))
        return
    if not atoms:
        print(translate("no_results_found"))
        return

    try:
        finish_launching()
        cmd.fetch(pdb_code)
        atom_ids = [a.id for a in atoms]
        selection_str = "id " + "+".join(map(str, atom_ids))
        cmd.select("selecionados", selection_str)
        cmd.show("sticks", "selecionados")
        cmd.zoom("selecionados")
        cmd.deselect()
        _pymol_enable_click_info()
        print(f"{translate('pymol_opened')} {pdb_code}.")
    except Exception as e:
        print(translate("pymol_visualization_unavailable"))
        logger.error(f"Error in PyMOL visualization: {e}")

def _pymol_enable_click_info() -> None:
    """Mostra informações em tempo real ao clicar (seleção pk1)."""
    if not pymol_available:
        return
    try:
        cmd.set("mouse_selection_mode", 0)  # type: ignore[call-arg]
        cmd.label(
            "pk1",
            "\"%s %s%s chain %s alt %s elem %s b=%.2f q=%.2f id=%s (%.3f,%.3f,%.3f)\" % "
            "(name, resn, resi, chain, alt, elem, b, q, ID, x, y, z)"
        )  # type: ignore[call-arg]
        cmd.set("label_color", "yellow")  # type: ignore[call-arg]
        cmd.set("label_size", 14)  # type: ignore[call-arg]
        cmd.set("auto_zoom", 0)  # type: ignore[call-arg]
        cmd.set("pick_tolerance", 1.0)  # type: ignore[call-arg]
    except Exception:
        pass
    _ensure_pymol_click_monitor()

def _pymol_set(name: str, value: t.Union[int, float, str]) -> None:
    """Wrapper para evitar warnings de tipagem e tolerar falhas."""
    try:
        cmd.set(name, value)  # type: ignore[arg-type]
    except Exception:
        pass

def _apply_pymol_style() -> None:
    """Estilo visual premium para competir com PyMOL presets."""
    if not pymol_available:
        return
    try:
        cmd.bg_color("black")  # fundo limpo
        _pymol_set("ray_opaque_background", 0)
        _pymol_set("antialias", 2)
        _pymol_set("ambient", 0.35)
        _pymol_set("specular", 0.5)
        _pymol_set("shininess", 40)
        _pymol_set("cartoon_fancy_helices", 1)
        _pymol_set("cartoon_smooth_loops", 1)
        _pymol_set("cartoon_highlight_color", "white")
        _pymol_set("stick_radius", 0.12)
        _pymol_set("sphere_scale", 0.25)
        _pymol_set("transparency", 0.05)
        _pymol_set("line_width", 2.0)
        _pymol_set("label_shadow_mode", 3)
        _pymol_set("ray_shadow", 0)
        _pymol_set("depth_cue", 1)
        _pymol_set("use_shaders", 0)  # desativa GLSL (evita erros GL/ES)
    except Exception:
        pass

def _format_pick_label(atom: t.Any) -> str:
    model = getattr(atom, "model", "")
    segi = getattr(atom, "segi", "") or ""
    chain = getattr(atom, "chain", "") or ""
    resn = getattr(atom, "resn", "") or ""
    resi = getattr(atom, "resi", "") or ""
    name = getattr(atom, "name", "") or ""
    segi_part = f"/{segi}" if segi else ""
    return f"/{model}{segi_part}/{chain}/{resn}`{resi}/{name}"

def _format_pick_details(atom: t.Any, coords: t.Optional[tuple[float, float, float]] = None) -> str:
    label = _format_pick_label(atom)
    elem = getattr(atom, "elem", "") or ""
    alt = getattr(atom, "alt", "") or ""
    b = getattr(atom, "b", None)
    q = getattr(atom, "q", None)
    atom_id = getattr(atom, "id", None)
    if coords is None:
        try:
            coords = cmd.get_atom_coords("pk1")  # type: ignore[call-arg]
        except Exception:
            coords = None
    if coords is not None:
        x, y, z = coords
        coord_txt = f"xyz=({x:.3f},{y:.3f},{z:.3f})"
    else:
        coord_txt = "xyz=(?, ?, ?)"
    return (
        f"{label} | elem={elem or '?'} | alt={alt or '-'} | "
        f"b={b:.2f} | occ={q:.2f} | id={atom_id} | {coord_txt}"
    )

def _spec_of(a: Atomo) -> str:
    return f"{a.cadeia},{a.sequencia}{a.icode},{a.nome}"

def _is_grid_selection(name: str) -> bool:
    return name.startswith("grid_") and name.endswith("_sel")

def _dump_grid_selection(sel_name: str, max_inline: int = GRID_DUMP_MAX_INLINE, auto_save: bool = GRID_DUMP_AUTO_SAVE) -> None:
    try:
        global _last_grid_dump
        if cmd.count_atoms(sel_name) <= 0:  # type: ignore[call-arg]
            print(f"[PyMOL] Grid {sel_name}: vazio.")
            return

        data: t.List[dict[str, t.Any]] = []
        space: dict[str, t.Any] = {"out": data}
        cmd.iterate_state(  # type: ignore[call-arg]
            1,
            sel_name,
            "out.append({'model': model, 'chain': chain, 'resn': resn, 'resi': resi, "
            "'name': name, 'elem': elem, 'alt': alt, 'b': b, 'q': q, 'id': ID, "
            "'x': x, 'y': y, 'z': z})",
            space=space,
        )

        if not data:
            print(f"[PyMOL] Grid {sel_name}: nenhum átomo.")
            return

        # Organizar: ordenar por cadeia/resíduo/átomo
        data.sort(key=_grid_sort_key)
        _last_grid_dump = {"sel": sel_name, "data": data}

        residues = {}
        elements = {}
        for a in data:
            res_key = (a.get("chain", ""), a.get("resn", ""), a.get("resi", ""))
            residues[res_key] = residues.get(res_key, 0) + 1
            elem = a.get("elem", "?")
            elements[elem] = elements.get(elem, 0) + 1

        print(f"[PyMOL] Grid {sel_name}")
        print(f"- Atomos: {len(data)}")
        print(f"- Residuos: {len(residues)}")
        print(f"- Elementos: {', '.join(f'{k}:{v}' for k, v in sorted(elements.items()))}")
        print("---- Detalhes por atomo ----")

        lines = _format_grid_lines(data)

        if len(lines) <= max_inline:
            for L in lines:
                print(L)
        else:
            for L in lines[:max_inline]:
                print(L)
            print(f"[PyMOL] Mostrando {max_inline}/{len(lines)} átomos. Use 'Exportar último grid clicado' para salvar tudo.")

        if auto_save:
            try:
                csv_name = f"{sel_name}_atoms.csv"
                pd.DataFrame(data).to_csv(csv_name, index=False)
                print(f"[PyMOL] CSV salvo em {csv_name}.")
            except Exception:
                pass
    except Exception as e:
        logger.error(f"Erro ao listar grid {sel_name}: {e}")

def _summarize_grid_selection(sel_name: str, max_residues: int = 25) -> None:
    try:
        if cmd.count_atoms(sel_name) <= 0:  # type: ignore[call-arg]
            print(f"[PyMOL] Grid {sel_name}: vazio.")
            return

        data: t.List[dict[str, t.Any]] = []
        space: dict[str, t.Any] = {"out": data}
        cmd.iterate_state(  # type: ignore[call-arg]
            1,
            sel_name,
            "out.append({'chain': chain, 'resn': resn, 'resi': resi, 'name': name, 'elem': elem})",
            space=space,
        )

        residues: dict[tuple[str, str, str], int] = {}
        elements: dict[str, int] = {}
        for a in data:
            res_key = (a.get("chain", ""), a.get("resn", ""), a.get("resi", ""))
            residues[res_key] = residues.get(res_key, 0) + 1
            elem = a.get("elem", "?")
            elements[elem] = elements.get(elem, 0) + 1

        print(f"[PyMOL] Conteudo do {sel_name}:")
        base = sel_name[:-4] if sel_name.endswith("_sel") else sel_name
        if base in _grid_bounds:
            g = _grid_bounds[base]
            vol = max((g.nxmax - g.nxmin) * (g.nymax - g.nymin) * (g.nzmax - g.nzmin), 0.0)
            print(f"- Bounds: x[{g.nxmin:.3f},{g.nxmax:.3f}] y[{g.nymin:.3f},{g.nymax:.3f}] z[{g.nzmin:.3f},{g.nzmax:.3f}] vol={vol:.3f}")
        if sel_name in _grid_phys_stats:
            p = _grid_phys_stats.get(sel_name, {})
            print(f"- vdW: mode={p.get('vdw_mode','?')} occ={p.get('vdw_occupancy_fraction',0.0):.4f} vol={p.get('vdw_volume_A3',0.0):.2f} A3")
        print(f"- Atomos: {len(data)} | Residuos: {len(residues)} | Elementos: {', '.join(f'{k}:{v}' for k, v in sorted(elements.items()))}")
        print("- Residuos (amostra):")
        res_list = sorted(residues.items(), key=lambda kv: (kv[0][0], kv[0][2], kv[0][1]))
        for (chain, resn, resi), count in res_list[:max_residues]:
            print(f"  {chain:1s} {resn:3s} {resi:>4s}  atoms={count}")
        if len(res_list) > max_residues:
            print(f"  ... (+{len(res_list) - max_residues} residuos)")
    except Exception as e:
        logger.error(f"Erro ao resumir grid {sel_name}: {e}")

def _grid_sort_key(a: dict[str, t.Any]) -> tuple[t.Any, ...]:
    resi_raw = str(a.get("resi", "0")).replace("`", "")
    try:
        resi = int(resi_raw)
    except Exception:
        resi = 0
    return (a.get("chain", ""), resi, a.get("resn", ""), a.get("name", ""))

def _format_grid_lines(data: t.List[dict[str, t.Any]]) -> t.List[str]:
    lines: t.List[str] = []
    for a in data:
        chain = a.get("chain", "")
        resn = a.get("resn", "")
        resi = str(a.get("resi", ""))
        name = a.get("name", "")
        elem = a.get("elem", "?")
        alt = a.get("alt", "-") or "-"
        b = float(a.get("b", 0.0))
        q = float(a.get("q", 0.0))
        atom_id = a.get("id", "?")
        x = float(a.get("x", 0.0))
        y = float(a.get("y", 0.0))
        z = float(a.get("z", 0.0))
        lines.append(
            f"{chain:1s} {resn:3s} {resi:>4s} {name:>4s} "
            f"elem={elem:>2s} alt={alt:>1s} b={b:6.2f} occ={q:4.2f} "
            f"id={atom_id} xyz=({x:8.3f},{y:8.3f},{z:8.3f})"
        )
    return lines

def _pymol_grid_monitor(stop_event: threading.Event) -> None:
    global _pymol_last_grid_sel
    while not stop_event.is_set():
        try:
            try:
                active = cmd.get_names("selections", enabled_only=1)  # type: ignore[call-arg]
            except Exception:
                active = cmd.get_names("selections")  # type: ignore[call-arg]
            active_grids = [n for n in active if _is_grid_selection(n)]
            if active_grids:
                sel = sorted(active_grids)[0]
                if sel != _pymol_last_grid_sel:
                    _pymol_last_grid_sel = sel
                    limpar_tela()
                    _dump_grid_selection(sel, max_inline=GRID_DUMP_MAX_INLINE, auto_save=GRID_DUMP_AUTO_SAVE)
            else:
                if _pymol_last_grid_sel is not None:
                    _pymol_last_grid_sel = None
                    limpar_tela()
        except Exception:
            pass
        time.sleep(0.3)

def _pymol_click_monitor(stop_event: threading.Event) -> None:
    global _pymol_last_pick
    while not stop_event.is_set():
        try:
            try:
                selections = cmd.get_names("selections")  # type: ignore[call-arg]
            except Exception:
                selections = []
            if "pk1" not in selections:
                _pymol_last_pick = None
                time.sleep(0.2)
                continue

            if cmd.count_atoms("pk1") <= 0:  # type: ignore[call-arg]
                _pymol_last_pick = None
                time.sleep(0.2)
                continue

            m = cmd.get_model("pk1")  # type: ignore[call-arg]
            if m and getattr(m, "atom", None):
                a = m.atom[0]
                try:
                    coords = cmd.get_atom_coords("pk1")  # type: ignore[call-arg]
                except Exception:
                    coords = None
                pick = _format_pick_label(a)
                if pick and pick != _pymol_last_pick:
                    _pymol_last_pick = pick
                    limpar_tela()
                    print(f"[PyMOL] Click: {_format_pick_details(a, coords)}")
                    try:
                        atom_id = int(getattr(a, "id", -1))
                    except Exception:
                        atom_id = -1
                    if atom_id in _grid_atom_index:
                        grids = sorted(_grid_atom_index.get(atom_id, set()))
                        print(f"[PyMOL] Grid(s) do clique: {', '.join(grids)}")
                        if grids:
                            _summarize_grid_selection(grids[0])
                            if CLICK_VERBOSE:
                                _dump_grid_selection(grids[0], max_inline=50, auto_save=False)
                    else:
                        print("[PyMOL] Clique fora de qualquer grid mapeado.")
            else:
                _pymol_last_pick = None
        except Exception:
            pass
        time.sleep(0.2)

def _ensure_pymol_click_monitor() -> None:
    global _pymol_click_thread, _pymol_grid_thread, _pymol_click_stop
    if not (_pymol_click_thread and _pymol_click_thread.is_alive()):
        _pymol_click_stop = threading.Event()
        _pymol_click_thread = threading.Thread(
            target=_pymol_click_monitor,
            args=(_pymol_click_stop,),
            daemon=True,
            name="pymol_click_monitor",
        )
        _pymol_click_thread.start()
    if not (_pymol_grid_thread and _pymol_grid_thread.is_alive()):
        if _pymol_click_stop is None:
            _pymol_click_stop = threading.Event()
        _pymol_grid_thread = threading.Thread(
            target=_pymol_grid_monitor,
            args=(_pymol_click_stop,),
            daemon=True,
            name="pymol_grid_monitor",
        )
        _pymol_grid_thread.start()

def filtrar_por_cadeia(matriz_pdb: MatrizPDB, cadeia: str) -> MatrizPDB:
    c = (cadeia or "").strip()
    return [a for a in matriz_pdb if a.cadeia == c]

class Grid:
    def __init__(
        self,
        i: int,
        j: int,
        k: int,
        nxmin: float,
        nxmax: float,
        nymin: float,
        nymax: float,
        nzmin: float,
        nzmax: float,
        is_last_x: bool,
        is_last_y: bool,
        is_last_z: bool,
    ):
        self.i, self.j, self.k = int(i), int(j), int(k)
        self.nxmin, self.nxmax = float(nxmin), float(nxmax)
        self.nymin, self.nymax = float(nymin), float(nymax)
        self.nzmin, self.nzmax = float(nzmin), float(nzmax)
        self.is_last_x = bool(is_last_x)
        self.is_last_y = bool(is_last_y)
        self.is_last_z = bool(is_last_z)

    def __repr__(self) -> str:
        return (
            f"Grid({self.i}, {self.j}, {self.k}): x({self.nxmin}, {self.nxmax}), "
            f"y({self.nymin}, {self.nymax}), z({self.nzmin}, {self.nzmax})"
        )

def calcular_grids(dimensoes: dict[str, tuple[float, float]], nx: int, ny: int, nz: int) -> t.List[Grid]:
    minx, maxx = dimensoes["x"]
    miny, maxy = dimensoes["y"]
    minz, maxz = dimensoes["z"]

    xs = np.linspace(minx, maxx, nx + 1, dtype=float)
    ys = np.linspace(miny, maxy, ny + 1, dtype=float)
    zs = np.linspace(minz, maxz, nz + 1, dtype=float)

    grids: t.List[Grid] = []
    for i in range(nx):
        gx0, gx1 = float(xs[i]), float(xs[i + 1])
        for j in range(ny):
            gy0, gy1 = float(ys[j]), float(ys[j + 1])
            for k in range(nz):
                gz0, gz1 = float(zs[k]), float(zs[k + 1])
                grids.append(
                    Grid(
                        i, j, k,
                        gx0, gx1, gy0, gy1, gz0, gz1,
                        is_last_x=(i == nx - 1),
                        is_last_y=(j == ny - 1),
                        is_last_z=(k == nz - 1),
                    )
                )
    return grids

def buscar_aminoacidos_no_grid(matriz_pdb: MatrizPDB, grid: Grid, incluir_aminoacidos_completos: bool = False) -> tuple[MatrizPDB, t.List[int]]:
    eps = 1e-9

    def _in_range(v: float, vmin: float, vmax: float, inclusive_max: bool) -> bool:
        if inclusive_max:
            return (vmin - eps) <= v <= (vmax + eps)
        return (vmin - eps) <= v < (vmax - eps)

    resultados = [
        a for a in matriz_pdb
        if _in_range(a.x, grid.nxmin, grid.nxmax, grid.is_last_x)
        and _in_range(a.y, grid.nymin, grid.nymax, grid.is_last_y)
        and _in_range(a.z, grid.nzmin, grid.nzmax, grid.is_last_z)
    ]

    if not resultados:
        return [], []

    if incluir_aminoacidos_completos:
        keys = {(a.residuo, a.cadeia, a.sequencia, a.icode) for a in resultados}
        resultados = [a for a in matriz_pdb if (a.residuo, a.cadeia, a.sequencia, a.icode) in keys]

    seqs = list({a.sequencia for a in resultados})
    return resultados, seqs

def _grid_summary(atoms: MatrizPDB, detailed: bool = False, phys: t.Optional[dict[str, t.Union[str, float]]] = None) -> str:
    if not atoms:
        return translate("grid_not_found")

    residues = {(a.residuo, a.sequencia, a.cadeia, a.icode) for a in atoms}
    chains = {a.cadeia for a in atoms}
    counts = _grid_classification_counts(atoms)

    b_vals = [a.b_factor for a in atoms]
    b_mean = float(np.mean(b_vals)) if b_vals else 0.0
    b_std = float(np.std(b_vals)) if b_vals else 0.0
    altlocs = len({a.altloc for a in atoms if a.altloc})
    models = len({a.modelo for a in atoms})

    elems: dict[str, int] = {}
    for a in atoms:
        e = _infer_element(a.nome, a.tipo_atomo)
        elems[e] = elems.get(e, 0) + 1

    parts = [
        f"Atoms: {len(atoms)}",
        f"Residues: {len(residues)}",
        f"Chains: {len(chains)}",
        f"Protein: {counts['protein']}",
        f"Ligand: {counts['ligand']}",
        f"Water: {counts['water']}",
        f"Ion: {counts['ion']}",
        f"B-factor mean: {b_mean:.2f}",
        f"B-factor std: {b_std:.2f}",
        f"AltLocs: {altlocs}",
        f"Models: {models}",
    ]

    if phys:
        parts.append(f"vdW volume: {phys.get('vdw_volume_A3', 0.0):.2f} A3")
        parts.append(f"vdW occupancy: {phys.get('vdw_occupancy_fraction', 0.0):.4f}")
        if "vdw_mode" in phys:
            parts.append(f"vdW mode: {phys.get('vdw_mode')}")

    if detailed:
        top_elems = ", ".join(f"{k}:{v}" for k, v in sorted(elems.items()))
        res_list = sorted(residues, key=lambda x: (x[2], x[1], x[0], x[3]))
        sample = "; ".join([f"{r[0]}{r[1]}{r[3]}:{r[2]}" for r in res_list[:20]])
        parts.append(f"Elements: {top_elems}")
        parts.append(f"Residues(sample): {sample}{' ...' if len(res_list) > 20 else ''}")

    return " | ".join(parts)

def create_cgo_box(grid: Grid) -> t.List[t.Any]:
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
        VERTEX, grid.nxmin, grid.nymax, grid.nzmin, VERTEX, grid.nxmin, grid.nymax, grid.nzmax,
        END
    ]

def _compute_physical_grid_stats(
    atoms: MatrizPDB,
    dims: dict[str, tuple[float, float]],
    nx: int,
    ny: int,
    nz: int,
    voxel: float,
    mode: str = "union",
) -> dict[tuple[int, int, int], dict[str, t.Union[str, float]]]:
    minx, maxx = dims["x"]
    miny, maxy = dims["y"]
    minz, maxz = dims["z"]

    dx = (maxx - minx) / nx if nx > 0 else 1.0
    dy = (maxy - miny) / ny if ny > 0 else 1.0
    dz = (maxz - minz) / nz if nz > 0 else 1.0

    nxv = max(int(math.ceil((maxx - minx) / voxel)), 1)
    nyv = max(int(math.ceil((maxy - miny) / voxel)), 1)
    nzv = max(int(math.ceil((maxz - minz) / voxel)), 1)
    voxel_vol = voxel ** 3

    grid_counts: dict[tuple[int, int, int], int] = {}
    seen: set[int] = set()
    mode = (mode or "union").strip().lower()

    for a in atoms:
        elem = _infer_element(a.nome, a.tipo_atomo)
        r = _vdw_radius(elem)
        r2 = r * r

        ix0 = max(int(math.floor((a.x - r - minx) / voxel)), 0)
        iy0 = max(int(math.floor((a.y - r - miny) / voxel)), 0)
        iz0 = max(int(math.floor((a.z - r - minz) / voxel)), 0)
        ix1 = min(int(math.floor((a.x + r - minx) / voxel)), nxv - 1)
        iy1 = min(int(math.floor((a.y + r - miny) / voxel)), nyv - 1)
        iz1 = min(int(math.floor((a.z + r - minz) / voxel)), nzv - 1)

        for ix in range(ix0, ix1 + 1):
            cx = minx + (ix + 0.5) * voxel
            dx2 = (cx - a.x) ** 2
            if dx2 > r2:
                continue
            for iy in range(iy0, iy1 + 1):
                cy = miny + (iy + 0.5) * voxel
                dy2 = (cy - a.y) ** 2
                if dx2 + dy2 > r2:
                    continue
                for iz in range(iz0, iz1 + 1):
                    cz = minz + (iz + 0.5) * voxel
                    if dx2 + dy2 + (cz - a.z) ** 2 > r2:
                        continue

                    lin = ix + nxv * (iy + nyv * iz)
                    if mode == "union":
                        if lin in seen:
                            continue
                        seen.add(lin)

                    gi = min(max(int((cx - minx) / dx), 0), nx - 1)
                    gj = min(max(int((cy - miny) / dy), 0), ny - 1)
                    gk = min(max(int((cz - minz) / dz), 0), nz - 1)
                    key = (gi, gj, gk)
                    grid_counts[key] = grid_counts.get(key, 0) + 1

    stats: dict[tuple[int, int, int], dict[str, t.Union[str, float]]] = {}
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                count = grid_counts.get((i, j, k), 0)
                grid_vol = dx * dy * dz
                occ_vol = count * voxel_vol
                stats[(i, j, k)] = {
                    "vdw_voxels": float(count),
                    "vdw_volume_A3": float(occ_vol),
                    "vdw_occupancy_fraction": float(occ_vol / grid_vol) if grid_vol > 0 else 0.0,
                    "vdw_mode": mode,
                }
    return stats

def _write_xyz_simple(res_atoms: t.List[Atomo], filepath: str) -> None:
    with open(filepath, "w", encoding="utf-8") as f:
        f.write(f"{len(res_atoms)}\n")
        f.write(f"Generated by {APP_FULL_NAME} — residue export\n")
        for a in res_atoms:
            elem = _infer_element(a.nome, a.tipo_atomo)
            f.write(f"{elem:2s} {a.x: .6f} {a.y: .6f} {a.z: .6f}\n")

def _write_xyz_annotated(res_atoms: t.List[Atomo], filepath: str) -> None:
    """XYZ com colunas extras por linha."""
    with open(filepath, "w", encoding="utf-8") as f:
        f.write(f"{len(res_atoms)}\n")
        f.write("grid/residue export — element x y z atom_name resname chain resid bfactor atom_id\n")
        for a in res_atoms:
            elem = _infer_element(a.nome, a.tipo_atomo)
            f.write(
                f"{elem:2s} {a.x: .6f} {a.y: .6f} {a.z: .6f} "
                f"{a.nome} {a.residuo} {a.cadeia} {a.sequencia:d} {a.b_factor:.2f} {a.id:d}\n"
            )

def _write_extxyz(res_atoms: t.List[Atomo], filepath: str) -> None:
    """
    EXTXYZ com linha de propriedades — compatível com ASE/OVITO.
    Campos: species, pos(3), atom_name, resname, chain, resid, bfactor, atom_id
    """
    with open(filepath, "w", encoding="utf-8") as f:
        f.write(f"{len(res_atoms)}\n")
        f.write("Properties=species:S:1:pos:R:3:atom_name:S:1:resname:S:1:chain:S:1:resid:I:1:bfactor:R:1:atom_id:I:1\n")
        for a in res_atoms:
            elem = _infer_element(a.nome, a.tipo_atomo)
            f.write(
                f"{elem:2s} {a.x: .6f} {a.y: .6f} {a.z: .6f} "
                f"{a.nome} {a.residuo} {a.cadeia} {a.sequencia:d} {a.b_factor:.2f} {a.id:d}\n"
            )

def export_residues_to_xyz(
    result_atoms: t.List[Atomo],
    pdb_code: str,
    outdir: str = "xyz_exports",
    make_stack: bool = False,
    spacing: float = 5.0,
    fmt: str = "extxyz",  # 'xyz' | 'annot' | 'extxyz'
) -> str:
    """
    Gera:
      (1) um arquivo por resíduo (coords originais) no formato escolhido
      (2) opcional: STACK_<pdb>.xyz (resíduos enfileirados ao longo de +X; formato XYZ simples)
    Retorna o caminho do CSV-mestre com o catálogo da exportação.
    """
    fmt = (fmt or "extxyz").lower().strip()
    os.makedirs(outdir, exist_ok=True)
    groups = _group_atoms_by_residue(result_atoms)

    def _write_residue_file(res_atoms: t.List[Atomo], path: str) -> None:
        if fmt == "xyz":
            _write_xyz_simple(res_atoms, path)
        elif fmt in ("annot", "annotated", "xyz_annot"):
            _write_xyz_annotated(res_atoms, path)
        else:
            _write_extxyz(res_atoms, path)

    rows: t.List[dict[str, t.Union[str, int]]] = []
    for (resn, chain, seq, icode), atoms in groups.items():
        atoms_sorted = sorted(atoms, key=lambda a: a.id)
        suffix = {"xyz": "xyz", "annot": "annot.xyz", "annotated": "annot.xyz", "xyz_annot": "annot.xyz"}.get(fmt, "extxyz")
        icode_tag = f"{icode}" if icode else ""
        filename = f"{pdb_code}_{chain}_{resn}{seq}{icode_tag}.{suffix}"
        path = os.path.join(outdir, filename)
        _write_residue_file(atoms_sorted, path)
        rows.append({
            "PDB": pdb_code,
            "Residue": resn,
            "Chain": chain,
            "Seq": int(seq),
            "ICode": icode,
            "Atoms": int(len(atoms_sorted)),
            "XYZ_File": path,
            "Stack_File": "",
        })

    stack_path = ""
    if make_stack and groups:
        stack_filename = f"STACK_{pdb_code}.xyz"
        stack_path = os.path.join(outdir, stack_filename)

        ordered = sorted(groups.items(), key=lambda kv: (kv[0][1], kv[0][2], kv[0][3], kv[0][0]))
        total_atoms = sum(len(v) for _, v in ordered)

        lines: t.List[str] = []
        offset_x = 0.0
        spacing = float(spacing) if spacing and spacing > 0 else 5.0

        for (resn, chain, seq, icode), atoms in ordered:
            atoms_sorted = sorted(atoms, key=lambda a: a.id)
            cx, cy, cz = _centroid(atoms_sorted)
            for a in atoms_sorted:
                elem = _infer_element(a.nome, a.tipo_atomo)
                x = a.x - cx + offset_x
                y = a.y - cy
                z = a.z - cz
                lines.append(f"{elem:2s} {x: .6f} {y: .6f} {z: .6f}\n")
            offset_x += spacing

        with open(stack_path, "w", encoding="utf-8") as f:
            f.write(f"{total_atoms}\n")
            f.write(f"Stacked residues along +X; spacing = {spacing:.3f} Å\n")
            for L in lines:
                f.write(L)

        for r in rows:
            r["Stack_File"] = stack_path

    csv_path = os.path.join(outdir, f"{pdb_code}_residue_xyz_index.csv")
    pd.DataFrame(rows).sort_values(by=["Chain", "Seq", "Residue"]).to_csv(csv_path, index=False)
    return csv_path

def subdividir_estrutura_em_grids(matriz_pdb: MatrizPDB, pdb_code: str) -> None:
    # Seleção de cadeia (opcional)
    while True:
        esc = _safe_input(translate("analyze_specific_chain")).strip().upper()
        if esc in ["S", "Y", "YES"]:
            cadeia_sel = _safe_input(translate("enter_chain_identifier")).strip().upper()
            matriz_pdb_filtrada = filtrar_por_cadeia(matriz_pdb, cadeia_sel)
            if not matriz_pdb_filtrada:
                print(translate("chain_not_found").format(cadeia=cadeia_sel))
                _safe_input(translate("press_enter"))
                continue
            matriz_pdb_para_grids = matriz_pdb_filtrada
            print(translate("chain_selected").format(cadeia=cadeia_sel))
            break
        elif esc in ["N", "NO", ""]:
            print(translate("analyzing_all_chains"))
            cadeia_sel = None
            matriz_pdb_para_grids = matriz_pdb
            break
        else:
            print(translate("invalid_option"))

    incluir_residuos = _yes_default(_safe_input(translate("grid_residue_mode_prompt")), default=False)

    dimensoes = calcular_dimensoes(matriz_pdb_para_grids)

    try:
        nx = int(_safe_input(translate("enter_num_grids_x")))
        ny = int(_safe_input(translate("enter_num_grids_y")))
        nz = int(_safe_input(translate("enter_num_grids_z")))
    except ValueError:
        print(translate("invalid_grid_values"))
        return

    if nx <= 0 or ny <= 0 or nz <= 0:
        print(translate("grid_numbers_must_be_positive"))
        return

    phys_stats: t.Optional[dict[tuple[int, int, int], dict[str, t.Union[str, float]]]] = None
    if _yes(_safe_input(translate("grid_physical_prompt"))):
        mode_in = _safe_input(translate("grid_phys_mode_prompt")).strip()
        if mode_in == "2":
            phys_mode = "sum"
        else:
            phys_mode = "union"

        pad = 0.0
        if _yes_default(_safe_input(translate("grid_expand_bounds_prompt")), default=True):
            extra_in = _safe_input(translate("grid_expand_amount_prompt")).strip()
            try:
                extra = float(extra_in) if extra_in else 0.0
            except ValueError:
                extra = 0.0
            pad = max(_max_vdw_radius(matriz_pdb_para_grids) + extra, 0.0)
            dimensoes = _expand_dims(dimensoes, pad)

        voxel_in = _safe_input(translate("grid_voxel_size_prompt")).strip()
        try:
            voxel = float(voxel_in) if voxel_in else 0.5
        except ValueError:
            voxel = 0.5
        if voxel <= 0:
            voxel = 0.5
        print(translate("grid_physical_note"))
        phys_stats = _compute_physical_grid_stats(matriz_pdb_para_grids, dimensoes, nx, ny, nz, voxel, phys_mode)

    grids = calcular_grids(dimensoes, nx, ny, nz)
    print(f"\n{translate('total_subboxes_generated')}: {len(grids)}\n")

    if not pymol_available:
        print(translate("pymol_visualization_unavailable"))
    else:
        try:
            finish_launching()
            cmd.fetch(pdb_code)
            _pymol_enable_click_info()
            _apply_pymol_style()

            if cadeia_sel:
                cmd.remove(f"not chain {cadeia_sel}")
                print(translate("visualizing_chain_in_pymol").format(cadeia=cadeia_sel))

            colors = ["red", "green", "yellow", "orange", "purple", "cyan", "magenta", "blue"]
            global _grid_atom_index, _grid_bounds, _grid_phys_stats
            _grid_atom_index = {}
            _grid_bounds = {}
            _grid_phys_stats = {}
            for g in grids:
                grid_name = f"grid_{g.i}_{g.j}_{g.k}"
                _grid_bounds[grid_name] = g
                if phys_stats and (g.i, g.j, g.k) in phys_stats:
                    _grid_phys_stats[f"{grid_name}_sel"] = phys_stats[(g.i, g.j, g.k)]
                cmd.load_cgo(create_cgo_box(g), grid_name)

                resultados, _ = buscar_aminoacidos_no_grid(matriz_pdb_para_grids, g, incluir_aminoacidos_completos=incluir_residuos)
                if resultados:
                    atom_ids = [a.id for a in resultados]
                    selection_str = "id " + "+".join(map(str, atom_ids))
                    sel_name = f"{grid_name}_sel"
                    cmd.select(sel_name, selection_str)
                    color = colors[(g.i + g.j + g.k) % len(colors)]
                    cmd.color(color, sel_name)
                    if phys_stats and (g.i, g.j, g.k) in phys_stats:
                        _grid_phys_stats[sel_name] = phys_stats[(g.i, g.j, g.k)]
                    for aid in atom_ids:
                        try:
                            _grid_atom_index.setdefault(int(aid), set()).add(sel_name)
                        except Exception:
                            pass

                    residuos_unicos = {(a.residuo, a.sequencia, a.cadeia) for a in resultados}
                    print(f"{translate('aminoacids_in_grid')} {grid_name}:")
                    for res, seq, ch in sorted(residuos_unicos, key=lambda x: x[1]):
                        nome_res = AMINOACID_CODES.get(res, res)
                        print(f"- {translate('residue')}: {nome_res} ({res}), {translate('sequence')}: {seq}, {translate('chain')}: {ch}")
                    print("")

            cmd.show("cgo", "grid_*")
            cmd.show("cartoon", "all")
            cmd.zoom("all")

            print(f"\n{translate('instructions_pymol_grid')}")
            print(translate("show_sticks_example"))
            print(translate("zoom_example"))
        except Exception as e:
            print(translate("pymol_visualization_unavailable"))
            logger.error(f"Error in PyMOL visualization: {e}")

    # -------- Inspect grid --------
    grid_inspect = _safe_input(translate("grid_inspect_prompt"))
    if grid_inspect.strip():
        try:
            gi, gj, gk = [int(x.strip()) for x in grid_inspect.split(",")]
            g_candidates = [g for g in grids if g.i == gi and g.j == gj and g.k == gk]
            if not g_candidates:
                print(translate("grid_not_found"))
            else:
                gsel = g_candidates[0]
                res_atoms, _ = buscar_aminoacidos_no_grid(matriz_pdb_para_grids, gsel, incluir_aminoacidos_completos=incluir_residuos)
                detail_in = _safe_input(translate("grid_inspect_depth")).strip()
                detailed = (detail_in == "2")
                phys = phys_stats.get((gi, gj, gk)) if phys_stats else None
                print(_grid_summary(res_atoms, detailed=detailed, phys=phys))
        except Exception as e:
            print(translate("grid_not_found"))
            logger.error(e)

    # -------- Export XYZ from a specific grid --------
    grid_str = _safe_input(translate("grid_export_which"))
    if grid_str.strip():
        try:
            gi, gj, gk = [int(x.strip()) for x in grid_str.split(",")]
            g_candidates = [g for g in grids if g.i == gi and g.j == gj and g.k == gk]
            if not g_candidates:
                print(translate("grid_not_found"))
                return
            gsel = g_candidates[0]
            res_atoms, _ = buscar_aminoacidos_no_grid(matriz_pdb_para_grids, gsel, incluir_aminoacidos_completos=incluir_residuos)
            if not res_atoms:
                print(translate("grid_not_found"))
                return

            if _yes(_safe_input(translate("export_xyz_prompt"))):
                fmt_in = _safe_input(translate("export_xyz_format")).strip()
                if fmt_in == "1":
                    fmt = "xyz"
                elif fmt_in == "2":
                    fmt = "annot"
                else:
                    fmt = "extxyz"

                outdir = _safe_input(translate("export_xyz_folder")).strip() or "xyz_exports"
                do_stack = _yes(_safe_input(translate("export_stack_prompt")))
                spacing_in = _safe_input(translate("export_stack_spacing")).strip()
                spacing = float(spacing_in) if spacing_in else 5.0

                csv_path = export_residues_to_xyz(res_atoms, pdb_code, outdir, do_stack, spacing, fmt)
                print(f"{translate('export_xyz_done')} {csv_path}")
        except Exception as e:
            print(translate("grid_not_found"))
            logger.error(e)

    # -------- Export grid stats --------
    if _yes(_safe_input(translate("grid_export_stats"))):
        stats_name = _safe_input(translate("grid_stats_filename")).strip() or "grid_stats.csv"
        rows: t.List[dict[str, t.Union[int, float, str]]] = []
        for g in grids:
            res_atoms, _ = buscar_aminoacidos_no_grid(matriz_pdb_para_grids, g, incluir_aminoacidos_completos=incluir_residuos)
            n_atoms = len(res_atoms)
            vol = max((g.nxmax - g.nxmin) * (g.nymax - g.nymin) * (g.nzmax - g.nzmin), 0.0)
            density = (n_atoms / vol) if vol > 0 else 0.0
            b_mean = float(np.mean([a.b_factor for a in res_atoms])) if n_atoms > 0 else 0.0

            row: dict[str, t.Union[int, float, str]] = {
                "i": g.i, "j": g.j, "k": g.k,
                "xmin": g.nxmin, "xmax": g.nxmax,
                "ymin": g.nymin, "ymax": g.nymax,
                "zmin": g.nzmin, "zmax": g.nzmax,
                "atoms": n_atoms,
                "volume_A3": float(vol),
                "density_atoms_per_A3": float(density),
                "bfactor_mean": float(b_mean),
            }
            if phys_stats:
                p: dict[str, t.Union[str, float]] = phys_stats.get((g.i, g.j, g.k), {})
                row["vdw_voxels"] = float(p.get("vdw_voxels", 0.0))
                row["vdw_volume_A3"] = float(p.get("vdw_volume_A3", 0.0))
                row["vdw_occupancy_fraction"] = float(p.get("vdw_occupancy_fraction", 0.0))
                if "vdw_mode" in p:
                    row["vdw_mode"] = str(p.get("vdw_mode"))
            rows.append(row)

        pd.DataFrame(rows).to_csv(stats_name, index=False)
        print(f"{translate('results_exported')} {stats_name}")

def exportar_resultados(
    resultados: t.Sequence[t.Union[Atomo, tuple[Atomo, Atomo, float]]],
    filename: str,
    formato: str = "csv",
    modo: str = "full",
) -> None:
    try:
        fmt = (formato or "").strip().lower()
        mode = (modo or "full").strip().lower()

        if fmt in ("csv", "tsv", "excel", "xlsx", "xls", "json"):
            rows = _build_export_rows(resultados, mode)
            if not rows:
                print(translate("no_results_found"))
                return
            df = pd.DataFrame(rows)

            if fmt in ("excel", "xlsx", "xls"):
                if not filename.lower().endswith((".xlsx", ".xls")):
                    filename += ".xlsx"
                df.to_excel(filename, index=False)
            elif fmt == "tsv":
                if not filename.lower().endswith(".tsv"):
                    filename += ".tsv"
                df.to_csv(filename, index=False, sep="\t")
            elif fmt == "json":
                if not filename.lower().endswith(".json"):
                    filename += ".json"
                df.to_json(filename, orient="records", force_ascii=False, indent=2)
            else:
                if not filename.lower().endswith(".csv"):
                    filename += ".csv"
                df.to_csv(filename, index=False)

        elif fmt in ("pdb", "cif"):
            atom_list = [a for a in resultados if isinstance(a, Atomo)]
            if len(atom_list) != len(resultados):
                print(translate("unsupported_format"))
                return
            if mode not in ("full", "compact", "xyz"):
                print(translate("unsupported_format"))
                return
            if fmt == "pdb":
                if not filename.lower().endswith(".pdb"):
                    filename += ".pdb"
                _export_atoms_pdb(atom_list, filename)
            else:
                if not filename.lower().endswith(".cif"):
                    filename += ".cif"
                _export_atoms_cif(atom_list, filename)
        else:
            print(translate("unsupported_format"))
            return

        print(f"{translate('results_exported')} {filename}")
    except (IOError, ValueError) as e:
        logger.error(f"{translate('error_exporting_results')}: {e}")
        print(f"{translate('error_exporting_results')}: {e}")

def _build_export_rows(
    resultados: t.Sequence[t.Union[Atomo, tuple[Atomo, Atomo, float]]],
    mode: str,
) -> t.List[dict[str, t.Union[int, float, str]]]:
    mode = (mode or "full").strip().lower()
    atoms = [a for a in resultados if isinstance(a, Atomo)]
    tuples_ = [a for a in resultados if not isinstance(a, Atomo)]

    if tuples_ and mode not in ("full", "interactions"):
        mode = "interactions"

    rows: t.List[dict[str, t.Union[int, float, str]]] = []

    if mode == "compact":
        for a in atoms:
            rows.append({
                "ID": a.id, "Átomo": a.nome, "Resíduo": a.residuo,
                "Cadeia": a.cadeia, "Sequência": a.sequencia,
                "X": a.x, "Y": a.y, "Z": a.z, "Elemento": _infer_element(a.nome, a.tipo_atomo)
            })

    elif mode == "composition":
        counts: dict[str, int] = {}
        for a in atoms:
            elem = _infer_element(a.nome, a.tipo_atomo)
            counts[elem] = counts.get(elem, 0) + 1
        for elem, n in sorted(counts.items()):
            rows.append({"Elemento": elem, "Contagem": int(n)})

    elif mode == "xyz":
        for a in atoms:
            elem = _infer_element(a.nome, a.tipo_atomo)
            rows.append({"Elemento": elem, "X": a.x, "Y": a.y, "Z": a.z})

    elif mode == "residue_summary":
        residue_groups: dict[tuple[str, str, int, str], t.List[Atomo]] = {}
        for a in atoms:
            key = (a.residuo, a.cadeia, a.sequencia, a.icode)
            residue_groups.setdefault(key, []).append(a)

        for (resn, chain, seq, icode), group in residue_groups.items():
            cx, cy, cz = _centroid(group)
            b_mean = sum(a.b_factor for a in group) / max(len(group), 1)
            xs = [a.x for a in group]
            ys = [a.y for a in group]
            zs = [a.z for a in group]
            b_std = float(np.std([a.b_factor for a in group])) if len(group) > 1 else 0.0

            rows.append({
                "Resíduo": resn, "Cadeia": chain, "Sequência": seq, "ICode": icode,
                "N_átomos": len(group), "Cx": cx, "Cy": cy, "Cz": cz,
                "X_min": min(xs), "X_max": max(xs),
                "Y_min": min(ys), "Y_max": max(ys),
                "Z_min": min(zs), "Z_max": max(zs),
                "B-factor_médio": float(b_mean),
                "B-factor_std": float(b_std),
            })

    elif mode == "chain_summary":
        chain_groups: dict[str, t.List[Atomo]] = {}
        for a in atoms:
            chain_groups.setdefault(a.cadeia, []).append(a)

        for chain, group in chain_groups.items():
            xs = [a.x for a in group]
            ys = [a.y for a in group]
            zs = [a.z for a in group]
            residues = {(a.residuo, a.sequencia, a.icode) for a in group}
            cx, cy, cz = _centroid(group)
            b_mean = sum(a.b_factor for a in group) / max(len(group), 1)
            b_std = float(np.std([a.b_factor for a in group])) if len(group) > 1 else 0.0

            rows.append({
                "Cadeia": chain,
                "N_átomos": len(group),
                "N_resíduos": len(residues),
                "X_min": min(xs), "X_max": max(xs),
                "Y_min": min(ys), "Y_max": max(ys),
                "Z_min": min(zs), "Z_max": max(zs),
                "Cx": cx, "Cy": cy, "Cz": cz,
                "B-factor_médio": float(b_mean),
                "B-factor_std": float(b_std),
            })

    elif mode == "interactions":
        for a in tuples_:
            a1, a2, d = t.cast(tuple[Atomo, Atomo, float], a)
            rows.append({
                "Atomo1_ID": a1.id, "Atomo1_Nome": a1.nome,
                "Atomo2_ID": a2.id, "Atomo2_Nome": a2.nome,
                "Distância": float(d),
            })

    else:  # full
        for a in atoms:
            rows.append({
                "ID": a.id,
                "Nome": a.nome,
                "Resíduo": AMINOACID_CODES.get(a.residuo, a.residuo),
                "Código Resíduo": a.residuo,
                "Cadeia": a.cadeia,
                "Sequência": a.sequencia,
                "ICode": a.icode,
                "AltLoc": a.altloc,
                "X": a.x, "Y": a.y, "Z": a.z,
                "B-factor": a.b_factor,
                "Ocupação": a.ocupacao,
                "HetFlag": a.hetflag,
                "Modelo": a.modelo,
                "Tipo Atomo": a.tipo_atomo,
                "Elemento": _infer_element(a.nome, a.tipo_atomo),
            })

    return rows

def export_last_grid_dump() -> None:
    if not _last_grid_dump:
        print(translate("export_last_grid_none"))
        return
    sel = _last_grid_dump.get("sel", "grid_sel")
    data = t.cast(t.List[dict[str, t.Any]], _last_grid_dump.get("data", []))
    if not data:
        print(translate("export_last_grid_none"))
        return

    fmt = _safe_input(translate("export_last_grid_format")).strip().lower() or "csv"
    default_name = f"{sel}_atoms.{fmt}"
    filename = _safe_input(translate("export_last_grid_filename")).strip() or default_name

    if fmt == "txt":
        lines = _format_grid_lines(sorted(data, key=_grid_sort_key))
        with open(filename, "w", encoding="utf-8") as f:
            f.write("\n".join(lines) + "\n")
        print(f"{translate('results_exported')} {filename}")
        return

    # default csv
    pd.DataFrame(data).to_csv(filename, index=False)
    print(f"{translate('results_exported')} {filename}")

def _parse_csv_list(s: str) -> t.List[str]:
    if not s:
        return []
    parts = [p.strip() for p in s.replace(";", ",").split(",")]
    return [p for p in parts if p]

def _parse_range(s: str) -> t.Optional[tuple[float, float]]:
    if not s:
        return None
    s = s.replace(" ", "")
    if "-" in s:
        a, b = s.split("-", 1)
    elif "," in s:
        a, b = s.split(",", 1)
    else:
        return None
    try:
        lo = float(a)
        hi = float(b)
    except ValueError:
        return None
    if lo > hi:
        lo, hi = hi, lo
    return (lo, hi)

def _parse_xyz_box(s: str) -> t.Optional[tuple[float, float, float, float, float, float]]:
    if not s:
        return None
    parts = [p.strip() for p in s.replace(";", ",").split(",")]
    if len(parts) != 6:
        return None
    try:
        vals = [float(p) for p in parts]
    except ValueError:
        return None
    return (vals[0], vals[1], vals[2], vals[3], vals[4], vals[5])

def _apply_export_filters(
    atoms: t.Sequence[Atomo],
    chains: set[str],
    residues: set[str],
    seq_range: t.Optional[tuple[float, float]],
    atom_names: set[str],
    elements: set[str],
    b_range: t.Optional[tuple[float, float]],
    xyz_box: t.Optional[tuple[float, float, float, float, float, float]],
) -> t.List[Atomo]:
    out: t.List[Atomo] = []
    for a in atoms:
        if chains and a.cadeia.upper() not in chains:
            continue
        if residues and a.residuo.upper() not in residues:
            continue
        if atom_names and a.nome.upper() not in atom_names:
            continue
        if elements:
            elem = _infer_element(a.nome, a.tipo_atomo).upper()
            if elem not in elements:
                continue
        if seq_range:
            lo, hi = seq_range
            if not (lo <= a.sequencia <= hi):
                continue
        if b_range:
            lo, hi = b_range
            if not (lo <= a.b_factor <= hi):
                continue
        if xyz_box:
            x0, x1, y0, y1, z0, z1 = xyz_box
            if not (x0 <= a.x <= x1 and y0 <= a.y <= y1 and z0 <= a.z <= z1):
                continue
        out.append(a)
    return out

def _format_pdb_atom_name(name: str, element: str) -> str:
    n = name.strip()
    elem = element.strip()
    if len(elem) == 1 and len(n) < 4:
        return f" {n:<3}"
    return f"{n:<4}"[:4]

def _export_atoms_pdb(atoms: t.Sequence[Atomo], filename: str) -> None:
    with open(filename, "w", encoding="utf-8") as f:
        for a in atoms:
            elem = _infer_element(a.nome, a.tipo_atomo).strip().rjust(2)
            name = _format_pdb_atom_name(a.nome, elem)
            resn = a.residuo.strip().rjust(3)[:3]
            chain = (a.cadeia or "A")[:1]
            seq = a.sequencia if a.sequencia is not None else 1
            icode = (a.icode or " ")[:1]
            altloc = (a.altloc or " ")[:1]
            occ = a.ocupacao if a.ocupacao is not None else 1.00
            b = a.b_factor if a.b_factor is not None else 0.00
            record = "HETATM" if (a.hetflag or " ").strip() not in ("", " ") else "ATOM  "
            f.write(
                f"{record}{a.id:5d} {name}{altloc}{resn} {chain}{seq:4d}{icode}"
                f"   {a.x:8.3f}{a.y:8.3f}{a.z:8.3f}{occ:6.2f}{b:6.2f}          {elem}\n"
            )
        f.write("END\n")

def _export_atoms_cif(atoms: t.Sequence[Atomo], filename: str) -> None:
    with open(filename, "w", encoding="utf-8") as f:
        f.write("data_export\n")
        f.write("loop_\n")
        f.write("_atom_site.group_PDB\n")
        f.write("_atom_site.id\n")
        f.write("_atom_site.type_symbol\n")
        f.write("_atom_site.label_atom_id\n")
        f.write("_atom_site.label_comp_id\n")
        f.write("_atom_site.label_asym_id\n")
        f.write("_atom_site.label_seq_id\n")
        f.write("_atom_site.label_alt_id\n")
        f.write("_atom_site.Cartn_x\n")
        f.write("_atom_site.Cartn_y\n")
        f.write("_atom_site.Cartn_z\n")
        f.write("_atom_site.occupancy\n")
        f.write("_atom_site.B_iso_or_equiv\n")
        for a in atoms:
            elem = _infer_element(a.nome, a.tipo_atomo)
            chain = (a.cadeia or "A")[:1]
            seq = a.sequencia if a.sequencia is not None else 1
            altloc = (a.altloc or ".")[:1]
            occ = a.ocupacao if a.ocupacao is not None else 1.0
            group = "HETATM" if (a.hetflag or " ").strip() not in ("", " ") else "ATOM"
            f.write(
                f"{group} {a.id} {elem} {a.nome} {a.residuo} {chain} {seq} {altloc} "
                f"{a.x:.3f} {a.y:.3f} {a.z:.3f} {occ:.2f} {a.b_factor:.2f}\n"
            )

def visualizar_estruturas_secundarias(pdb_code: str) -> None:
    if not pymol_available:
        print(translate("pymol_visualization_unavailable"))
        return
    try:
        finish_launching()
        cmd.fetch(pdb_code)
        _pymol_enable_click_info()
        _apply_pymol_style()
        cmd.show("cartoon", "all")
        cmd.color("gray", "all")
        cmd.select("helices", "ss h")
        cmd.color("red", "helices")
        cmd.select("sheets", "ss s")
        cmd.color("yellow", "sheets")
        cmd.select("turns", "ss l+")
        cmd.color("green", "turns")
        cmd.zoom("all")
        print(translate("secondary_structures_highlighted"))
    except Exception as e:
        print(translate("pymol_visualization_unavailable"))
        logger.error(f"Error in PyMOL visualization: {e}")

def visualizar_ligantes(pdb_code: str, ligantes: t.List[str]) -> None:
    if not ligantes:
        print(translate("no_ligands_found"))
        return
    if not pymol_available:
        print(translate("pymol_visualization_unavailable"))
        return
    try:
        finish_launching()
        cmd.fetch(pdb_code)
        _pymol_enable_click_info()
        _apply_pymol_style()
        cmd.show("cartoon", "polymer")
        cmd.color("gray", "polymer")
        lig_str = "+".join(ligantes)
        cmd.select("ligantes", f"resn {lig_str}")
        cmd.show("sticks", "ligantes")
        cmd.color("cyan", "ligantes")
        cmd.zoom("ligantes")
        print(translate("ligands_highlighted"))
    except Exception as e:
        print(translate("pymol_visualization_unavailable"))
        logger.error(f"Error in PyMOL visualization: {e}")

def colorir_por_bfactor(pdb_code: str) -> None:
    if not pymol_available:
        print(translate("pymol_visualization_unavailable"))
        return
    try:
        finish_launching()
        cmd.fetch(pdb_code)
        _pymol_enable_click_info()
        _apply_pymol_style()
        cmd.show("cartoon", "all")
        cmd.spectrum("b", "blue_white_red", "all")
        cmd.zoom("all")
        print(translate("structure_colored_by_bfactor"))
    except Exception as e:
        print(translate("pymol_visualization_unavailable"))
        logger.error(f"Error in PyMOL visualization: {e}")

def fechar_pymol() -> None:
    if not pymol_available:
        print(translate("pymol_visualization_unavailable"))
        return
    try:
        global _pymol_click_stop, _grid_atom_index, _grid_bounds, _grid_phys_stats, _last_grid_dump
        if _pymol_click_stop:
            _pymol_click_stop.set()
            _pymol_click_stop = None
        _grid_atom_index = {}
        _grid_bounds = {}
        _grid_phys_stats = {}
        _last_grid_dump = None
        # Hide the GUI and clear objects without terminating this process.
        try:
            cmd.do("window hide")
        except Exception:
            try:
                cmd.window("hide")  # type: ignore[call-arg]
            except Exception:
                pass
        try:
            cmd.delete("all")
        except Exception:
            pass
        try:
            cmd.reinitialize()
        except Exception:
            pass
        print(translate("pymol_closed"))
    except Exception as e:
        print(translate("pymol_close_failed"))
        logger.error(f"Error closing PyMOL: {e}")

def _toggle_click_mode() -> None:
    global CLICK_VERBOSE
    CLICK_VERBOSE = not CLICK_VERBOSE
    mode = translate('mode_detailed') if CLICK_VERBOSE else translate('mode_summary')
    print(translate('click_mode_now').format(modo=mode, mode=mode))

def _find_atom_by_spec(matriz_pdb: MatrizPDB, spec: str) -> t.Optional[Atomo]:
    """
    spec: "A,50,CA" (chain,residue,atom)
    Aceita espaços, insertion code (ex: 50A), e fallback para N,C,O se o átomo não existir.
    """
    try:
        parts = [p.strip() for p in spec.split(",")]
        if len(parts) != 3:
            return None
        chain, resi_raw, name = parts

        # parse resi + icode
        digits = "".join(ch for ch in resi_raw if ch.isdigit() or ch == "-")
        icode = "".join(ch for ch in resi_raw if ch.isalpha())
        resi = int(digits) if digits else None

        def match_atom(target_name: str) -> t.Optional[Atomo]:
            for a in matriz_pdb:
                if a.cadeia.upper() != chain.upper():
                    continue
                if resi is not None and a.sequencia != resi:
                    continue
                if icode and a.icode.upper() != icode.upper():
                    continue
                if a.nome.upper() == target_name.upper():
                    return a
            return None

        # 1) tenta nome solicitado
        hit = match_atom(name)
        if hit:
            return hit

        # 2) fallback se CA não existe: tenta backbone N,C,O
        backbone = ["N", "C", "O"]
        for alt in backbone:
            hit = match_atom(alt)
            if hit:
                return hit
        return None
    except Exception:
        return None

def _tunneling_path(
    matriz_pdb: MatrizPDB,
    start: Atomo,
    end: Atomo,
    cutoff: float = 4.5
) -> TunnelResult:
    """
    Grafo: nós = átomos; arestas se dist < cutoff.
    Peso = dist * fator_elemento (metais=0.5, aromáticos=0.8, padrão=1.0)
    Dijkstra com registro de hops.
    """
    try:
        coords = np.array([(a.x, a.y, a.z) for a in matriz_pdb], dtype=float)
        tree = KDTree(coords)
    except Exception:
        return TunnelResult(False, "KDTree falhou", "", "", cutoff, 0.0, 0.0, 0, 0, [], [])

    idx_start = matriz_pdb.index(start)
    idx_end = matriz_pdb.index(end)

    elem_factor: dict[str, float] = {
        "CU": 0.5, "FE": 0.5, "ZN": 0.6, "NI": 0.6, "MN": 0.6,
        "C": 1.0, "N": 1.0, "O": 1.0, "S": 0.9
    }
    aromatic_res = {"PHE", "TYR", "TRP", "HIS"}

    def node_weight(idx: int) -> float:
        a = matriz_pdb[idx]
        elem = _infer_element(a.nome, a.tipo_atomo).upper()
        base = elem_factor.get(elem, 1.0)
        if a.residuo.upper() in aromatic_res:
            base *= 0.8
        return base

    n = len(matriz_pdb)
    dist_arr = np.full(n, np.inf)
    prev: list[int | None] = [None] * n
    prev_cost = np.full(n, np.inf)
    prev_dist = np.full(n, np.inf)
    prev_factor = np.full(n, np.inf)
    visited = np.zeros(n, dtype=bool)

    dist_arr[idx_start] = 0.0
    import heapq
    pq: list[tuple[float, int]] = [(0.0, idx_start)]

    while pq:
        d, u = heapq.heappop(pq)
        if visited[u]:
            continue
        visited[u] = True
        if u == idx_end:
            break
        neigh = tree.query_ball_point(coords[u], cutoff)
        w_u = node_weight(u)
        for v in neigh:
            if v == u:
                continue
            dv = float(np.linalg.norm(coords[u] - coords[v]))
            w_v = node_weight(v)
            factor = 0.5 * (w_u + w_v)
            w = dv * factor
            nd = d + w
            if nd < dist_arr[v]:
                dist_arr[v] = nd
                prev[v] = u
                prev_cost[v] = w
                prev_dist[v] = dv
                prev_factor[v] = factor
                heapq.heappush(pq, (nd, v))

    direct = start.distancia_ate(end)
    if not np.isfinite(dist_arr[idx_end]):
        return TunnelResult(
            False,
            "Sem caminho dentro do cutoff",
            _spec_of(start),
            _spec_of(end),
            cutoff,
            direct,
            np.inf,
            n,
            0,
            [],
            []
        )

    path_idx: t.List[int] = []
    cur = idx_end
    while cur is not None:
        path_idx.append(cur)
        cur = prev[cur]
    path_idx.reverse()

    hops: list[TunnelHop] = []
    cumulative = 0.0
    for i in range(len(path_idx) - 1):
        a_idx = path_idx[i]
        b_idx = path_idx[i + 1]
        dist_hop = float(prev_dist[b_idx])
        cost_hop = float(prev_cost[b_idx])
        factor = float(prev_factor[b_idx])
        cumulative += cost_hop
        a = matriz_pdb[a_idx]
        b = matriz_pdb[b_idx]
        hops.append(
            TunnelHop(
                i=i,
                a_id=a.id,
                b_id=b.id,
                a_spec=_spec_of(a),
                b_spec=_spec_of(b),
                dist_A=dist_hop,
                factor=factor,
                cost=cost_hop,
                cumulative=cumulative,
            )
        )

    nodes = [matriz_pdb[i] for i in path_idx]
    return TunnelResult(
        True,
        "ok",
        _spec_of(start),
        _spec_of(end),
        cutoff,
        direct,
        cumulative,
        len(nodes),
        len(hops),
        nodes,
        hops,
    )

def _prompt_tunneling(matriz_pdb: MatrizPDB) -> None:
    a_spec = _safe_input(translate("tunnel_prompt_start")).strip()
    b_spec = _safe_input(translate("tunnel_prompt_end")).strip()
    a_atom = _find_atom_by_spec(matriz_pdb, a_spec)
    b_atom = _find_atom_by_spec(matriz_pdb, b_spec)
    if not a_atom or not b_atom:
        print(translate("tunnel_not_found"))
        return
    res = _tunneling_path(matriz_pdb, a_atom, b_atom)
    if not res.ok:
        print(res.reason)
        return
    print(translate("tunnel_path_found").format(steps=len(res.nodes), weight=res.total_cost))
    for hop in res.hops:
        print(
            f"hop {hop.i}: {hop.a_spec} -> {hop.b_spec} "
            f"dist={hop.dist_A:.2f}Å factor={hop.factor:.2f} "
            f"cost={hop.cost:.2f} cum={hop.cumulative:.2f}"
        )
    print("Nós (ordem):")
    for a in res.nodes:
        print(f"{_spec_of(a)} xyz=({a.x:.3f},{a.y:.3f},{a.z:.3f})")
def limpar_tela() -> None:
    try:
        if os.name == "nt":
            os.system("cls")
        else:
            # ANSI clear screen + scrollback (VSCode terminal keeps scrollback otherwise)
            print("\033[2J\033[H\033[3J", end="")
    except Exception:
        os.system("cls" if os.name == "nt" else "clear")

def exibir_header(header: t.List[str]) -> None:
    limpar_tela()
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(APP_FULL_NAME)
    if header:
        print("\n".join(header))
    else:
        print(translate("no_header"))
    print(f"Timestamp: {timestamp}")

def menu_principal() -> None:
    global language, AMINOACID_CODES

    limpar_tela()
    print(translate("language_choice"))
    print(translate("language_options"))
    lang_choice = _safe_input(translate("choice")).strip()

    if lang_choice == "2":
        language = "en"
        AMINOACID_CODES = AMINOACID_CODES_EN
    else:
        language = "pt"
        AMINOACID_CODES = AMINOACID_CODES_PT

    print(translate("welcome"))

    while True:
        limpar_tela()
        print(translate("start_mode"))
        print(translate("start_mode_options"))
        modo = _safe_input(translate("choice")).strip()

        pdb_content = ""
        formato = ""
        pdb_code = ""

        if modo == "2":
            while True:
                path = selecionar_arquivo_local_gui()
                if not path:
                    path = _safe_input(translate("enter_local_path")).strip()
                if _back(path):
                    break
                loaded = carregar_arquivo_local(path)
                if loaded is None:
                    continue
                pdb_content, formato, pdb_code = loaded
                break
            if _back(path):
                continue
        else:
            pdb_code = _safe_input(translate("enter_pdb_code")).strip().upper()
            if (pdb_code == "SAIR" and language == "pt") or (pdb_code == "EXIT" and language == "en"):
                print(translate("exit_message"))
                return

            resultado_download = download_pdb(pdb_code)
            if resultado_download is None:
                print(translate("invalid_pdb"))
                continue
            pdb_content, formato = resultado_download

        estrutura = ler_estrutura_pdb(pdb_content, formato, pdb_code)
        if estrutura is None:
            print(translate("fail_read_structure"))
            continue

        header = ler_cabecalho_de_estrutura(estrutura)
        matriz_pdb = converter_estrutura_para_atomos(estrutura)
        if not matriz_pdb:
            print(translate("fail_extract_atoms"))
            continue

        ligantes = buscar_ligantes(estrutura)
        historico_acoes: t.List[str] = [f"{translate('pdb_loaded')} {pdb_code}."]

        while True:
            exibir_header(header)
            escolha = _safe_input(translate("menu_options")).strip()

            if escolha == "1":
                while True:
                    exibir_header(header)
                    criterio = _safe_input(translate("enter_search_field")).strip()
                    if _back(criterio):
                        break
                    valor_raw = _safe_input(translate("enter_value_for").format(criterio=criterio))
                    valor: t.Union[int, float, str] = valor_raw

                    try:
                        fieldmap = FIELD_MAP_EN if language == "en" else FIELD_MAP_PT
                        attr = fieldmap.get(criterio.lower(), criterio)
                        if attr in ["x", "y", "z"]:
                            valor = float(valor_raw)
                        elif attr in ["id", "sequencia", "sequence"]:
                            valor = int(valor_raw)
                    except ValueError:
                        print(translate("invalid_value_for_criterion"))
                        continue

                    resultados = buscar_por_criterio(matriz_pdb, criterio, valor)
                    limpar_tela()
                    exibir_header(header)
                    imprimir_resultados(resultados)
                    historico_acoes.append(translate("search_by_criterion_done").format(criterio=criterio, valor=valor))
                    _safe_input(translate("press_enter"))
                    if not _repeat_prompt():
                        break

            elif escolha == "2":
                while True:
                    exibir_header(header)
                    eixo_in = _safe_input(translate("enter_axis")).strip()
                    if _back(eixo_in):
                        break
                    eixo = _validate_axis(eixo_in)
                    if eixo is None:
                        print(translate("invalid_attribute"))
                        continue
                    try:
                        inicio = float(_safe_input(translate("enter_start_value")))
                        fim = float(_safe_input(translate("enter_end_value")))
                    except ValueError:
                        print(translate("enter_numeric_values"))
                        continue

                    resultados = buscar_por_intervalo(matriz_pdb, eixo, inicio, fim)
                    limpar_tela()
                    exibir_header(header)
                    imprimir_resultados(resultados)
                    historico_acoes.append(translate("search_by_range_done").format(eixo=eixo, inicio=inicio, fim=fim))
                    _safe_input(translate("press_enter"))
                    if not _repeat_prompt():
                        break

            elif escolha == "3":
                while True:
                    exibir_header(header)
                    dimensoes = calcular_dimensoes(matriz_pdb)
                    imprimir_dimensoes(dimensoes)
                    historico_acoes.append(translate("structure_dimensions_calculated"))
                    _safe_input(translate("press_enter"))
                    if not _repeat_prompt():
                        break

            elif escolha == "4":
                while True:
                    exibir_header(header)
                    dist_in = _safe_input(translate("enter_max_distance")).strip()
                    if _back(dist_in):
                        break
                    try:
                        dmax = float(dist_in)
                    except ValueError:
                        print(translate("enter_numeric_value_for_distance"))
                        continue

                    ignore_h = _yes(_safe_input(translate("distance_ignore_h")))
                    exclude_same = _yes(_safe_input(translate("distance_exclude_same_residue")))
                    resultados = identificar_interacoes(matriz_pdb, dmax, ignore_h, exclude_same)

                    limpar_tela()
                    exibir_header(header)
                    imprimir_resultados(resultados)
                    historico_acoes.append(translate("search_by_max_distance_done").format(distancia=dmax))
                    _safe_input(translate("press_enter"))
                    if not _repeat_prompt():
                        break

            elif escolha == "5":
                while True:
                    exibir_header(header)
                    try:
                        ix = float(_safe_input(translate("enter_start_x")))
                        fx = float(_safe_input(translate("enter_end_x")))
                        iy = float(_safe_input(translate("enter_start_y")))
                        fy = float(_safe_input(translate("enter_end_y")))
                        iz = float(_safe_input(translate("enter_start_z")))
                        fz = float(_safe_input(translate("enter_end_z")))
                    except ValueError:
                        print(translate("enter_numeric_values"))
                        continue

                    incluir = _yes(_safe_input(translate("include_complete_aminoacids")))
                    resultados, _ = buscar_por_intervalo_simultaneo(matriz_pdb, (ix, iy, iz), (fx, fy, fz), incluir)

                    limpar_tela()
                    exibir_header(header)
                    imprimir_resultados(resultados)

                    if resultados and _yes(_safe_input(translate("export_xyz_prompt"))):
                        fmt_in = _safe_input(translate("export_xyz_format")).strip()
                        if fmt_in == "1":
                            fmt = "xyz"
                        elif fmt_in == "2":
                            fmt = "annot"
                        else:
                            fmt = "extxyz"
                        outdir = _safe_input(translate("export_xyz_folder")).strip() or "xyz_exports"
                        do_stack = _yes(_safe_input(translate("export_stack_prompt")))
                        spacing_in = _safe_input(translate("export_stack_spacing")).strip()
                        spacing = float(spacing_in) if spacing_in else 5.0
                        csv_path = export_residues_to_xyz(resultados, pdb_code, outdir, do_stack, spacing, fmt)
                        print(f"{translate('export_xyz_done')} {csv_path}")

                    historico_acoes.append(translate("simultaneous_range_search_done").format(incluir=incluir))
                    _safe_input(translate("press_enter"))
                    if not _repeat_prompt():
                        break

            elif escolha == "6":
                while True:
                    exibir_header(header)
                    subdividir_estrutura_em_grids(matriz_pdb, pdb_code)
                    historico_acoes.append(translate("structure_subdivided_into_grids"))
                    _safe_input(translate("press_enter"))
                    if not _repeat_prompt():
                        break

            elif escolha == "7":
                while True:
                    exibir_header(header)
                    salvar_arquivo(pdb_code, pdb_content, formato)
                    historico_acoes.append(translate("file_saved_locally").format(filename=f"{pdb_code}.{formato}"))
                    _safe_input(translate("press_enter"))
                    if not _repeat_prompt():
                        break

            elif escolha == "8":
                while True:
                    exibir_header(header)
                    print(f"\n{translate('action_history')}:")
                    for acao in historico_acoes:
                        print(f"- {acao}")
                    _safe_input(translate("press_enter"))
                    if not _repeat_prompt():
                        break

            elif escolha == "9":
                while True:
                    exibir_header(header)
                    filename = _safe_input(translate("enter_filename_for_export"))
                    mode_choice = _safe_input(translate("export_mode_prompt")).strip()
                    mode_map = {
                        "1": "full",
                        "2": "compact",
                        "3": "residue_summary",
                        "4": "chain_summary",
                        "5": "composition",
                        "6": "xyz",
                    }
                    export_mode = mode_map.get(mode_choice, "full")
                    formato_export = _safe_input(translate("enter_export_format")).strip().lower()

                    export_atoms = list(matriz_pdb)
                    if _yes(_safe_input(translate("export_filter_prompt"))):
                        chains = {c.strip().upper() for c in _parse_csv_list(_safe_input(translate("export_filter_chain")))}
                        residues = {r.strip().upper() for r in _parse_csv_list(_safe_input(translate("export_filter_residue")))}
                        seq_range = _parse_range(_safe_input(translate("export_filter_seq")).strip())
                        atom_names = {n.strip().upper() for n in _parse_csv_list(_safe_input(translate("export_filter_atom")))}
                        elements = {e.strip().upper() for e in _parse_csv_list(_safe_input(translate("export_filter_element")))}
                        b_range = _parse_range(_safe_input(translate("export_filter_bfactor")).strip())
                        xyz_box = _parse_xyz_box(_safe_input(translate("export_filter_coords")).strip())

                        export_atoms = _apply_export_filters(export_atoms, chains, residues, seq_range, atom_names, elements, b_range, xyz_box)
                        print(translate("export_filter_applied") + str(len(export_atoms)))

                    exportar_resultados(export_atoms, filename, formato_export, export_mode)
                    historico_acoes.append(translate("results_exported_to_file").format(filename=filename))
                    _safe_input(translate("press_enter"))
                    if not _repeat_prompt():
                        break

            elif escolha == "10":
                while True:
                    exibir_header(header)
                    print(f"\n{translate('about_project')}:")
                    print(translate("project_info"))
                    _safe_input(translate("press_enter"))
                    if not _repeat_prompt():
                        break

            elif escolha == "11":
                while True:
                    visualizar_estruturas_secundarias(pdb_code)
                    historico_acoes.append(translate("secondary_structures_visualized"))
                    _safe_input(translate("press_enter"))
                    if not _repeat_prompt():
                        break

            elif escolha == "12":
                while True:
                    visualizar_ligantes(pdb_code, ligantes)
                    historico_acoes.append(translate("ligands_visualized"))
                    _safe_input(translate("press_enter"))
                    if not _repeat_prompt():
                        break

            elif escolha == "13":
                while True:
                    colorir_por_bfactor(pdb_code)
                    historico_acoes.append(translate("structure_colored_by_bfactor_done"))
                    _safe_input(translate("press_enter"))
                    if not _repeat_prompt():
                        break

            elif escolha == "14":
                limpar_tela()
                print(translate("returning_to_start"))
                break

            elif escolha == "15":
                fechar_pymol()
                _safe_input(translate("press_enter"))

            elif escolha == "16":
                while True:
                    exibir_header(header)
                    export_last_grid_dump()
                    _safe_input(translate("press_enter"))
                    if not _repeat_prompt():
                        break

            elif escolha == "17":
                while True:
                    limpar_tela()
                    _toggle_click_mode()
                    _safe_input(translate("press_enter"))
                    if not _repeat_prompt():
                        break

            elif escolha == "18":
                while True:
                    exibir_header(header)
                    _prompt_tunneling(matriz_pdb)
                    _safe_input(translate("press_enter"))
                    if not _repeat_prompt():
                        break

            elif escolha == "19":
                print(translate("exit_message"))
                return

            else:
                print(translate("invalid_option"))

if __name__ == "__main__":
    menu_principal()
