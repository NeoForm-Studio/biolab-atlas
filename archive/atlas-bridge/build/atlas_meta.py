# -*- coding: utf-8 -*-
"""
atlas_meta.py  —  Curated knowledge for BioLab Atlas. EDIT THIS to add genomes, genes,
expression, or phenotype, then run `python build_atlas.py` to rebuild ../index.html.

INTEGRITY RULES (do not break — see CONTEXT.md):
  - expr values must be REAL DESeq2 output from a REAL GEO series. Never invent them.
  - PHENO values must come from published literature (cited in CONTEXT.md). Never estimate.
  - Sequences / GC / identity / presence are computed from RefSeq by build_atlas.py, not hand-entered.
"""

SPECIES = {'DR1': {'name': 'D. radiodurans',
         'strain': 'R1',
         'kind': 'radiotolerant',
         'acc': 'GCF_000008565.1',
         'lin': 'Deinococcus–Thermus',
         'dom': 'Bacteria'},
 'DDES': {'name': 'D. deserti',
          'strain': 'VCD115',
          'kind': 'radiotolerant',
          'acc': 'GCF_000020685.1',
          'lin': 'Deinococcus–Thermus',
          'dom': 'Bacteria'},
 'DGOB': {'name': 'D. gobiensis',
          'strain': 'I-0',
          'kind': 'radiotolerant',
          'acc': 'GCF_000252445.1',
          'lin': 'Deinococcus–Thermus',
          'dom': 'Bacteria'},
 'DGEO': {'name': 'D. geothermalis',
          'strain': 'DSM 11300',
          'kind': 'radiotolerant',
          'acc': 'GCF_000196275.1',
          'lin': 'Deinococcus–Thermus',
          'dom': 'Bacteria'},
 'DPRO': {'name': 'D. proteolyticus',
          'strain': 'MRP',
          'kind': 'radiotolerant',
          'acc': 'GCF_000190555.1',
          'lin': 'Deinococcus–Thermus',
          'dom': 'Bacteria'},
 'TTH': {'name': 'T. thermophilus',
         'strain': 'HB8',
         'kind': 'sensitive control',
         'acc': 'GCF_000091545.1',
         'lin': 'Deinococcus–Thermus',
         'dom': 'Bacteria'},
 'RXY': {'name': 'R. xylanophilus',
         'strain': 'DSM 9941',
         'kind': 'radiotolerant',
         'acc': 'GCF_019448335.1',
         'lin': 'Actinobacteria',
         'dom': 'Bacteria'},
 'KRAD': {'name': 'K. radiotolerans',
          'strain': 'SRS30216',
          'kind': 'radiotolerant',
          'acc': 'GCF_000017305.1',
          'lin': 'Actinobacteria',
          'dom': 'Bacteria'},
 'CHRO': {'name': 'Chroococcidiopsis',
          'strain': 'PCC 7203',
          'kind': 'radiotolerant',
          'acc': 'GCF_000317125.1',
          'lin': 'Cyanobacteria',
          'dom': 'Bacteria'},
 'HYM': {'name': 'H. yonginensis',
         'strain': '—',
         'kind': 'radiotolerant',
         'acc': 'GCF_027625995.1',
         'lin': 'Bacteroidetes',
         'dom': 'Bacteria'},
 'TGAM': {'name': 'T. gammatolerans',
          'strain': 'EJ3',
          'kind': 'radiotolerant',
          'acc': 'GCF_000022365.1',
          'lin': 'Euryarchaeota',
          'dom': 'Archaea',
          'distant': True},
 'HALO': {'name': 'H. salinarum',
          'strain': 'NRC-1',
          'kind': 'radiotolerant',
          'acc': 'GCF_000006805.1',
          'lin': 'Euryarchaeota',
          'dom': 'Archaea',
          'distant': True},
 'PFUR': {'name': 'P. furiosus',
          'strain': 'DSM 3638',
          'kind': 'radiotolerant',
          'acc': 'GCF_008245085.1',
          'lin': 'Euryarchaeota',
          'dom': 'Archaea',
          'distant': True},
 'SHEW': {'name': 'S. oneidensis',
          'strain': 'MR-1',
          'kind': 'sensitive control',
          'acc': 'GCF_000146165.2',
          'lin': 'Proteobacteria',
          'dom': 'Bacteria'},
 'EC': {'name': 'E. coli',
        'strain': 'K-12',
        'kind': 'sensitive control',
        'acc': 'GCF_000005845.2',
        'lin': 'Proteobacteria',
        'dom': 'Bacteria'}}

# Order genes appear in the UI
ANCHORS = ['recA',
 'ddrA',
 'ddrB',
 'pprA',
 'uvrA',
 'uvrB',
 'ddrC',
 'ddrD',
 'irrE',
 'recN',
 'recF',
 'recQ',
 'ruvB',
 'mutS',
 'gyrA',
 'polA',
 'lexA',
 'tkt']

# Gene-symbol synonyms across annotations
SYN = {'irrE': ['irrE', 'pprI'], 'tkt': ['tkt', 'tktA']}

# Measured radiation survival (published literature; cite in CONTEXT.md)
PHENO = {'DR1': {'short': '~5–12 kGy',
         'full': 'Survives 5,000 Gy with almost no loss of viability; ~10% survive 12,000 Gy. '
                 'Roughly 50× more radiation-resistant than E. coli, ~1,000× more than humans.'},
 'DDES': {'short': 'up to ~15 kGy',
          'full': 'Isolated from the Sahara after a 15 kGy gamma blast; tolerates acute doses '
                  'around 15,000 Gy, with strong UV and desiccation resistance.'},
 'DGOB': {'short': '>15 kGy',
          'full': 'Resistant to >15,000 Gy gamma and >600 J/m² UV — among the most '
                  'radiation-tolerant organisms known. Nearly half its genome is '
                  'unknown-function.'},
 'DGEO': {'short': 'multi-kGy',
          'full': 'A thermophilic, radiation- and desiccation-resistant Deinococcus from hot '
                  'springs; survives several kGy of ionizing radiation, in the Deinococcus '
                  'high-tolerance range.'},
 'DPRO': {'short': 'multi-kGy',
          'full': 'A radiation-resistant Deinococcus (strain MRP); survives several kGy of '
                  'ionizing radiation, in the Deinococcus high-tolerance range.'},
 'TTH': {'short': 'lethal by a few hundred Gy',
         'full': 'A Deinococcus–Thermus relative that is radiation-SENSITIVE — killed by a few '
                 'hundred Gy, like E. coli. An in-phylum negative control: close relatives need '
                 'not be resistant.'},
 'RXY': {'short': '>25 kGy',
         'full': 'A thermophilic actinobacterium among the most radioresistant organisms known '
                 '(survives >25,000 Gy); pairs trehalose with manganese for antioxidant '
                 'protection.'},
 'KRAD': {'short': '~10–12 kGy (D10)',
          'full': 'An actinobacterium with γ-radiation D10 ≈ 10–12 kGy, approaching D. radiodurans '
                  '— radioresistance in a separate phylum.'},
 'CHRO': {'short': '~15 kGy',
          'full': 'A desiccation-tolerant cyanobacterium that survives ~15,000 Gy ionizing '
                  'radiation and high UV — photosynthetic radioresistance.'},
 'HYM': {'short': '~10–30 kGy (genus)',
         'full': 'A red-pigmented Bacteroidetes bacterium; the genus is recovered from desert '
                 'soils after 10–30 kGy gamma — radioresistance in yet another phylum, far from '
                 'Deinococcus.'},
 'TGAM': {'short': 'up to ~30 kGy',
          'full': 'A hyperthermophilic archaeon often called the most radiation-resistant organism '
                  'known: 5,000 Gy with no loss of viability, D10 ≈ 7.5 kGy, surviving up to '
                  '30,000 Gy. NOTE: archaea name most repair genes differently, so its by-symbol '
                  'toolkit count here is a severe undercount.'},
 'HALO': {'short': '~5 kGy (D10)',
          'full': 'A haloarchaeon with D10 ≈ 5 kGy and a high manganese/iron ratio like D. '
                  'radiodurans; resists desiccation, UV and IR. NOTE: archaeal repair genes are '
                  'differently named — the by-symbol count undercounts it.'},
 'PFUR': {'short': 'multi-kGy',
          'full': 'A hyperthermophilic archaeon with radioresistant behaviour (survives multi-kGy '
                  'IR), below T. gammatolerans. NOTE: archaeal repair genes are differently named '
                  '— the by-symbol count undercounts it.'},
 'SHEW': {'short': 'lethal at low doses',
          'full': 'Radiation-SENSITIVE reference (strain MR-1): high intracellular iron and low '
                  'manganese — the opposite chemistry of the radioresistant organisms, a clean '
                  'test of the manganese-antioxidant model.'},
 'EC': {'short': 'lethal by a few hundred Gy',
        'full': 'Radiation-sensitive reference: ~10% survive only 150 Gy, and 200–800 Gy is '
                'lethal. (For scale, 5–10 Gy kills a human.)'}}

# Curated per-gene knowledge. expr = REAL DESeq2 (DR1 / GSE176207) or None.
META = {'recA': {'scope': 'cons',
          'pathway': 'Homologous recombination',
          'plain': 'An enzyme that searches the genome for matching DNA and swaps strands between '
                   'the pieces — the search-and-match engine of repair.',
          'why': 'After radiation shatters the chromosome into hundreds of fragments, RecA pairs '
                 'each broken end with an intact matching copy so the genome can be rebuilt '
                 'accurately. Central to how D. radiodurans survives doses lethal to most life.',
          'tags': ['homologous recombination', 'double-strand break repair', 'strand exchange'],
          'role': 'RecA recombinase — ATP-dependent homologous strand exchange.',
          'expr': {'fc': 2.13, 'padj': 0.36, 'base': 27.7}},
 'ddrA': {'scope': 'deino',
          'pathway': 'DNA-end protection',
          'plain': 'A protective protein that caps the loose single-stranded ends of broken DNA so '
                   "they aren't chewed away.",
          'why': 'Right after irradiation it shields exposed DNA ends, buying time for repair. One '
                 "of the first Deinococcus 'mystery genes' shown by cross-species comparison to be "
                 'a Rad52-like repair factor.',
          'tags': ['double-strand break repair', 'single-strand DNA end protection'],
          'role': 'DdrA — Rad52-like 3′ ssDNA-end protection.',
          'expr': {'fc': 3.08, 'padj': 0.087, 'base': 38.3}},
 'ddrB': {'scope': 'deino',
          'pathway': 'Single-strand annealing',
          'plain': 'A protein that coats single-stranded DNA and helps matching pieces find each '
                   'other and rejoin.',
          'why': 'Drives a fast, early route to rejoin broken DNA (single-strand annealing) '
                 'especially important in Deinococcus.',
          'tags': ['single-strand DNA binding',
                   'single-strand annealing',
                   'double-strand break repair'],
          'role': 'DdrB — ssDNA-binding, single-strand annealing.',
          'expr': {'fc': 2.61, 'padj': 0.2, 'base': 153.3}},
 'pprA': {'scope': 'deino',
          'pathway': 'Break-end binding / ligation',
          'plain': 'A protein that grabs the two ends of a double-strand break and helps seal them '
                   'shut.',
          'why': 'Stimulates the DNA-ligase step that re-seals breaks; a signature Deinococcus '
                 'radiation factor — removing it sharply lowers radiation survival.',
          'tags': ['double-strand break repair', 'DNA ligation'],
          'role': 'PprA — binds DSB ends, stimulates ligase.',
          'expr': {'fc': 2.86, 'padj': 0.23, 'base': 31.5}},
 'uvrA': {'scope': 'cons',
          'pathway': 'Nucleotide excision repair',
          'plain': 'The damage scout of a repair crew that cuts bulky lesions out of DNA.',
          'why': 'Recognizes UV- and chemical-damaged spots and flags them for removal — first '
                 'step of nucleotide excision repair.',
          'tags': ['nucleotide excision repair', 'damage recognition'],
          'role': 'UvrA — UvrABC damage recognition (2 paralogs).',
          'expr': {'fc': 1.84,
                   'padj': 0.39,
                   'base': 61.7,
                   'note': 'two paralogs in R1 (+1.84, +2.39)'}},
 'uvrB': {'scope': 'cons',
          'pathway': 'Nucleotide excision repair',
          'plain': 'The verifier and unwinding partner of UvrA that confirms a lesion and sets up '
                   'the cut.',
          'why': 'Works with UvrA to excise damaged DNA so it can be resynthesized cleanly.',
          'tags': ['nucleotide excision repair', 'helicase'],
          'role': 'UvrB — UvrABC helicase / lesion verification.',
          'expr': {'fc': 2.43, 'padj': 0.16, 'base': 57.6}},
 'ddrC': {'scope': 'deino',
          'pathway': 'Nucleoid protection',
          'plain': 'A protein that packs and protects the genome while the cell is recovering.',
          'why': 'Helps condense and shield DNA during repair, limiting further damage.',
          'tags': ['nucleoid protection', 'genome organization'],
          'role': 'DdrC — nucleoid-associated genome protection.',
          'expr': None},
 'ddrD': {'scope': 'deino',
          'pathway': 'Radiation-response (uncharacterized)',
          'plain': 'A strongly radiation-induced protein whose exact job is still unknown.',
          'why': 'Switches on hard after irradiation, making it a prime candidate for an '
                 'uncharacterized repair function — the kind this atlas exists to investigate.',
          'tags': ['radiation-induced', 'uncharacterized'],
          'role': 'DdrD — RDR-induced; role being characterized.',
          'expr': None},
 'irrE': {'scope': 'deino',
          'pathway': 'Master regulation',
          'plain': 'The master switch that turns the whole radiation-response program on.',
          'why': 'After damage, IrrE cuts a repressor (DdrO), releasing the entire '
                 'radiation/desiccation gene set at once. A top-level controller found in '
                 'Deinococcus and relatives.',
          'tags': ['regulation', 'master switch', 'protease'],
          'role': 'IrrE/PprI — protease master regulator.',
          'expr': None},
 'recN': {'scope': 'cons',
          'pathway': 'Break-end tethering',
          'plain': "A protein that holds the two ends of a break together so they don't drift "
                   'apart.',
          'why': 'Tethers and aligns double-strand break ends so repair enzymes can rejoin them.',
          'tags': ['double-strand break repair', 'DNA end tethering'],
          'role': 'RecN — SMC-like DNA-end tethering.',
          'expr': None},
 'recF': {'scope': 'cons',
          'pathway': 'RecA loading (RecFOR)',
          'plain': 'A helper that loads the RecA repair enzyme onto DNA at the right spot.',
          'why': 'Part of the RecFOR pathway that starts recombination at single-strand gaps.',
          'tags': ['homologous recombination', 'RecA loading'],
          'role': 'RecF — RecFOR-pathway RecA loading.',
          'expr': None},
 'recQ': {'scope': 'cons',
          'pathway': 'Recombination processing',
          'plain': 'A helicase — a motor that unwinds the DNA double helix.',
          'why': 'Opens and untangles DNA and repair intermediates so they can be processed '
                 'correctly.',
          'tags': ['helicase', 'homologous recombination'],
          'role': 'RecQ — helicase, intermediate processing.',
          'expr': None},
 'ruvB': {'scope': 'cons',
          'pathway': 'Holliday-junction resolution',
          'plain': 'A molecular motor that moves the point where two DNA strands cross over.',
          'why': 'Drives the late step of recombination that resolves crossed-over DNA into two '
                 'clean copies.',
          'tags': ['homologous recombination', 'Holliday junction'],
          'role': 'RuvB — Holliday-junction branch migration.',
          'expr': None},
 'mutS': {'scope': 'cons',
          'pathway': 'Mismatch repair',
          'plain': 'A proofreader that spots mismatched or mispaired DNA bases.',
          'why': 'Catches errors so repair and recombination stay accurate.',
          'tags': ['mismatch repair', 'damage recognition'],
          'role': 'MutS — mismatch recognition.',
          'expr': None},
 'gyrA': {'scope': 'cons',
          'pathway': 'DNA topology',
          'plain': 'Part of DNA gyrase, the machine that manages the twist and tension of the '
                   'chromosome.',
          'why': 'Keeps DNA topology in check during the heavy unwinding and resynthesis of '
                 'repair.',
          'tags': ['DNA topology', 'supercoiling'],
          'role': 'GyrA — DNA gyrase subunit A.',
          'expr': None},
 'polA': {'scope': 'cons',
          'pathway': 'Gap filling / excision repair',
          'plain': 'A DNA-copying enzyme that fills in gaps and replaces removed stretches.',
          'why': 'After damaged DNA is cut out, it resynthesizes the missing piece.',
          'tags': ['excision repair', 'DNA synthesis'],
          'role': 'PolA — DNA polymerase I.',
          'expr': None},
 'lexA': {'scope': 'cons',
          'pathway': 'SOS regulation',
          'plain': "A switch that keeps emergency repair genes off until they're needed.",
          'why': "Represses the bacterial SOS response; when damage hits it's cleaved and those "
                 'genes switch on. (D. radiodurans leans more on the IrrE/DdrO switch.)',
          'tags': ['regulation', 'SOS response'],
          'role': 'LexA — SOS-response repressor.',
          'expr': None},
 'tkt': {'scope': 'cons',
         'pathway': 'Metabolic recovery',
         'plain': 'A metabolic enzyme in the pentose-phosphate pathway.',
         'why': "Regenerates NADPH — the cell's antioxidant and building-block fuel — powering the "
                'energy-intensive recovery after irradiation.',
         'tags': ['metabolism', 'metabolic recovery'],
         'role': 'Transketolase — PPP, NADPH regeneration.',
         'expr': None}}
