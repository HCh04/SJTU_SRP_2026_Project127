"""
build_glcnac_model.py

Loads MEWpy's bundled enzyme-constrained yeast model (GECKO-formatted,
via reframed), adds the two GlcNAc-specific reactions, and applies real
kcat values from Student A's data — per request for an
enzyme/kcat-constrained model.

Per Student A: the dephosphorylation step (UDP-GlcNAc -> free GlcNAc) is
intentionally gene-less and unconstrained —native yeast has no dedicated
enzyme for it (routes UDP-GlcNAc into chitin/glycans instead). This is a
confirmed modeling decision, not a placeholder awaiting an enzyme.
"""

import csv
from pathlib import Path

from mewpy.model.gecko import GeckoModel
from mewpy.simulation import get_simulator

UDPGLCNAC_MET = 's_1544_c'
UDP_MET = 's_1538_c'
GLCNAC_MET = 's_4017_c'
DEPHOSPHO_RXN_ID = 'r_GLCNAC_1'
EXPORT_RXN_ID = 'r_GLCNAC_EXCHANGE'
BIOMASS_RXN = 'r_2111'
MODULE_DIR = Path(__file__).resolve().parent
KCAT_FILE = MODULE_DIR / 'data' / 'customKcats_GlcNAc.tsv'


def load_kcat_data(path=KCAT_FILE):
    """Reads the custom kcat file, returns a list of dicts.
    Rows with structure_status == 'REJECT' are excluded by default."""
    rows = []
    with open(path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            if row['structure_status'] == 'REJECT':
                print(f"  [SKIPPED - low confidence] {row['gene']} ({row['uniprot']}): {row['notes']}")
                continue
            rows.append(row)
    return rows


def build_glcnac_model(kcat_path=KCAT_FILE):
    """
    Returns
    -------
    mewpy.simulation.reframed.GeckoSimulation
        Enzyme-constrained model with GlcNAc reactions + real kcats applied.
    """
    model = GeckoModel('single-pool')
    sim = get_simulator(model)

    sim.add_metabolite(
        GLCNAC_MET, formula='C8H15NO6',
        name='N-Acetyl-D-glucosamine [cytoplasm]', compartment='c',
    )

    sim.add_reaction(
        rxn_id=DEPHOSPHO_RXN_ID,
        name='UDP-GlcNAc dephosphorylation (unconstrained)',
        stoichiometry={UDPGLCNAC_MET: -1, GLCNAC_MET: 1, UDP_MET: 1},
        reversible=False, lb=0, ub=1000, gpr=None,
    )

    sim.add_reaction(
        rxn_id=EXPORT_RXN_ID,
        name='GlcNAc exchange',
        stoichiometry={GLCNAC_MET: -1},
        reversible=True, lb=-1000, ub=1000, gpr=None,
    )

    print("Applying real kcat values:")
    kcat_rows = load_kcat_data(kcat_path)
    applied = []
    not_found = []
    for row in kcat_rows:
        uniprot = row['uniprot']
        new_kcat = float(row['kcat_s1'])
        try:
            existing = sim.get_Kcats(uniprot)
        except ValueError:
            not_found.append((row['gene'], uniprot))
            continue
        for rxn_id in existing:
            if rxn_id.startswith('draw_prot_'):
                continue
            sim.set_Kcat(uniprot, rxn_id, new_kcat)
            applied.append((row['gene'], uniprot, rxn_id, new_kcat))

    for gene, uniprot, rxn_id, kcat in applied:
        print(f"  [APPLIED] {gene} ({uniprot}) -> {rxn_id}: kcat = {kcat}")
    for gene, uniprot in not_found:
        print(f"  [NOT FOUND IN MODEL] {gene} ({uniprot})")

    return sim


if __name__ == '__main__':
    sim = build_glcnac_model()
    print(f"\n{DEPHOSPHO_RXN_ID} present:", DEPHOSPHO_RXN_ID in sim.reactions)
    print(f"{EXPORT_RXN_ID} present:", EXPORT_RXN_ID in sim.reactions)

    sim.objective = BIOMASS_RXN
    result = sim.simulate()
    print(f"Growth rate: {result.objective_value:.5f} /hour")
