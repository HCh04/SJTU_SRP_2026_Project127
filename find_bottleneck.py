"""
find_bottleneck.py

Runs FVA on the product reaction at a fixed growth rate, using MEWpy's
enzyme-constrained model.



from build_glcnac_model import build_glcnac_model, BIOMASS_RXN, EXPORT_RXN_ID


def find_bottleneck(sim, growth_rate, product_rxn=EXPORT_RXN_ID):
    sim.set_reaction_bounds(BIOMASS_RXN, growth_rate, 10000)
    fva = sim.FVA(reactions=[product_rxn], obj_frac=0.9)
    return fva


if __name__ == '__main__':
    sim = build_glcnac_model()
    result = find_bottleneck(sim, growth_rate=0.1)
    print(result)
"""

"""
find_bottleneck.py

Identifies the current rate-limiting enzyme(s) by:
1. Running FVA on the product reaction (original check)
2. Maximizing product flux and ranking candidate proteins by how much
   of their "draw_prot_" usage reaction flux is consumed -- the
   proteins with the highest usage are the current bottlenecks.
"""

from build_glcnac_model import build_glcnac_model, BIOMASS_RXN, EXPORT_RXN_ID

# Same 10 real, model-present proteins used in suggest_targets.py
CANDIDATE_PROTEINS = {
    'HXK1': 'P04806', 'HXK2': 'P04807', 'PGI1': 'P12709', 'UGP1': 'P32861',
    'GFA1': 'P14742', 'GNA1': 'P43577', 'AGM1': 'P38628', 'UAP1': 'P43123',
    'CHS3': 'P29465', 'CHS2': 'P14180',
}


def find_bottleneck(sim, growth_rate, product_rxn=EXPORT_RXN_ID):
    """Return the product flux range with biomass fixed at ``growth_rate``."""
    sim.set_reaction_bounds(BIOMASS_RXN, growth_rate, growth_rate)
    fva = sim.FVA(reactions=[product_rxn], obj_frac=0.0)
    return fva


def rank_protein_usage(sim, top_n=3):
    """
    Maximizes product flux, then ranks candidate proteins by their
    draw_prot_ usage flux. Returns the top_n most-used proteins --
    these are the current bottleneck 'neighborhood'.
    """
    sim.objective = EXPORT_RXN_ID
    result = sim.simulate()

    usage = {}
    for gene, uniprot in CANDIDATE_PROTEINS.items():
        rxn_id = f'draw_prot_{uniprot}'
        flux = abs(result.fluxes.get(rxn_id, 0.0))
        usage[uniprot] = {'gene': gene, 'flux': flux}

    ranked = sorted(usage.items(), key=lambda kv: kv[1]['flux'], reverse=True)
    top = ranked[:top_n]

    print(f"Top {top_n} bottleneck candidates by enzyme usage:")
    for uniprot, info in top:
        print(f"  {info['gene']} ({uniprot}): usage flux = {info['flux']:.6g}")

    return [uniprot for uniprot, _ in top]


if __name__ == '__main__':
    sim = build_glcnac_model()
    result = find_bottleneck(sim, growth_rate=0.1)
    print(result)
    print()
    bottleneck_proteins = rank_protein_usage(sim, top_n=3)
    print("Bottleneck protein IDs:", bottleneck_proteins)
