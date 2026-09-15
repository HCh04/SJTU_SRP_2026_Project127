"""
suggest_targets.py

Enzyme-aware strain design using MEWpy's GeckoOUProblem — searches over
enzyme over/under-expression (not just knockouts) to maximize GlcNAc
production while keeping the cell viable. This replaces plain OptKnock,
since we now have real enzyme constraints to search over.
"""

"""
from mewpy.problems.gecko import GeckoOUProblem
from mewpy.optimization.evaluation import BPCY, WYIELD
from mewpy.optimization import EA

from build_glcnac_model import build_glcnac_model, BIOMASS_RXN, EXPORT_RXN_ID


def suggest_targets(sim, max_generations=2):
    f1 = BPCY(BIOMASS_RXN, EXPORT_RXN_ID, method='lMOMA')
    f2 = WYIELD(BIOMASS_RXN, EXPORT_RXN_ID)
    problem = GeckoOUProblem(sim.model, fevaluation=[f1, f2])
    ea = EA(problem, max_generations=max_generations)
    ea.run()
    return ea


if __name__ == '__main__':
    sim = build_glcnac_model()
    ea = suggest_targets(sim, max_generations=2)  # small run for a quick local test
    print(ea)
"""

from build_glcnac_model import BIOMASS_RXN, EXPORT_RXN_ID, build_glcnac_model
from mewpy.optimization import EA
from mewpy.optimization.evaluation import BPCY, WYIELD
from mewpy.problems.gecko import GeckoOUProblem

def suggest_targets(
    sim,
    max_generations=3,
    pop_size=20,
    use_multiprocessing=False,
):
    # 1. Fast evaluation using FBA instead of heavy lMOMA
    f1 = BPCY(BIOMASS_RXN, EXPORT_RXN_ID, method='FBA')
    f2 = WYIELD(BIOMASS_RXN, EXPORT_RXN_ID)

    # 2. Restrict search space -- must use real UniProt IDs, not gene symbols
    target_enzymes = [
        'P04806', 'P04807',  # HXK1, HXK2
        'P12709',            # PGI1
        'P32861',            # UGP1
        'P14742',            # GFA1
        'P43577',            # GNA1
        'P38628',            # AGM1
        'P43123',            # UAP1
        'P29465', 'P14180',  # CHS3, CHS2
        'P23797',            # GPI12
        # GPI3/GPI12/ALG7/DOG1 are NOT in this bundled model's proteome
        # (confirmed earlier via get_Kcats) -- omitted, would be ignored anyway
    ]

    valid_targets = [e for e in target_enzymes if e in sim.model.proteins]
    if not valid_targets:
        valid_targets = None
        print("WARNING: no valid targets matched -- falling back to full search space")
    else:
        print(f"Restricting search to {len(valid_targets)} proteins: {valid_targets}")

    problem = GeckoOUProblem(
        sim.model,
        fevaluation=[f1, f2],
        target=valid_targets,   # correct kwarg name
    )

    ea = EA(
        problem,
        max_generations=max_generations,
        population_size=pop_size,
        mp=use_multiprocessing,
    )

    ea.run()
    return ea


if __name__ == '__main__':
    sim = build_glcnac_model()
    ea = suggest_targets(
        sim,
        max_generations=3,
        pop_size=20,
        use_multiprocessing=False,
    )
    results = ea.dataframe()
    print(results.to_string())
    results.to_csv('suggest_targets_results.csv', index=False)
