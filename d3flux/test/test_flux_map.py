import pytest
from cobra.core import Metabolite, Reaction, Model
from d3flux import flux_map

@pytest.fixture()
def simple_model():

    model = Model('simple_model')

    A = Metabolite('A')
    B = Metabolite('B')
    C = Metabolite('C')
    D = Metabolite('D')
    E = Metabolite('E')
    P = Metabolite('P')

    R1 = Reaction('R1')
    R2 = Reaction('R2')
    R3 = Reaction('R3')
    R4 = Reaction('R4')
    R5 = Reaction('R5')
    R6 = Reaction('R6')
    R7 = Reaction('R7')
    R8 = Reaction('R8')
    R9 = Reaction('R9')
    R10 = Reaction('R10')

    model.add_metabolites([A, B, C, D, E, P])
    model.add_reactions([R1, R2, R3, R4, R5, R6, R7, R8, R9, R10])

    model.reactions.R1.build_reaction_from_string('--> A')
    model.reactions.R2.build_reaction_from_string('<--> B')
    model.reactions.R3.build_reaction_from_string('P -->')
    model.reactions.R4.build_reaction_from_string('E -->')
    model.reactions.R5.build_reaction_from_string('A --> B')
    model.reactions.R6.build_reaction_from_string('A --> C')
    model.reactions.R7.build_reaction_from_string('A --> D')
    model.reactions.R8.build_reaction_from_string('B <--> C')
    model.reactions.R9.build_reaction_from_string('B --> P')
    model.reactions.R10.build_reaction_from_string('C + D --> E + P')

    return model


def test_flux_map(simple_model):
    svg = flux_map(simple_model, display_name_format=lambda x: str(x.id),
                   figsize=(300,250), flux_dict={rxn.id: None for rxn in
                                                 simple_model.reactions})

    assert svg is not None


def test_flux_map_optimize(simple_model):
    simple_model.objective = simple_model.reactions.R4
    simple_model.optimize()
    svg = flux_map(simple_model, figsize=(300,250))

    assert svg is not None

def test_flux_map_knockout(simple_model):
    simple_model.objective = simple_model.reactions.R4
    simple_model.reactions.R8.knock_out()
    simple_model.optimize()
    svg = flux_map(simple_model, figsize=(300,250))

    assert svg is not None

