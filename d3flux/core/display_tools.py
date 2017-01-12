import pandas as pd
from functools import reduce
# TODO: write docstrings

def metabolite_summary(met):
    return pd.Series({r.id : (r.x * r.metabolites[met]) 
                      for r in met.reactions}, name='flux')

common_ox_cofactors = ['nad_c', 'nadp_c', 'q8_c']

def redox_summary(cobra_model, tol=1E-8, ox_cofactors=None):

    if ox_cofactors == None: ox_cofactors = common_ox_cofactors
    elif not ox_cofactors: return pd.Series(), pd.Series()

    redox_series = reduce(lambda x, y: x.add(y, fill_value=0), (
            metabolite_summary(cobra_model.metabolites.get_by_id(cofactor))
            for cofactor in ox_cofactors))
    oxidizing = redox_series[redox_series > tol]
    reducing = redox_series[redox_series < -tol]
    return oxidizing, reducing

def color_redox_rxns(cobra_model, reset_groups=True, color_knockouts=True,
                     starting_group=1, **kwargs):
    """Add group info to the cobra_model according to the results of
    `redox_summary`. 

    cobra_model: a cobra.Model object

    reset_groups: bool
        whether or not to delete group info from model reactions prior to
        assigning new groups.

    color_knockouts: bool
        whether or not to conditionally format reactions which are currently
        knocked out. (i.e., rxn.bounds = (0,0))

    starting_group: int
        To use different colors, start from a group other than 1. (Highest
        color is 8)

    Additional kwargs are passed directly to the `redox_summary` function call:

    tol: int, default 1E-8
        flux tolerance above which to include the reaction
    
    ox_cofactors: list, default is common_ox_cofactors
        Cofactors which are produced or consumed in several reactions to use as
        the basis for coloring the flux model.

    """
    # Pop existing group info
    if reset_groups:
        for rxn in cobra_model.reactions:
            try: del rxn.notes['map_info']['group']
            except KeyError: pass

    # Assign group for knocked out reactions
    if color_knockouts:
        for rxn in cobra_model.reactions.query(
                lambda r: r.bounds == (0, 0), None):
            rxn.notes['map_info']['group'] = 'ko'

    # Assign new group info based on redox balance
    oxidizing, reducing = redox_summary(cobra_model, **kwargs)
    for rxn, flux in oxidizing.items():
        cobra_model.reactions.get_by_id(
            rxn).notes['map_info']['group'] = starting_group
    for rxn, flux in reducing.items():
        cobra_model.reactions.get_by_id(
            rxn).notes['map_info']['group'] = starting_group + 1

    return cobra_model

                                            
def update_cofactors(cobra_model, cofactor_list):
    """ Given a model and list of cofactors to display, update the map_info
    field of each reaction containing the given cofactors. Updates the model
    inplace.
    
    cobra_model : a cobra.Model object
    cofactor_list : a list of strings indicating the desired cofactor IDs to show.

    """

    def rxn_cofactor_update(rxn, cofactor):
    
        if 'map_info' not in rxn.notes:
            rxn.notes['map_info'] = {}

        try:
            if cofactor in rxn.notes['map_info']['cofactors']: return
            rxn.notes['map_info']['cofactors'].update({cofactor: {}})
        except KeyError:
            rxn.notes['map_info']['cofactors'] = {cofactor: {}}

    for met_id in cofactor_list:
        this_met = cobra_model.metabolites.get_by_id(met_id)
        for rxn in this_met.reactions: rxn_cofactor_update(rxn, met_id)

common_cofactors = ['coa', 'nadh', 'nad', 'nadph', 'nadp', 'atp', 'adp', 'amp',
                    'q8', 'q8h2', 'pi', 'co2', 'h2o', 'h', 'o2', 'h2', 'nh4']

