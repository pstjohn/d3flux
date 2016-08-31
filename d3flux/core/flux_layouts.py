"""
A D3-based network layout engine for plotting model solutions in the jupyter window
"""

import os
import json
import itertools
import re

from jinja2 import Environment, FileSystemLoader
from IPython.display import HTML

from ..io.json import _to_dict

env = Environment(loader=FileSystemLoader(
    os.path.join(os.path.dirname(__file__), 'templates')))


def flux_map(cobra_model, 
             excluded_metabolites=None, excluded_reactions=None,
             excluded_compartments=None, display_name_format=True,
             overwrite_reversibility=True, **kwargs):
    """Create a flux map representation of the cobra.Model, including or
    excluding the given metabolites, reactions, and/or compartments.

    excluded_metabolites:
        A list of metabolites to not include. Probably better to simply set
        met.notes['map_info']['hidden'] = True

    excluded_reactions:
        A list of reactions not to display.

    excluded_compartment:
        Hide metabolites not in a compartment

    diplay_name_format: Bool or function
        How to format the metabolite names in map.display_name.

    Additional kwargs are passed directly to `render_model`:

    background_template: 
        filename for an SVG to render behind the flux figure.  Useful for
        compartments or layout guides.

    custom_css:
        Additional CSS to embed in the figure. Use HTML inspector to show
        labels and classes applied to reactions and metabolites.

    figure_id:
        Each figure in the page requires a unique ID, which can be passed or
        generated automatically.

    hide_unused:
        whether or not to show metabolites and reactions with zero flux.

    hide_unused_cofactors:
        similar to hide_unused, but only hide cofactor nodes for reactions with
        0 flux.

    overwrite_reversibility:
        whether or not to overwrite the `reversibility` attribute in map_info.
        Useful in ensuring reversible, knocked-out reactions are rendered
        appropriately. Defaults to True.

    figsize:
        size, in pixels, of the generated SVG window. Defaults to 1024x768.

    fontsize:
        text size, in pt. Defaults to 12
    
    """
    # build cofactor metabolites from strings
    cobra_metabolites = []
    if excluded_metabolites:
        compartments = cobra_model.get_compartments()
        metabolite_list = [
            cf + '_' + co for cf, co in itertools.product(
                excluded_metabolites, compartments)] + excluded_metabolites
        for cofactor in metabolite_list:
            try: 
                cobra_metabolites += [
                    cobra_model.metabolites.get_by_id(cofactor)]
            except KeyError: pass
            # what if its already a cobra metabolite?

    cobra_rxns = []
    if excluded_reactions:
        for rxnid in excluded_reactions:
            try: 
                cobra_rxns += [
                    cobra_model.reactions.get_by_id(rxnid)]
            except KeyError: pass


    # Exclude metabolites and reactions in the given comparment
    excluded_metabolites = set(cobra_metabolites)
    excluded_reactions = set(cobra_rxns)

    if excluded_compartments:
        met_compartments = set((m for m in cobra_model.metabolites.iquery(
            lambda x: set(x.compartment).intersection(
                set(excluded_compartments)))))
        excluded_metabolites |= met_compartments

        # Do I want to redo this not to include excluded metabolites? 
        rxn_compartments = set((r for r in cobra_model.reactions.iquery(
            lambda x: set(x.compartments).intersection(
                set(excluded_compartments)))))
        excluded_reactions |= rxn_compartments

    for reaction in excluded_reactions:
        reaction.notes['map_info'] = {'hidden' : True}

    for metabolite in excluded_metabolites:
        metabolite.notes['map_info'] = {'hidden' : True}

    def is_hidden(obj):
        try: return bool(obj.notes['map_info']['hidden'])
        except KeyError: return False

    for reaction in cobra_model.reactions:
        
        if overwrite_reversibility:
            try:
                reaction.notes['map_info']['reversibility'] = \
                reaction.reversibility
            except KeyError:
                reaction.notes['map_info'] = {
                    'reversibility': reaction.reversibility}

        # Hide reactions if all of their products or reactants are hidden
        if (all([is_hidden(met) for met in reaction.reactants]) or
            all([is_hidden(met) for met in reaction.products])):

            try:
                reaction.notes['map_info']['hidden'] = True
            except KeyError:
                reaction.notes['map_info'] = {'hidden': True}

    # Add diplay names to the cobra metabolites accoring to the
    # display_name_format function
    if display_name_format:

        # Handle the case for a default display name formatter.
        if display_name_format == True:
            display_name_format = (
                lambda met: re.sub('__[D,L]', '', met.id[:-2].upper()))

        for met in cobra_model.metabolites:

            if 'map_info' not in met.notes:
                met.notes['map_info'] = {}

            # Don't overwrite existing display names
            if 'display_name' not in met.notes['map_info']:
                met.notes['map_info']['display_name'] = (
                    display_name_format(met))
        
    return render_model(cobra_model, **kwargs)


def create_model_json(cobra_model):
    """ Convert a cobra.Model object to a json string for d3. Adds flux
    information if the model has been solved
    
    """
    # Add flux info
    for reaction in cobra_model.reactions:
        try:
            # If I'm styling reaction knockouts, don't set the flux for a
            # knocked out reaction
            try: 
                if reaction.notes['map_info']['group'] == 'ko': 
                    # Delete the flux key, if it exists
                    try:
                        reaction.notes['map_info']['flux'] = 0
                    except Exception: pass

                    # Onto the next reaction
                    continue

            # Onto the next reaction
            except KeyError: pass

            reaction.notes['map_info']['flux'] = reaction.x

        except KeyError:
            # Create a new "map_info" object
            reaction.notes['map_info'] = {'flux' : reaction.x}

        except AttributeError:
            # Model likely hasn't been solved, get out now
            return json.dumps(_to_dict(cobra_model), allow_nan=False)

    for metabolite in cobra_model.metabolites:
        try:
            carried_flux = sum([abs(r.x * r.metabolites[metabolite]) for r in
                                metabolite.reactions])/2
            metabolite.notes['map_info']['flux'] = carried_flux

        except KeyError:
            # Create a new "map_info" object
            metabolite.notes['map_info'] = {'flux' : carried_flux}

        except AttributeError:
            # Model likely hasn't been solved, get out now
            return json.dumps(_to_dict(cobra_model), allow_nan=False)

    return json.dumps(_to_dict(cobra_model), allow_nan=False)



def render_model(cobra_model, background_template=None, custom_css=None,
                 figure_id=None, hide_unused=None, hide_unused_cofactors=None,
                 figsize=None, label=None,
                 fontsize=None):
    """ Render a cobra.Model object in the current window 
    
    Parameters:

    background_template: 
        filename for an SVG to render behind the flux figure.  Useful for
        compartments or layout guides.

    custom_css:
        Additional CSS to embed in the figure. Use HTML inspector to show
        labels and classes applied to reactions and metabolites.

    figure_id:
        Each figure in the page requires a unique ID, which can be passed or
        generated automatically.

    hide_unused:
        whether or not to show metabolites and reactions with zero flux.

    hide_unused_cofactors:
        similar to hide_unused, but only hide cofactor nodes for reactions with
        0 flux.

    figsize:
        size, in pixels, of the generated SVG window. Defaults to 1024x768.

    fontsize:
        text size, in pt. Defaults to 12
        

    """

    # Increment figure counter

    # Get figure name and JSON string for the cobra model
    if not figure_id:
        render_model._fignum += 1
        figure_id = 'd3flux{:0>3d}'.format(render_model._fignum)


    if not figsize:
        figsize = (1028, 768)

    modeljson = create_model_json(cobra_model)

    if not hide_unused:
        hide_unused = "false"
    else: hide_unused = "true"

    if not hide_unused_cofactors:
        hide_unused_cofactors = "false"
    else: hide_unused_cofactors = "true"

    # Handle custom CSS
    if not custom_css: 
        custom_css = ''

    if not fontsize: fontsize = 12

    # Handle background template
    if not background_template:
        background_svg = ''
        no_background = "true"
    else:
        from IPython.display import SVG
        background_svg = SVG(background_template).data
        no_background = "false"

    # Initialize the jinja templates
    template_css = env.get_template('network_style.css')
    template_html = env.get_template('output_template.html')
    template_js = env.get_template('d3flux.js')

    # Render the jinja templates with the given variables
    css = template_css.render(figure_id=figure_id, figwidth=figsize[0],
                              figheight=figsize[1])

    js = template_js.render(figure_id=figure_id, modeljson=modeljson,
                            no_background=no_background,
                            hide_unused=hide_unused,
                            hide_unused_cofactors=hide_unused_cofactors,
                            figwidth=figsize[0], figheight=figsize[1],
                            fontsize=fontsize)

    html = template_html.render(figure_id=figure_id,
                                background_svg=background_svg,
                                default_style=css, custom_css=custom_css,
                                javascript_source=js)

    # compile and return HTML
    return HTML(html)

# Initialize figure counter
render_model._fignum = 0


