"""
A D3-based network layout engine for plotting model solutions in the jupyter
window
"""

import os
import json
import itertools
import re

from jinja2 import Environment, FileSystemLoader
from IPython.display import HTML
from csscompressor import compress

from cobra.io.json import model_to_dict
from cobra.exceptions import OptimizationError

import d3flux

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

    default_flux_width:
        If reaction fluxes are missing, the default thickness to use for
        connecting arrows.

    flux_dict:
        A dictionary-like object containing the desired fluxes for each
        reaction in the model

    metabolite_dict:
        A dictionary-like object containing the desired carried fluxes for each
        metabolite in the model

    flowLayout:
        whether or not to use the webcola's flow layout to force a heirarchical
        diagram

    """

    # Initialize empty map_info field in object notes
    for obj in itertools.chain([cobra_model], cobra_model.metabolites,
                               cobra_model.reactions):
        if 'map_info' not in obj.notes:
            obj.notes['map_info'] = {}

    # build cofactor metabolites from strings
    cobra_metabolites = []
    if excluded_metabolites:
        compartments = set((m.compartment for m in cobra_model.metabolites))
        metabolite_list = [
            cf + '_' + co for cf, co in itertools.product(
                excluded_metabolites, compartments)] + excluded_metabolites
        for cofactor in metabolite_list:
            try:
                cobra_metabolites += [
                    cobra_model.metabolites.get_by_id(cofactor)]
            except KeyError:
                pass
            # what if its already a cobra metabolite?

    cobra_rxns = []
    if excluded_reactions:
        for rxnid in excluded_reactions:
            try:
                cobra_rxns += [
                    cobra_model.reactions.get_by_id(rxnid)]
            except KeyError:
                pass

    # Exclude metabolites and reactions in the given comparment
    excluded_metabolites = set(cobra_metabolites)
    excluded_reactions = set(cobra_rxns)

    if excluded_compartments:
        met_compartments = set((m for m in cobra_model.metabolites.query(
            lambda x: set(x.compartment).intersection(
                set(excluded_compartments)), None)))
        excluded_metabolites |= met_compartments

        # Do I want to redo this not to include excluded metabolites?
        # rxn_compartments = set((r for r in cobra_model.reactions.query(
        #     lambda x: set(x.compartments).intersection(
        #         set(excluded_compartments)), None)))
        # excluded_reactions |= rxn_compartments

    # for reaction in excluded_reactions:
    #     reaction.notes['map_info'] = {'hidden': True}

    for metabolite in excluded_metabolites:
        metabolite.notes['map_info'] = {'hidden': True}

    def is_hidden(obj):
        try:
            return bool(obj.notes['map_info']['hidden'])
        except KeyError:
            return False

    for reaction in cobra_model.reactions:

        if overwrite_reversibility:
            reaction.notes['map_info']['reversibility'] = \
                bool(reaction.reversibility)

        # Unless 'hidden' specifically set to False, hide the reaction if all
        # the reactants or products are hidden (excluding cofactors)
        if ('hidden' not in reaction.notes['map_info']) & ~is_hidden(reaction):

            # Hide reactions if all of their products or reactants are hidden.
            # Don't include cofactor metabolites in this calculation.
            if 'cofactors' in reaction.notes['map_info']:
                cofactors = reaction.notes['map_info']['cofactors'].keys()
            else:
                cofactors = {}

            if (all([is_hidden(met) for met in reaction.reactants
                     if met.id not in cofactors]) or
                all([is_hidden(met) for met in reaction.products
                     if met.id not in cofactors])):

                reaction.notes['map_info']['hidden'] = True

    # Add diplay names to the cobra metabolites accoring to the
    # display_name_format function
    if display_name_format:

        # Handle the case for a default display name formatter. This is
        # optimized for models using the typical bigg_id naming convention,
        # ending with 'ID_c' compartment identifier.
        if display_name_format is True:
            display_name_format = (
                lambda met: re.sub('__[D,L]', '', met.id[:-2].upper()))

        for met in cobra_model.metabolites:

            # Don't overwrite existing display names
            if 'display_name' not in met.notes['map_info']:
                met.notes['map_info']['display_name'] = (
                    display_name_format(met))

    # Append model's map_info kwargs
    render_kwargs = dict(cobra_model.notes['map_info'])
    render_kwargs.update(kwargs)

    return render_model(cobra_model, **render_kwargs)


def create_model_json(cobra_model, flux_dict=None, metabolite_dict=None):
    """ Convert a cobra.Model object to a json string for d3. Adds flux
    information if the model has been solved

    flux_dict: dict-like
        Contains an external setting of the flux solution that should be
        plotted.

    metabolite_dict:
        A dictionary-like object containing the desired carried fluxes for each
        metabolite in the model

    """
    def get_flux(reaction):
        if flux_dict is not None:
            return flux_dict[reaction.id]
        else:

            try:
                return reaction.flux

            except OptimizationError:
                # The model hasn't been solved, so we just throw in a None
                return None


    # Add flux info
    for reaction in cobra_model.reactions:

        # If I'm styling reaction knockouts, don't set the flux for a
        # knocked out reaction
        if reaction.lower_bound == reaction.upper_bound == 0:
            reaction.notes['map_info']['group'] = 'ko'
            
            # Delete the flux key, if it exists
            try:
                del reaction.notes['map_info']['flux']
            except KeyError:
                pass

        else: 
            try:
                if abs(get_flux(reaction)) < 1E-8:
                    reaction.notes['map_info']['flux'] = 0.
                else:
                    reaction.notes['map_info']['flux'] = get_flux(reaction)
            except (KeyError, TypeError):
                if 'flux' in reaction.notes['map_info']:
                    del reaction.notes['map_info']['flux']

            # cobrapy doesn't track contexted changes to the notes field. So if
            # a reaction is set to the 'ko' group, reset it if it doens't match
            # the bounds requirements
            if 'group' in reaction.notes['map_info']:
                if reaction.notes['map_info']['group'] == 'ko':
                    del reaction.notes['map_info']['group']

    def get_met_flux(metabolite):
        if metabolite_dict is not None:
            return metabolite_dict[metabolite.id]
        else:
            return sum([abs(get_flux(r) * r.metabolites[metabolite]) for r in
                        metabolite.reactions]) / 2

    for metabolite in cobra_model.metabolites:

        try:
            del metabolite.notes['map_info']['flux']

        except KeyError:
            pass

        try:
            carried_flux = get_met_flux(metabolite)
            if carried_flux > 1E-8:
                metabolite.notes['map_info']['flux'] = carried_flux
            else:
                metabolite.notes['map_info']['flux'] = 0.

        except Exception:
            pass

    return json.dumps(model_to_dict(cobra_model), allow_nan=False)


def render_model(cobra_model, background_template=None, custom_css=None,
                 figure_id=None, hide_unused=None, hide_unused_cofactors=None,
                 inactive_alpha=1., figsize=None, label=None, fontsize=None,
                 default_flux_width=2.5, flux_dict=None, metabolite_dict=None,
                 svg_scale=100, flowLayout=False):
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

    inactive_alpha:
        Alpha value with which to color reactions and nodes without any carried
        flux. Defaults to 1.

    figsize:
        size, in pixels, of the generated SVG window. Defaults to 1024x768.

    fontsize:
        text size, in pt. Defaults to 12

    default_flux_width:
        If reaction fluxes are missing, the default thickness to use for
        connecting arrows.

    flux_dict:
        A dictionary-like object containing the desired fluxes for each
        reaction in the model

    metabolite_dict:
        A dictionary-like object containing the desired carried fluxes for each
        metabolite in the model

    """

    # Increment figure counter

    # Get figure name and JSON string for the cobra model
    if not figure_id:
        render_model._fignum += 1
        figure_id = 'd3flux{:0>3d}'.format(render_model._fignum)

    if not figsize:
        figsize = (1028, 768)

    modeljson = create_model_json(cobra_model, flux_dict, metabolite_dict)

    if not hide_unused:
        hide_unused = "false"
    else:
        hide_unused = "true"

    if not hide_unused_cofactors:
        hide_unused_cofactors = "false"
    else:
        hide_unused_cofactors = "true"

    # Handle custom CSS
    if not custom_css:
        custom_css = ''

    if not fontsize:
        fontsize = 12

    # Handle background template
    if not background_template:
        background_svg = ''
        no_background = "true"
    else:
        from IPython.display import SVG
        background_svg = SVG(background_template).data
        no_background = "false"

    # Initialize the jinja templates
    env = Environment(loader=FileSystemLoader(
        os.path.join(os.path.dirname(d3flux.__file__), 'templates')))
    
    template_css = env.get_template('network_style.css')
    template_html = env.get_template('output_template.html')
    template_js = env.get_template('d3flux.js')

    # Render the jinja templates with the given variables
    css = template_css.render(inactive_alpha=inactive_alpha,
                              fontsize=fontsize, cf_fontsize=0.8 * fontsize)

    js = template_js.render(figure_id=figure_id, modeljson=modeljson,
                            no_background=no_background,
                            hide_unused=hide_unused,
                            hide_unused_cofactors=hide_unused_cofactors,
                            figwidth=figsize[0], figheight=figsize[1],
                            css=compress(css + custom_css),
                            default_flux_width=default_flux_width,
                            svg_scale=svg_scale, flowLayout=flowLayout)

    html = template_html.render(figure_id=figure_id,
                                background_svg=background_svg,
                                javascript_source=js)

    # compile and return HTML
    return HTML(html)


# Initialize figure counter
render_model._fignum = 0
