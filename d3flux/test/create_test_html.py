from cobra.io import load_json_model
from d3flux import flux_map

new_model = load_json_model('simple_model_no_layout.json')

new_model.objective = new_model.reactions.R4
# new_model.reactions.R8.knock_out()
new_model.reactions.R6.knock_out()
new_model.optimize()

new_model.metabolites.B.notes['map_info']['align'] = 'left'
new_model.metabolites.P.notes['map_info']['align'] = 'lower center'
new_model.metabolites.C.notes['map_info']['align'] = 'upper right'

new_model.reactions.R5.bounds = (-1000, 1000)

html = flux_map(new_model, figsize=(900,250), inactive_alpha=0.5)
with open('test.html', 'w') as f:
    f.write('<!DOCTYPE html> <html> <head> <title>Test d3flux page</title> <script src="https://cdnjs.cloudflare.com/ajax/libs/require.js/2.3.2/require.min.js" integrity="sha256-Vjusm6Kh2U7/tb6jBh+MOfxnaf2TWsTph34bMKhC1Qc=" crossorigin="anonymous"></script> <script src="https://code.jquery.com/jquery-3.1.1.min.js" integrity="sha256-hVVnYaiADRTO2PzUGmuLJr8BLUSjGIZsDYGmIJLv2b8=" crossorigin="anonymous"></script> <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css" integrity="sha384-BVYiiSIFeK1dGmJRAkycuHAHRg32OmUcww7on3RYdg4Va+PmSTsz/K68vbdEjh4u" crossorigin="anonymous"> </head> <body>')
    f.write(html.data)
    f.write('</body> </html>')
