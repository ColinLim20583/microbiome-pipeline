import os
from bioblend.galaxy import GalaxyInstance

# Never commit credentials. Read them from the environment instead:
#   export GALAXY_URL=http://localhost:8080
#   export GALAXY_API_KEY=your-key-here
galaxy_url = os.environ.get("GALAXY_URL", "http://localhost:8080")
api_key = os.environ.get("GALAXY_API_KEY")
if not api_key:
    raise SystemExit(
        "GALAXY_API_KEY environment variable is not set. "
        "Set it before running: export GALAXY_API_KEY=your-key"
    )

gi = GalaxyInstance(url=galaxy_url, key=api_key)

# Upload a file to Galaxy history
history = gi.histories.create_history(name='PICRUSt2 run')
fasta = gi.tools.upload_file('/path/to/your_seqs.fasta', history['id'])
table = gi.tools.upload_file('/path/to/your_table.biom', history['id'])

# Find and run the PICRUSt2 tool
tools = gi.tools.get_tools(name='picrust2_pipeline')
picrust2_tool_id = tools[0]['id']

gi.tools.run_tool(history_id=history['id'], tool_id=picrust2_tool_id, inputs={
    'input_sequences': {'src': 'hda', 'id': fasta['outputs'][0]['id']},
    'input_table': {'src': 'hda', 'id': table['outputs'][0]['id']},
    'threads': '2',
    'stratified': True
})
