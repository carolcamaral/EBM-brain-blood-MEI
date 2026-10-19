import re, sys, nbformat as nbf
src = open('nb_source.py').read()
parts = re.split(r'^# %%(.*)$', src, flags=re.M)
cells = []
for i in range(1, len(parts), 2):
    kind, body = parts[i].strip(), parts[i+1].strip('\n')
    if kind == '[markdown]':
        cells.append(nbf.v4.new_markdown_cell('\n'.join(l[2:] if l.startswith('# ') else l.lstrip('#') for l in body.splitlines())))
    else:
        cells.append(nbf.v4.new_code_cell(body))
nb = nbf.v4.new_notebook(); nb.cells = cells
nb.metadata['kernelspec'] = {'name': 'python3', 'display_name': 'Python 3', 'language': 'python'}
out = sys.argv[1]
nbf.write(nb, out)
from nbclient import NotebookClient
NotebookClient(nb, timeout=600, kernel_name='python3', resources={'metadata': {'path': 'data'}}).execute()
nbf.write(nb, out)
print('executed ok')
