import os

STATE = 'state'
DIST_TO_START ='D'
TIME_TILL_NOW = 'T'


def save_forest(forest, nwk):
    os.makedirs(os.path.dirname(os.path.abspath(nwk)), exist_ok=True)
    with open(nwk, 'w+') as f:
        for tree in forest:
            nwk = tree.write(format=5, format_root_node=True)
            f.write('{}\n'.format(nwk))


def save_log(model, total_tips, T, u, log):
    res = model.get_epidemiological_parameters()
    res['tips'] = total_tips
    res['time'] = T
    res['hidden_trees'] = u
    os.makedirs(os.path.dirname(os.path.abspath(log)), exist_ok=True)
    with open(log, 'w+') as f:
        f.write('{}\n'.format(','.join(res.keys())))
        f.write('{}\n'.format(','.join(str(_) for _ in res.values())))
