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


def save_ltt(real_ltt, observed_ltt, ltt_file):
    os.makedirs(os.path.dirname(os.path.abspath(ltt_file)), exist_ok=True)
    with open(ltt_file, 'w+') as f:
        f.write('time\treal lineages\tobserved lineages\n')
        observed = 0
        real = 0
        for time in sorted(set(real_ltt.keys()) | set(observed_ltt.keys())):
            if time in observed_ltt:
                observed = observed_ltt[time]
            if time in real_ltt:
                real = real_ltt[time]
            f.write('{}\t{}\t{}\n'.format(time, real, observed))


def save_log(model, total_tips, T, u, log):
    res = model.get_epidemiological_parameters()
    res['tips'] = total_tips
    res['time'] = T
    res['hidden_trees'] = u
    os.makedirs(os.path.dirname(os.path.abspath(log)), exist_ok=True)
    with open(log, 'w+') as f:
        f.write('{}\n'.format(','.join(res.keys())))
        f.write('{}\n'.format(','.join(str(_) for _ in res.values())))
