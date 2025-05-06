import os

from treesimulator.mtbd_models import CTModel

STATE = 'state'
DIST_TO_START = 'D'
TIME_TILL_NOW = 'T'


def save_forest(forest, nwk, state_feature=STATE):
    os.makedirs(os.path.dirname(os.path.abspath(nwk)), exist_ok=True)
    with open(nwk, 'w+') as f:
        for tree in forest:
            features = [state_feature] if state_feature is not None else []
            nwk = tree.write(format=5, format_root_node=True, features=features)
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


def save_log(models, skyline_times, total_tips, T, u, log, kappa=0, observed_frequencies=None):
    os.makedirs(os.path.dirname(os.path.abspath(log)), exist_ok=True)
    if skyline_times is None:
        skyline_times = []
    skyline_times = [_ for _ in skyline_times if _ <= T]
    skyline_times += [T]
    is_ct = isinstance(models[0], CTModel)
    with open(log, 'w+') as f:
        states = models[0].states
        pis = [f'pi_{s}' for s in states]
        keys = sorted([k for k in models[0].get_epidemiological_parameters().keys() if k not in set(pis)],
                      key=lambda _: (_.split('_')[0], _))

        f.write('{}{},{}tips,hidden_trees,end_time\n'\
                .format(','.join(keys),
                        ',kappa' if is_ct else '',
                        (','.join(f'{pi},{pi}_observed' for pi in pis) + ',') if len(states) > 1 else ''))
        for model, end_time, obs in zip(models, skyline_times, observed_frequencies):
            tips = '' if end_time < T else total_tips
            params = model.get_epidemiological_parameters()
            f.write('{}{},{}{},{},{:g}\n'
                    .format(','.join(f'{params[k]:g}' for k in keys),
                            f',{kappa:g}' if is_ct else '',
                            (','.join(f'{pi:g},{o:g}' for (pi, o) in zip((params[_] for _ in pis), obs)) + ',')
                                                if len(states) > 1 else '',
                            tips, u, end_time))
            u = ''
