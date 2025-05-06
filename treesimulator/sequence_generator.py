import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns


def get_normalised_generator(frequencies, rate_matrix=None):
    """
    Calculates the normalised generator from the rate matrix and character state frequencies.

    :param frequencies: character state frequencies.
    :type frequencies: numpy.array
    :param rate_matrix: (optional) rate matrix (by default an all-equal-rate matrix is used)
    :type rate_matrix: numpy.ndarray
    :return: normalised generator 1/mu Q
    :rtype: numpy.ndarray
    """
    if rate_matrix is None:
        n = len(frequencies)
        rate_matrix = np.ones(shape=(n, n), dtype=float) - np.eye(n)
    generator = rate_matrix * frequencies
    generator -= np.diag(generator.sum(axis=1))
    mu = -generator.diagonal().dot(frequencies)
    generator /= mu
    return generator


def simulate_sequence(nucleotide_array, rate_matrix, max_time):
    """
    Simulates the tree evolution from a root over the given time based on the given model.

    :param max_time: float, time over which we generate a tree.
    :return: the simulated tree (ete3.Tree).
    """
    # evolve till the time is up, following Gillespie
    time = 0

    nt_rate_sums = rate_matrix.sum(axis=0)
    event_times = []
    nucleotide_arrays = [nucleotide_array]

    while time < max_time:
        # first we need to calculate rate sum
        rate_sums = nt_rate_sums[nucleotide_arrays[-1]]
        total_rate = rate_sums.sum()

        # now let us see when next event takes place
        time += np.random.exponential(1 / total_rate, 1)[0]

        # Check if the time is up
        if time > max_time:
            break
        event_times.append(time)

        # now let us see which event will happen
        random_event = np.random.uniform(0, 1, 1)[0] * total_rate

        next_nucleotide_array = np.copy(nucleotide_arrays[-1])

        for i in range(len(next_nucleotide_array)):
            # check if this nucleotide is to mutate
            if random_event < rate_sums[i]:
                transition_rates = rate_matrix[nucleotide_arrays[-1][i], :]
                for j in range(4):
                    if random_event < transition_rates[j]:
                        next_nucleotide_array[i] = j
                        # print('Nucleotide at position {i} i changed from {k} to {j}'.format(i=i, k=nucleotide_arrays[-1][i], j=j))
                        break
                    random_event -= transition_rates[j]
                break
            random_event -= rate_sums[i]
            if i == len(next_nucleotide_array) - 1:
                print('No event')
        nucleotide_arrays.append(next_nucleotide_array)

    return nucleotide_arrays, event_times


if __name__ == "__main__":
    rate_matrix = np.ones(shape=(4, 4), dtype=float) - np.eye(4)
    rate_matrix *= np.array([0.7, 0.1, 0.1, 0.1])

    # rate_matrix = get_normalised_generator(np.array([0.3, 0.2, 0.2, 0.3]))
    # rate_matrix -= rate_matrix * np.eye(len(rate_matrix))
    print(rate_matrix)
    N = 100
    nucleotide_array = np.array([0, 1, 2, 3] * N)
    M = len(nucleotide_array)
    nucleotide_arrays, event_times = simulate_sequence(nucleotide_array, rate_matrix, 10)
    print(event_times)
    nucleotide_matrix = np.array(nucleotide_arrays).transpose()
    n2color = {0: 'red', 1: 'green', 2: 'yellow', 3: 'blue'}
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1)

    def plot(ax):
        for i in range(M):
            prev_time = 0
            prev_nuc = nucleotide_matrix[i, 0]
            for time, nuc in zip(event_times, nucleotide_matrix[i, 1:]):
                if time == event_times[-1] or prev_nuc != nuc:
                    ax.plot([prev_time, time], [M - i, M - i], color=n2color[prev_nuc])
                    prev_time = time
                    prev_nuc = nuc
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.get_yaxis().set_visible(False)
        ax.get_yaxis().set_ticks([])

    nucleotide_matrix = nucleotide_matrix[nucleotide_matrix[:, 0].argsort()]
    plot(ax1)

    nucleotide_matrix = nucleotide_matrix[nucleotide_matrix[:, -len(event_times) // 4].argsort()]
    plot(ax2)

    nucleotide_matrix = nucleotide_matrix[nucleotide_matrix[:, -1].argsort()]
    plot(ax3)

    # sns.set(style="whitegrid")
    # sns.heatmap(nucleotide_matrix)
    plt.show()


