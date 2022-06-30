import pickle

# Plotting
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.ticker import PercentFormatter


# Load data (deserialize)
with open('tanim_freq.pickle', 'rb') as handle:
    tanimoto_freqs = pickle.load(handle)

def plot_scores(x, n_bins=50, cutoff = None):
  fig, axs = plt.subplots(1, 2, tight_layout=True, figsize=(10, 6))

  N, bins, patches = axs[0].hist(x, bins=n_bins)
  axs[0].set_xlabel("Scores")
  axs[0].set_ylabel('Frequency')



  fracs = N / N.max()
  norm = colors.Normalize(fracs.min(), fracs.max())

  for thisfrac, thispatch in zip(fracs, patches):
      color = plt.cm.viridis(norm(thisfrac))
      thispatch.set_facecolor(color)
  
  axs[1].boxplot(x)
  axs[1].set_ylabel("Scores")
  axs[1].set_xlabel('Boxplot')
  
  if cutoff:
    axs[0].axvline(x = cutoff, color = 'r', label = 'Score Cutoff')
    axs[1].axhline(y = cutoff, color = 'r', label = 'Score Cutoff')


plot_scores(tanimoto_freqs.keys(), tanimoto_freqs.values())

