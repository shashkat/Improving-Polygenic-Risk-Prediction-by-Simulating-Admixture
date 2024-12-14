- admixture_ld_plot_1.png was made using following settings:
    - num_bins = 100
    - bin_length = 10000
    - sns.set_theme(rc={'figure.figsize': (15,10)}, style = 'white')
    - sns.lineplot(df, x = 'distance', y = 'correlation_values', hue = 'label', errorbar= 'ci')
    - sh

- admixture_ld_plot_2.png was made using following settings:
    - same everything as above, just modified plot's y-axis to be log scaled.

- admixture_ld_plot_3.png was made using following settings:
    - num_bins = 1000
    - bin_length = 10000
    - also since it was tzking a lot of time when finding all the bin pairs with a certian gap, I only found a tenth of the initial number of bin pairs for a given gap.

- admixture_ld_plot_4.png was made using following settings:
    - num_bins = 150
    - bin_length = 10000
    - undid the part "tenth of the initial number of bin pairs" which was done in prev variant.
    - though the num_bins was 150, i plotted only for bin-distances till 100 in the image, as in near-end values of bin-distance, the data gets less and hence plot has more variance. 
    - also the xlabel was slightly inaccurate till now (bin_indices were called as distance). Corrected it.
    - log scale of y axis

- admixture_ld_plot_5.png was made using following settings:
    - same everything as above, just no log scale of y axis.