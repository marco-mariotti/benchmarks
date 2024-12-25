

rule graph_benchmarks:
    input:
        "{RESULTS_DIR}/collected_results.csv"
    output:
        "{RESULTS_DIR}/{benchmark}.pdf"
    run:
        import pandas as pd

        benchmark = wildcards.benchmark  # MaxRSS_kB or Elapsed
        f = input[0]

        # Load the uploaded file to inspect the data
        data = pd.read_csv(f)

        unique_libraries = data['Library'].nunique()
        palette = sns.color_palette("tab20", unique_libraries)  # Use a larger palette

        # Create a facet plot based on Genome and MaxLength
        g = sns.relplot(
            data=data,
            x='NumberRows',
            y=benchmark,
            hue='Library',
            col='Genome',
            row='MaxLength',
            kind='line',
            marker='o',
            height=4,
            aspect=1.5,
            palette=palette,
        )

        # Set log scale for the x-axis
        for ax in g.axes.flat:
            ax.set_xscale('log')

        # Adjust the layout and add titles
        g.set_axis_labels("Number of Rows (log scale)", benchmark)
        g.set_titles("Genome: {col_name} | MaxLength: {row_name}")
        g.tight_layout()

        # Show the plot
        plt.show()




