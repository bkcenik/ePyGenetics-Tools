# RetroToolsPy
## Plotting TE copy numbers and counts

Date created: September 25, 2023 2:29 PM
Last edited time: September 25, 2023 2:38 PM

This collection of python scripts is for taking a gtf file and a list of TEs and plotting either their size distribution median, or the size against the copy number.

`plotGeneCopiesVsMedianLength.py` calculates the number of copies of a retroelement and the median lengths for all copies, and makes a scatterplot, color coding them by family name.

```python
import sys
import matplotlib.pyplot as plt
import numpy as np

# Check if the correct number of arguments are provided
if len(sys.argv) != 3:
    print("Usage: python plotGeneCopiesVsMedianLength.py <path to input gtf> <path to text file containing genes>")
    sys.exit(1)

# Retrieve the GTF file path and the path to the text file containing gene names from command-line arguments
file_path = sys.argv[1]
gene_list_path = sys.argv[2]

# Read the list of gene names from the provided file
with open(gene_list_path, 'r') as f:
    genes = [line.strip() for line in f]

# Initialize a dictionary to store the sizes of the instances for each gene and their family_id
sizes_dict = {gene: {"sizes": [], "family_id": set()} for gene in genes}

# Open the GTF file and read it line by line
with open(file_path, 'r') as file:
    for line in file:
        # Skip comment lines
        if line.startswith('#'):
            continue

        # Split the line into columns
        columns = line.strip().split('\t')

        # Check whether the attributes field contains any of the specified gene names
        attributes = columns[8]

        for gene in genes:
            if gene in attributes:
                # Calculate the size of the instance and append it to the corresponding list in the dictionary
                start = int(columns[3])
                end = int(columns[4])
                size = end - start + 1
                sizes_dict[gene]["sizes"].append(size)

                # Extract and store the family_id for the gene
                family_id = attributes.split('family_id "')[1].split('";')[0]
                sizes_dict[gene]["family_id"].add(family_id)

# Prepare a color map for different family_ids
unique_family_ids = set(family_id for gene_data in sizes_dict.values() for family_id in gene_data["family_id"])
color_map = {family_id: plt.cm.tab10(i / len(unique_family_ids)) for i, family_id in enumerate(unique_family_ids)}

# Prepare data for plotting
num_copies = []
median_lengths = []
colors = []
labels = []

for gene, gene_data in sizes_dict.items():
    sizes = gene_data["sizes"]
    family_ids = gene_data["family_id"]

    if sizes:
        num_copies.append(len(sizes))
        median_lengths.append(np.median(sizes))
        colors.append(color_map[next(iter(family_ids))])  # Use color of the first family_id found for this gene
        labels.append(gene)

# Check if there is any data to plot
if not num_copies:
    print("No instances of the specified genes were found in the GTF file.")
    sys.exit(1)

# Create a scatter plot of number of copies against median length, colored by family_id
plt.figure(figsize=(10, 6))
for x, y, color, label in zip(num_copies, median_lengths, colors, labels):
    plt.scatter(x, y, color=color, label=label)

# Annotate each point with the gene name
for i, label in enumerate(labels):
    plt.annotate(label, (num_copies[i], median_lengths[i]))

plt.xlabel('Number of Copies')
plt.ylabel('Median Length (bp)')
plt.title('Number of Copies vs Median Length of Gene Instances (Colored by Family ID)')

# Create a legend mapping color to family_id
legend_handles = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=color_map[family_id], markersize=10) 
                  for family_id in unique_family_ids]
plt.legend(legend_handles, unique_family_ids, loc='upper right')

# Save the plot as a PDF file
plt.savefig('gene_copies_vs_median_length_colored.pdf', format='pdf', bbox_inches='tight')

# Display the plot
plt.show()
```

Run as `python plotGeneCopiesVsMedianLength.py <path to input gtf> <path to TE list>`.


The text file should be something like:

```
MMERGLN-int
RLTR1B-int
LTRIS2
IAPEY3-int
BGLII_Mus
ERVB7_3-LTR_MM
GSAT_MM
L1Md_Gf
```

So, for instance, only include the base name of the element, no semicolons or family IDs. The script will pull the family ID itself from the gtf.

`calcTEsizeAll.py` will calculate the median length for a list of elements and make box and whisker plots. it will adjust the scale to exclude outliers in the graph while maintaining the usage of them in computing the stats.

```python
import sys
import matplotlib.pyplot as plt
import numpy as np

# Check if the correct number of arguments are provided
if len(sys.argv) != 3:
    print("Usage: python calcTEsizeAll.py <path to input gtf> <path to text file containing genes>")
    sys.exit(1)

# Retrieve the GTF file path and the path to the text file containing gene names from command-line arguments
gtf_file_path = sys.argv[1]
gene_list_path = sys.argv[2]

# Read the list of gene names from the provided file
with open(gene_list_path, 'r') as f:
    genes = [line.strip() for line in f]

# Initialize a dictionary to store the sizes and family_id of the instances for each gene
data_dict = {gene: {"sizes": [], "family_id": set()} for gene in genes}

# Open the GTF file and read it line by line
with open(gtf_file_path, 'r') as file:
    for line in file:
        # Skip comment lines
        if line.startswith('#'):
            continue

        # Split the line into columns
        columns = line.strip().split('\t')

        # Check whether the attributes field contains any of the specified gene names
        attributes = columns[8]

        for gene in genes:
            if f'gene_id "{gene}"' in attributes:
                # Calculate the size of the instance and append it to the corresponding list in the dictionary
                start = int(columns[3])
                end = int(columns[4])
                size = end - start + 1
                data_dict[gene]["sizes"].append(size)

                # Extract and store the family_id for the gene
                family_id = attributes.split('family_id "')[1].split('";')[0]
                data_dict[gene]["family_id"].add(family_id)

# Prepare a color map for different family_ids
unique_family_ids = set(family_id for gene_data in data_dict.values() for family_id in gene_data["family_id"])
color_map = {family_id: plt.cm.tab10(i / len(unique_family_ids)) for i, family_id in enumerate(unique_family_ids)}

# Prepare data for plotting
data_to_plot = []
colors = []
for gene, gene_data in data_dict.items():
    if gene_data["sizes"]:
        data_to_plot.append(gene_data["sizes"])
        colors.append(color_map[next(iter(gene_data["family_id"]))])  # Use color of the first family_id found for this gene

# Check if there is any data to plot
if not data_to_plot:
    print("No instances of the specified genes were found in the GTF file.")
    sys.exit(1)

# Create box plots, color by family_id, and exclude outliers in the plot
plt.figure(figsize=(10, 6))
bp = plt.boxplot(data_to_plot, patch_artist=True, showfliers=False)  # The showfliers parameter is set to False to exclude outliers
for patch, color in zip(bp['boxes'], colors):
    patch.set_facecolor(color)

# Create a legend mapping color to family_id
legend_handles = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=color_map[family_id], markersize=10)
                  for family_id in unique_family_ids]
plt.legend(legend_handles, unique_family_ids, loc='upper right')

plt.ylabel('Size (bp)')
plt.xticks(ticks=range(1, len(genes) + 1), labels=genes, rotation=45)
plt.title('Size Distribution of Gene Instances (Colored by Family ID)')

# Save the plot as a PDF file
plt.savefig('gene_sizes_box_plot_colored_no_outliers.pdf', format='pdf', bbox_inches='tight')

# Display the plot
plt.show()
```

Run as `python calcTEsizeAll.py <path to input gtf> <path to TE list>`

Et voila!
