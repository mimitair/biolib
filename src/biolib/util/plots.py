import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3
from pathlib import Path

def plotVenn2(s1: set, s2: set, plot_title: str, label1: str, label2: str, out_file: Path = None) -> None:
    plt.figure(figsize=(10,10))
    plt.title(plot_title)

    venn2((s1, s2), set_labels = (label1, label2))

    plt.show()
    
    if out_file:
        plt.savefig(out_file, format='svg')

    return None

def plotVenn3(s1: set, s2: set, s3: set, plot_title: str, label1: str, label2: str, label3: str, out_file: Path = None) -> None:
    plt.figure(figsize=(10,10))
    plt.title(plot_title)

    venn3((s1, s2, s3), set_labels=(label1, label2, label3))

    plt.show()

    return None
    
    
