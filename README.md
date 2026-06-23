# scRG-NMTF: Robust Graph Regularized Multi-view Non-negative Matrix Tri-factorization

[![MATLAB](https://img.shields.io/badge/MATLAB-R2021b%2B-blue.svg)](https://www.mathworks.com/products/matlab.html)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)
[![Paper](https://img.shields.io/badge/Paper-Pattern%20Recognition%20(Under%20Review)-orange.svg)](#)

> Official MATLAB implementation of the paper: **"Robust Graph Regularized Multi-view Non-negative Matrix Tri-factorization for Integrative Single-cell Multi-omics Clustering"**

---

## 📖 Introduction

The rapid evolution of single-cell multi-omics sequencing offers profound insights into cellular heterogeneity. However, computational integration is severely hindered by extreme dropout sparsity, inherent modality gaps, and the prohibitive cost of analyzing massive datasets.

**scRG-NMTF** is a robust, highly scalable matrix factorization framework designed for integrative single-cell multi-omics clustering. It effectively bridges semantic modality gaps and reconstructs faithful latent representations from highly incomplete data.

### ✨ Key Features
- **Zero-loss ReLU-mirrored expansion:** Perfectly preserves dense features from standard dimensionality reductions while strictly satisfying NMF's non-negativity constraint.
- **Projective Tri-factorization:** Decouples modality-specific features from topological structures, effectively bridging severe modality gaps.
- **Feature-Level Dynamic Sub-manifold Regularization:** Reconstructs the Laplacian graph strictly on the low-dimensional feature subspace, drastically reducing the dynamic topological complexity from $\mathcal{O}(N^2)$ to $\mathcal{O}(K^2)$.
- **Strict Linear Scalability:** Achieves $\mathcal{O}(N)$ computational complexity, making it exceptionally scalable for massive single-cell atlases (100k+ cells).

---

## 🚀 Framework Overview

![scRG-NMTF Framework](https://via.placeholder.com/800x400.png?text=Please+Upload+Your+Framework+Image+Here)
*(Please replace the placeholder image link above with the actual link to your `Framework.png` after uploading it to this repository).*

---

## 📁 Repository Structure

```text
scRG-NMTF/
│
├── data/                       # Contains the lightweight demo dataset
│   └── PBMC-10k_GoldStandard.mat  # Pre-processed features (RNA + ATAC)
│
├── src/                        # Core optimization engine
│   └── scrgnmf_engine.m        # The Multiplicative Update Rules (MUR) implementation
│
├── utils/                      # Evaluation metrics and graph construction tools
│   ├── get_acc.m               # Accuracy (ACC)
│   ├── get_nmi.m               # Normalized Mutual Information (NMI)
│   ├── get_ari.m               # Adjusted Rand Index (ARI)
│   ├── get_ami.m               # Adjusted Mutual Information (AMI)
│   └── my_construct_knn_local.m # Fast KNN graph construction
│
├── demo_PBMC10k.m              # Main script to execute the pipeline on the demo data
├── plot_convergence.m          # Script to visualize the monotonic convergence curve
├── LICENSE                     # MIT License
└── README.md                   # This document

```

---

## ⚙️ Prerequisites

* **MATLAB:** R2021b or later is strongly recommended.
* **Toolboxes required:**
* Statistics and Machine Learning Toolbox (required for `pca`, `kmeans`, `knnsearch`, and `pdist2` functions).



---

## ⚡ Quick Start

We provide a lightweight, pre-processed version of the **PBMC-10k** dataset to help you quickly reproduce our algorithm and evaluate its performance.

**Step 1: Clone the repository**

```bash
git clone [https://github.com/](https://github.com/)[ys15802550432-wq]/scRG-NMTF.git
cd scRG-NMTF

```

**Step 2: Run the clustering algorithm**
Open MATLAB, navigate to the `scRG-NMTF` folder, and run the main demo script in the Command Window:

```matlab
>> demo_PBMC10k

```

This script will automatically load the dataset, execute the scRG-NMTF optimization engine for 5 independent runs, and output the clustering evaluation metrics (ACC, NMI, ARI, AMI).

**Step 3: Plot the convergence curve**
After the demo finishes successfully, run the visualization script to verify the strictly monotonic convergence behavior of our algorithm:

```matlab
>> plot_convergence

```

---

## 📊 Data Availability

To facilitate a quick evaluation of our algorithm, we provide a pre-processed demo dataset (`data/PBMC-10k_GoldStandard.mat`) containing the dimensionality-reduced features of the PBMC 10k dataset.

For those who wish to reproduce the pipeline from scratch on all 8 benchmark datasets, the raw single-cell count matrices can be openly accessed from their original sources:

* **10x Genomics Datasets** (PBMC 10k, PBMC 3k, Human Brain 3k, Lymphoma 14k, Mouse Brain 5k): [10x Genomics Datasets Portal](https://www.10xgenomics.com/resources/datasets)
* **SHARE-seq Mouse Skin:** [NCBI GEO (GSE140203)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE140203)
* **CITE-seq Datasets** (BMMNC 10k, PBMC 1.4k): [NCBI GEO (GSE100866)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE100866)

---

## 📝 Citation

If you find this code, the dataset, or our framework useful for your research, please consider citing our paper:

```bibtex
@article{yuan2024scrgnmtf,
  title={Robust Graph Regularized Multi-view Non-negative Matrix Tri-factorization for Integrative Single-cell Multi-omics Clustering},
  author={Yuan, Sheng and Che, Hangjun and Su, Jingfeng and Wang, Huiwei and Wang, Yadi and Liu, Cheng and Leung, Man-Fai},
  journal={Pattern Recognition},
  year={2026},
  note={Under Review}
}

```

---

## ✉️ Contact

For any questions, issues, or collaborations, please feel free to open an issue or contact:

* **Sheng Yuan** - [ys15802550432@email.swu.edu.cn]
* **Hangjun Che** - [hjche123@swu.edu.cn]
