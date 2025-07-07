# Implicit SINDy- Handling Rational Nonlinearity

This repository contains all the code, notebooks, and symbolic tools for performing Implicit SINDy (SINDy-PI).
---

## 📁 Folder Structure


</pre><code> 
Chapter-6-Implicit-SINDy/
├── Examples_Related/ # Output files (LaTeX, CSV) for discovered models
├── Example_1_MMKinetics.ipynb # Michaelis–Menten example 
├── Example_2_Microbial.ipynb # Microbial regulatory model 
├── implicit_to_explicit.py # Converts implicit models to explicit symbolic rational form
├── sindypi_functions.py # SINDy-PI pipeline, thresholding, evaluation
├── rescaling_utils.py # Tools for rescaling explicit rational form
├── steady_states.py # Steady-state computation
├── plot_utils.py # Basic plotting utilities
├── network_utils.py # Interaction network visualization using Jacobians
├── environment.yml # Reproducible Python environment
└── README.md # This file
 <code></pre>
---


## 🧪 Environment Setup

To reproduce results or run notebooks, install [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or [Anaconda](https://www.anaconda.com/).

Then execute in terminal or Git CMD:

```bash
conda env create -f environment.yml
conda activate implicit_sindy
```
This will install all required dependencies into a new environment.

---

🚀 Running the Code
Launch Jupyter after activating the environment:

```bash
jupyter notebook
```

Then open and run any of the notebooks.

---


## 📦 Dependencies
This project uses Python 3.10 with the following key packages:

- pysindy
- sympy
- scikit-learn, scipy, numpy, pandas
- matplotlib, seaborn, plotly, scikit-learn
- graphviz, tqdm

All dependencies are listed in `environment.yml`.

---

## 📬 Contact

For questions or suggestions, please open an issue on the GitHub repository:
👉 [Alka-CBhub/Chapter-6-Implicit-SINDy](https://github.com/Alka-CBhub/Chapter-6-Implicit-SINDy)

---