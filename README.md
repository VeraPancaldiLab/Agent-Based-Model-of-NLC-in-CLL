# An Agent-Based Model of Monocyte Differentiation into Tumor-Associated Macrophages in Chronic Lymphocytic Leukemia

## Summary

This repository contains data related to the publication [An Agent-Based Model of Monocyte Differentiation into Tumor-Associated Macrophages in Chronic Lymphocytic Leukemia](https://www.biorxiv.org/content/10.1101/2021.12.17.473137v3) by Verstraete et al. 

It presents an agent-based model (ABM) describing the differentiation dynamics of monocytes into tumour-associated macrophages upon contact with cancer B cells in the context of chronic lymphocytic leukemia. The goal of this ABM is to reproduce dynamics compatible with an *in vitro* experimental setting of monocytes differentiation in presence of leukemic (B-CLL).

The model can be run online at [https://www.netlogoweb.org/](https://www.netlogoweb.org/). To do this: 
- download the NetLogo model [here](https://github.com/VeraPancaldiLab/Agent-Based-Model-of-NLC-in-CLL/blob/main/netlogo/ABM_NLC_CLL.nlogo) (`right click` the `Raw` button at the top of the file, select `Save Link As…`, choose the location on your computer where you want to save the file, and select `Save`).
- go at [https://www.netlogoweb.org/](https://www.netlogoweb.org/) and click on NetLogo Web
- upload the model
- click `setup` and `go`.

## Project organisation

- **netlogo** source file for the agent-based model.
- **openmole** source files for exploring the model on the OpenMOLE platform https://openmole.org
- **figures** exploration outputs and source files used to generate the figures and analysis of the related publication

![fig1](https://raw.githubusercontent.com/VeraPancaldiLab/Agent-Based-Model-of-NLC-in-CLL/master/images/Figure1-edit2-01.png)
_**Fig.1 Experimental setups and datasets from in vitro PBMC cultures from CLL patients.** A) Experimental set-up. B) Visualization of NLCs at 10 days of in vitro culture from two different patients. C) Time course datasets produced from the PBMC autologous cultures from 9 patients. D) Heterologous co-cultures._

![fig2](https://raw.githubusercontent.com/VeraPancaldiLab/Agent-Based-Model-of-NLC-in-CLL/master/images/Figure2-edit2.png)
_**Fig. 2. ABM representations.** A) Netlogo simulation of 5000 cells. B) Schematic diagram of the agents’ states and behaviors._

## Videos

- A video of a simulation of 5000 cells can be found at https://doi.org/10.6084/m9.figshare.21592368.v1.
- A video of a live-imaging microscopy capture attempting to acquire cell pattern formations is available at https://doi.org/10.6084/m9.figshare.22110926.v1. 
