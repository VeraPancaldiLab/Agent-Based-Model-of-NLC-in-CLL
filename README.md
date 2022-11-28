# An Agent-Based Model of Monocyte Differentiation into Tumor-Associated Macrophages in Chronic Lymphocytic Leukemia

## Summary

This repository contains data related to the publication [An Agent-Based Model of Monocyte Differentiation into Tumor-Associated Macrophages in Chronic Lymphocytic Leukemia](https://www.biorxiv.org/content/10.1101/2021.12.17.473137v3) by Verstraete et al. 

It presents an agent-based model (ABM) describing the differentiation dynamics of monocytes into tumour-associated macrophages upon contact with cancer B cells in the context of chronic lymphocytic leukemia. The goal of this ABM is to reproduce dynamics compatible with an *in vitro* experimental setting of monocytes differentiation in presence of leukemic (B-CLL).

<img src="https://github.com/VeraPancaldiLab/Agent-Based-Model-of-NLC-in-CLL/tree/main/images/Figure1-edit2-01.png" width="50%" height="50%">
_**Fig.1 Experimental setups and datasets from in vitro PBMC cultures from CLL patients.** A) Experimental set-up. B) Visualization of NLCs at 10 days of in vitro culture from two different patients. C) Time course datasets produced from the PBMC autologous cultures from 9 patients. D) Heterologous co-cultures._

<img src="[https://user-images.githubusercontent.com/16319829/81180309-2b51f000-8fee-11ea-8a78-ddfe8c3412a7.png]([https://github.com/VeraPancaldiLab/Agent-Based-Model-of-NLC-in-CLL/tree/main/images/Figure1-edit2-01.png](https://github.com/VeraPancaldiLab/Agent-Based-Model-of-NLC-in-CLL/tree/main/images/Figure2-edit2.png))" width=50% height=50%>
_**Fig. 2. ABM representations.** A) Netlogo simulation of 5000 cells. B) Schematic diagram of the agents’ states and behaviors._


The model can be run online at [https://www.netlogoweb.org/](https://www.netlogoweb.org/). To do this: 
- download the NetLogo model [here](https://github.com/VeraPancaldiLab/Agent-Based-Model-of-NLC-in-CLL/blob/main/netlogo/ABM_NLC_CLL.nlogo) (`right click` the `Raw` button at the top of the file, select `Save Link As…`, choose the location on your computer where you want to save the file, and select `Save`).
- go at [https://www.netlogoweb.org/](https://www.netlogoweb.org/) and click on NetLogo Web
- upload the model
- click `setup` and `go`.

## Project organisation

- **netlogo** source file for the agent-based model.
- **openmole** source files for exploring the model on the OpenMOLE platform https://openmole.org
- **figures** exploration outputs and source files used to generate the figures and analysis of the related publication
