# iLipidome package

- [Overview](#overview)
- [How to install](#how-to-install)
- [System Requirements](#system-requirements)
- [Installation Guide](#installation-guide)
- [Quick Example](#quick-example)
- [License](#license)

# Overview
This tutorial presents a series of ``iLipidome`` functions that facilitate a comprehensive comparison of lipid profiles using a novel substructure-based approach. iLipidome is an innovative method that leverages the lipid biosynthetic network to analyze lipidomics data, taking into account the interdependence and interconnectedness of measured lipids. It provides 'Lipid Substructure Analysis' functionality, allowing users to decompose lipids into substructures, convert lipid expression into substructure expression, reconstruct the lipid biosynthetic network, and identify significant altered lipid pathways and their genetic origins. iLipidome currently supports "two-group comparison" and performs substructure analysis based on fatty acids, lipid species, or lipid classes, enabling comprehensive comparisons of lipid profiles across different levels. Our goal is to empower researchers with a deeper understanding of the intricate changes in lipidomics observed across various samples.


# How to install

1. Make sure you have the `devtools` R package. If you do not already have it installed, install it using `install.packages("devtools")`.
2. Run `devtools::install_github("LewisLabUCSD/iLipidome-package")` in R to install the iLipidome package, and you're done!

# System Requirements
## Hardware requirements
To run example datasets with iLipidome, you only need a standard computer with sufficient RAM and R software version 4.0.0 or higher installed.

## Software requirements
### OS Requirements
The functions and example datasets have been tested on the following systems:
+ macOS: Ventura (13.0)

### R Dependencies
The version information about R, the OS and attached or loaded packages for `iLipidome` are listed below.

![image](readme_fig/try.png)

## Upload lipidomics data
iLipidome only requires users to upload a processed lipid expression table (data.frame) for analysis. The table should have lipids as rows and samples as columns. Lipid names should be placed in the first column, labeled as “feature”, and sample names should be in the first row. It is important to have a minimum of two samples in each group for accurate statistical calculations. Depending on the data source, preprocessing and normalization techniques like missing value imputation or log transformation may be necessary to improve analysis outcomes.

Lipid names in the table can be represented in two formats:
1. When the exact identity of FAs is unknown, the lipids can be represented using the following format:
[LipidClassAbbreviation]_[sum of FA chain length] : [sum of FA double bonds] ; [sum of FA oxygens]
For example, PC_34:1;0 or TAG_52:1;0
2. When the exact identity of FAs is known, the lipids can be represented using the following format:
[LipidClassAbbreviation]_[FA1 chain length] : [FA1 double bonds] ; [FA1 oxygens]_[FA2 chain length] : [FA2 double bonds] ; [FA2 oxygens]…
For example, PC_16:0;0_18:1;0 or TAG_16:0;0_18:0;0_18:1;0

You can refer to the ‘supported_lipid_class.csv’ file for the supported lipid classes, their abbreviations, and the corresponding number of FAs. Note that when using the exact identity format of FAs, we will verify if the fatty acid numbers match those recorded in the ‘supported_lipid_class.csv’ file. If they do not match, the analysis will be interrupted. Also, lipid classes with the same number of FAs (e.g., PC, PE) in the same pathways (e.g., Glycerophospholipid) should have a consistent lipid naming format. For example, PC_36:0;0 and PE_34:0;0 or PC_18:0;0_18:0;0 and PE_16:0;0_18:0;0. Additionally, dihydrosphingolipids (dh-) specify sphingolipids with sphingoid bases of 18:0:2 instead of 18:1:2.
