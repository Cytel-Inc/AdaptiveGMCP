# AdaptGMCP

## Installation

You can install the AdaptGMCP package in R directly from the GitHub repository using the `remotes` package:

```R
remotes::install_github("Cytel-Inc/AdaptiveGMCP", ref = "master")
```

This will download and install the latest version of the package from the master branch.

## Important Functions

The following are some of the key functions provided by the AdaptGMCP package:

- `simMAMSMEP`: Runs simulations for adaptive group sequential multiple testing procedures.
- `adaptGMCP_CER`: Performs analysis using the CER (Conditional Error Rate) method.
- `adaptGMCP_PC`: Performs analysis using the p-value combination method.

For detailed usage and arguments of each function, please refer to the help section in R after installing the package (e.g., `?simMAMSMEP`).

## Quick Start: Running `simMAMSMEP`

To help you get started, an example script demonstrating the use of the `simMAMSMEP` function is provided in the repository:

[internalData/MAMSMEP_Simulation_Example.R](internalData/MAMSMEP_Simulation_Example.R)

You can open and run this script in R to see a typical workflow.

For detailed documentation on all arguments, run:

```R
?simMAMSMEP
```

## Quick Start: Running `adaptGMCP_CER`

An example script for the `adaptGMCP_CER` function is available:

[internalData/AdaptGMCP_CER_Analysis_Example.R](internalData/AdaptGMCP_CER_Analysis_Example.R)

To use this function, open and run the script in R. This will demonstrate how to perform analysis using the Conditional Error Rate method.

For detailed documentation on all arguments, run:

```R
?adaptGMCP_CER
```

## Quick Start: Running `adaptGMCP_PC`

An example script for the `adaptGMCP_PC` function is available:

[internalData/AdaptGMCP_Analysis_Example.R](internalData/AdaptGMCP_Analysis_Example.R)

Open and run this script in R to see how to use the p-value combination method for your analysis.

For detailed documentation on all arguments, run:

```R
?adaptGMCP_PC
```
